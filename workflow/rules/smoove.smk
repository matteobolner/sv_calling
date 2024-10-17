rule smoove_call:
    input:
        bam=get_bam,
        bam_idx=get_bai,
        ref=config["ref_genome"],
        ref_idx=config["ref_genome_fai"],
        #dummy_header=config["dummy_header"],
    output:
        outdir=directory("data/sv_variants/smoove/called/{sample}/"),
        vcf=temp(
            "data/sv_variants/smoove/called/{sample}/{sample}-smoove.genotyped.vcf.gz"
        ),
        idx=temp(
            "data/sv_variants/smoove/called/{sample}/{sample}-smoove.genotyped.vcf.gz.csi"
        ),
    container:
        "docker://brentp/smoove:latest"
    threads: 4
    shell:
        # the code after || outputs the dummy vcf if no variants are found in the bam - this is necessary for further processing as empty files cause errors
        """smoove call --outdir {output.outdir} --name {wildcards.sample} --fasta {input.ref} -p {threads} --genotype {input.bam} """  #|| {{ sed 's/SS_69/{wildcards.sample}/g' {input.dummy_header} | bgzip -c > {output.vcf} && tabix -C -p vcf {output.vcf} ; }} """


rule smoove_merge:
    input:
        vcfs=expand(
            "data/sv_variants/smoove/called/{sample}/{sample}-smoove.genotyped.vcf.gz",
            sample=samples.index,
        ),
        ref=config["ref_genome"],
    output:
        dir=directory("data/sv_variants/smoove/merged/"),
        vcf="data/sv_variants/smoove/merged/all.sites.vcf.gz",
    container:
        "docker://brentp/smoove:latest"
    shell:
        "smoove merge --name all -f {input.ref} --outdir {output.dir} {input.vcfs}"


rule joint_genotype_svs:
    input:
        ref=config["ref_genome"],
        vcf="data/sv_variants/smoove/merged/all.sites.vcf.gz",
        bam=get_bam,
    output:
        "data/sv_variants/smoove/genotyped/{sample}-joint-smoove.genotyped.vcf.gz",
        "data/sv_variants/smoove/genotyped/{sample}-joint-smoove.genotyped.vcf.gz.csi",
    params:
        dir="data/sv_variants/smoove/genotyped/",
    container:
        "docker://brentp/smoove:latest"
    shell:
        "smoove genotype -d -x -p {threads} --name {wildcards.sample}-joint --outdir {params.dir} --fasta {input.ref} --vcf {input.vcf} {input.bam}"


rule merge_genotyped_samples:
    input:
        vcfs=expand(
            "data/sv_variants/smoove/genotyped/{sample}-joint-smoove.genotyped.vcf.gz",
            sample=samples.index,
        ),
        indexes=expand(
            "data/sv_variants/smoove/genotyped/{sample}-joint-smoove.genotyped.vcf.gz.csi",
            sample=samples.index,
        ),
    params:
        dir="data/sv_variants/smoove/genotyped/all_samples_genotyped",
    container:
        "docker://brentp/smoove:latest"
    output:
        file="data/sv_variants/smoove/genotyped/all_samples_genotyped.smoove.square.vcf.gz",
    shell:
        "smoove paste --name {params.dir} {input.vcfs}"


rule index_merged_genotyped_samples:
    input:
        "data/sv_variants/smoove/genotyped/all_samples_genotyped.smoove.square.vcf.gz",
    output:
        "data/sv_variants/smoove/genotyped/all_samples_genotyped.smoove.square.vcf.gz.tbi",
    params:
        "-p vcf",
    wrapper:
        "v4.5.0/bio/tabix/index"


rule annotate_with_gff_smoove:
    input:
        vcf="data/sv_variants/smoove/genotyped/all_samples_genotyped.smoove.square.vcf.gz",
        gff=config["ref_genome_gff"],
    output:
        vcf="data/sv_variants/smoove/annotated/annotated.vcf",
    container:
        "docker://brentp/smoove:latest"
    shell:
        "smoove annotate --gff {input.gff} {input.vcf} > {output.vcf}"

rule annotate_gene_variants_smoove_vep_cache:
    input:
        calls="data/sv_variants/smoove/annotated/annotated.vcf",
        plugins=config["vep_plugins_dir"],
        fasta=config["ref_genome"],
        fai=config["ref_genome_fai"],  # fasta index
        gff=config["ref_genome_gff"],
        #csi=config["ref_genome_gff_csi"],  # tabix index
    output:
        calls=temp("data/sv_variants/smoove/vep/ensembl_annotated.vcf"),
        stats="reports/vep_annotation.html",
    params:
        plugins=[],
        extra="--everything",
    threads: 4
    wrapper:
        "v4.5.0/bio/vep/annotate"

rule annotate_gene_variants_smoove_ncbi_gff:
    input:
        calls="data/sv_variants/smoove/annotated/annotated.vcf",
        cache=config["vep_cache_dir"],  # can be omitted if fasta and gff are specified
        plugins=config["vep_plugins_dir"],
    output:
        calls=temp("data/sv_variants/smoove/vep/ncbi_annotated.vcf"),
        stats="reports/vep_annotation.html",
    params:
        plugins=[],
        extra="--everything",
    threads: 4
    wrapper:
        "v4.5.0/bio/vep/annotate"

annotated_vcf = branch(
    config['vep_use_gff']
    then="data/sv_variants/smoove/vep/ncbi_annotated.vcf",
    otherwise="data/sv_variants/smoove/vep/ensembl_annotated.vcf"
)

rule bgzip_annotated_smoove:
    input:
        annotated_vcf,
    output:
        "data/sv_variants/smoove/vep/annotated.vcf.gz",
    params:
        extra="",  # optional
    threads: 1
    wrapper:
        "v4.5.0/bio/bgzip"


rule tabix_annotated_smoove:
    input:
        "data/sv_variants/smoove/vep/annotated.vcf.gz",
    output:
        "data/sv_variants/smoove/vep/annotated.vcf.gz.tbi",
    params:
        "-p vcf",
    wrapper:
        "v4.5.0/bio/tabix/index"

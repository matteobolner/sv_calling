rule smoove_call:
    input:
        bam="{bam_folder}/{sample}.bam",
        bam_idx="{bam_folder}/{sample}.bam.bai",
        ref=config["ref_genome"],
        ref_idx=config["ref_genome_fai"],
        dummy_header=config['dummy_header'],
    output:
        "{variants_folder}/smoove/called/{sample}/{sample}.disc.bam.orig.bam",
        temp("{variants_folder}/smoove/called/{sample}/{sample}.split.bam.orig.bam"),
        outdir=directory("{variants_folder}/smoove/called/{sample}"),
        vcf=temp("{variants_folder}/smoove/called/{sample}/{sample}-smoove.genotyped.vcf.gz"),
        idx=temp(
            "{variants_folder}/smoove/called/{sample}/{sample}-smoove.genotyped.vcf.gz.csi"
        ),
    threads: 1
    shell:
        """udocker -q run -v /home/bolner/work/:/home/bolner/work/ -w $PWD smoove smoove call --outdir {output.outdir} --name {wildcards.sample} --fasta {input.ref} -p {threads} --genotype {input.bam}  || {{ sed 's/SS_69/{wildcards.sample}/g' {input.dummy_header} | bgzip -c > {output.vcf} && tabix -C -p vcf {output.vcf} ; }} """


rule smoove_merge:
    input:
        vcfs=expand(
            "{variants_folder}/smoove/called/{sample}/{sample}-smoove.genotyped.vcf.gz",
            sample=samples.index,
        ),
        ref=config["ref_genome"],
    output:
        dir=directory("{variants_folder}/smoove/merged/"),
        vcf=temp("{variants_folder}/smoove/merged/all.sites.vcf.gz"),
    # conda:
    #    "../../../envs/smoove.yaml"
    shell:
        "udocker -q run -v /home/bolner/work/:/home/bolner/work/ -w $PWD  smoove smoove merge --name all -f {input.ref} --outdir {output.dir} {input.vcfs}"


rule joint_genotype_svs:
    input:
        ref=config["ref_genome"],
        vcf="{variants_folder}/smoove/merged/all.sites.vcf.gz",
        bam="{bam_folder}/{sample}.bam",
    output:
        temp("{variants_folder}/smoove/genotyped/{sample}-joint-smoove.genotyped.vcf.gz"),
        temp("{variants_folder}/smoove/genotyped/{sample}-joint-smoove.genotyped.vcf.gz.csi"),
    params:
        dir="{variants_folder}/smoove/genotyped/",
    # conda:
    #    "../../../envs/smoove.yaml"
    shell:
        "udocker -q run -v /home/bolner/work/:/home/bolner/work/ -w $PWD  smoove smoove genotype -d -x -p {threads} --name {wildcards.sample}-joint --outdir {params.dir} --fasta {input.ref} --vcf {input.vcf} {input.bam}"


rule merge_genotyped_samples:
    input:
        vcfs=expand(
            "{variants_folder}/smoove/genotyped/{sample}-joint-smoove.genotyped.vcf.gz",
            sample=samples.index,
        ),
        indexes=expand(
            "{variants_folder}/smoove/genotyped/{sample}-joint-smoove.genotyped.vcf.gz.csi",
            sample=samples.index,
        ),
    params:
        dir="{variants_folder}/smoove/genotyped/all_samples_genotyped",
    # conda:
    #    "../../../envs/smoove.yaml"
    output:
        file="{variants_folder}/smoove/genotyped/all_samples_genotyped.smoove.square.vcf.gz",
    shell:
        "udocker -q run -v /home/bolner/work/:/home/bolner/work/ -w $PWD  smoove smoove paste --name {params.dir} {input.vcfs}"


rule index_merged_genotyped_samples:
    input:
        vcf="{variants_folder}/smoove/genotyped/all_samples_genotyped.smoove.square.vcf.gz",
    output:
        index="{variants_folder}/smoove/genotyped/all_samples_genotyped.smoove.square.vcf.gz.tbi",
    shell:
        "tabix -f -p vcf {input.vcf}"


rule annotate_with_gff_smoove:
    input:
        vcf="{variants_folder}/smoove/genotyped/all_samples_genotyped.smoove.square.vcf.gz",
        gff=config["ref_genome_gff"],
    output:
        vcf="{variants_folder}/smoove/annotated/all_genes/annotated.vcf",
    shell:
        "udocker -q run -v /home/bolner/work/:/home/bolner/work/ -w $PWD  -v /home/bolner/work/resources/REFERENCE_GENOMES//sus_scrofa/GCF_000003025.6/:/home/bolner/work/resources/REFERENCE_GENOMES//sus_scrofa/GCF_000003025.6/ smoove smoove annotate --gff {input.gff} {input.vcf} > {output.vcf}"


# BND REMOVAL IS DONE BECAUSE VEP CANT PARSE THE ALT FIELD OF BND VARIANTS


rule remove_BND_before_VEP:
    input:
        vcf="{variants_folder}/smoove/annotated/all_genes/annotated.vcf",
    output:
        vcf="{variants_folder}/smoove/annotated/all_genes/annotated_no_BND.vcf",
    shell:
        "grep -v BND {input} > {output}"


rule annotate_gene_variants_smoove:
    input:
        calls="{variants_folder}/smoove/annotated/all_genes/annotated.vcf",  # .vcf, .vcf.gz or .bcf
        cache=config["vep_cache_dir"],  # can be omitted if fasta and gff are specified
        plugins=config["vep_plugins_dir"],
    output:
        calls=temp("{variants_folder}/smoove/vep/all_genes/annotated.vcf"),  # .vcf, .vcf.gz or .bcf
        stats="reports/variants/large/vep/all_genes.html",
    params:
        plugins=[],
        extra="--everything",
    log:
        "logs/vep/smoove/all_genes.log",
    threads: 4
    wrapper:
        "v3.13.6/bio/vep/annotate"


rule index_annotated_smoove:
    input:
        vcf="{variants_folder}/smoove/vep/all_genes/annotated.vcf",
    output:
        vcf="{variants_folder}/smoove/vep/all_genes/annotated.vcf.gz",
        idx="{variants_folder}/smoove/vep/all_genes/annotated.vcf.gz.tbi",
    shell:
        "bgzip -f {input.vcf} > {output.vcf} && tabix -p vcf {output.vcf}"


rule split_by_bed_smoove:
    input:
        vcf="{variants_folder}/smoove/vep/all_genes/annotated.vcf.gz",
        tbi="{variants_folder}/smoove/vep/all_genes/annotated.vcf.gz.tbi",
        bed="data/gene_info/{gene}/coordinates_padded.bed",
    output:
        vcf="{variants_folder}/smoove/vep/{gene}/annotated.vcf.gz",
        index="{variants_folder}/smoove/vep/{gene}/annotated.vcf.gz.tbi",
    shell:
        "tabix -R {input.bed} {input.vcf} -h | bgzip > {output.vcf} && tabix -p vcf {output.vcf}"

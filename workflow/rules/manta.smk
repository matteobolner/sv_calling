rule manta:
    input:
        ref=config['ref_genome'],
        samples=[expand("{bam_folder}/{{sample}}.bam", bam_folder=config['bam_folder'])],
        index=[expand("{bam_folder}/{{sample}}.bam.bai", bam_folder=config['bam_folder'])],
    output:
        vcf=temp(
            "{variants_folder}/manta/called/{sample}/{sample}.bcf"
        ),
        idx=temp(
            "{variants_folder}/manta/called/{sample}/{sample}.bcf.csi"
        ),

        cand_indel_vcf=expand("{variants_folder}/manta/called/{{sample}}/small_indels.vcf.gz", variants_folder=config['variants_folder']),
        cand_indel_idx=expand("{variants_folder}/manta/called/{sample}/small_indels.vcf.gz.tbi", variants_folder=config['variants_folder']),
        cand_sv_vcf=expand("{variants_folder}/manta/called/{sample}/cand_sv.vcf.gz", variants_folder=config['variants_folder']),
        cand_sv_idx=expand("{variants_folder}/manta/called/{sample}/cand_sv.vcf.gz.tbi", variants_folder=config['variants_folder']),
    params:
        extra_cfg="",  # optional
        extra_run="",  # optional
    log:
        "logs/manta/{sample}.log",
    threads: 2
    resources:
        mem_mb=4096,
    wrapper:
        "v4.5.0/bio/manta"

if config["vep_use_gff"] == False:

    rule get_vep_cache:
        output:
            directory(config["vep_cache_dir"]),
        params:
            species=config["ensembl_ref"]["species"],
            build=config["ensembl_ref"]["build"],
            release=config["ensembl_ref"]["release"],
        log:
            "logs/vep/cache.log",
        cache: "omit-software"  # save space and time with between workflow caching (see docs)
        wrapper:
            "v4.5.0/bio/vep/cache"


rule setup_ref_genome:
    output:
        fasta=config["ref_genome"],
        faidx=config["ref_genome_fai"],
    shell:
        "genomers --accession {config[genome_accession]} --name {config[genome_name]} -g | gunzip | bgzip -c > {output.fasta} && samtools faidx {output.fasta} -o {output.faidx}"


rule setup_ref_genome_gff:
    output:
        gff=config["ref_genome_gff"],
        csi=config["ref_genome_gff_csi"],
    shell:
        "genomers --accession {config[genome_accession]} --name {config[genome_name]} -gff | gunzip | grep -v "  #" | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > {output.gff} && tabix -p gff --csi {output.gff}"

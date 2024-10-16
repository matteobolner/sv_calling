if config['vep_use_gff']:
    rule get_vep_cache:
        output:
            directory(config['vep_cache_dir']),
        params:
            species=config['ensembl_ref']['species'],
            build=config['ensembl_ref']['build'],
            release=config['ensembl_ref']['release'],
        log:
            "logs/vep/cache.log",
        cache: "omit-software"  # save space and time with between workflow caching (see docs)
        wrapper:
            "v4.5.0/bio/vep/cache"

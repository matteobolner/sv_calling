rule mosdepth:
    input:
        bam=get_bam,
        bai=get_bai,
    output:
        "stats/mosdepth/{sample}.mosdepth.global.dist.txt",
        #"stats/mosdepth/{sample}.per-base.bed.gz",  # produced unless --no-per-base specified
        summary="stats/mosdepth/{sample}.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/mosdepth/{sample}.log",
    params:
        extra="--fast-mode --no-per-base",  # optional
    # additional decompression threads through `--threads`
    threads: 4  # This value - 1 will be sent to `--threads`
    wrapper:
        "v4.5.0/bio/mosdepth"

rule multiqc_mosdepth:
    input:
        expand("stats/mosdepth/{sample}.mosdepth.global.dist.txt", sample=samples.index),
    output:
        "reports/multiqc/mosdepth/report.html",
        directory("reports/multiqc/mosdepth/report_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc/mosdepth.log",
    wrapper:
        "v4.5.0/bio/multiqc"

rule mosdepth:
    input:
        bam=get_bam,
        bai=get_bai,
    output:
        "stats/mosdepth/{sample}.mosdepth.global.dist.txt",
        summary="stats/mosdepth/{sample}.mosdepth.summary.txt",
    log:
        "logs/mosdepth/{sample}.log",
    params:
        extra="--fast-mode --no-per-base",
    threads: 4
    wrapper:
        "v4.5.0/bio/mosdepth"


rule multiqc_mosdepth:
    input:
        expand("stats/mosdepth/{sample}.mosdepth.global.dist.txt", sample=samples.index),
    output:
        "reports/multiqc/mosdepth/report.html",
        directory("reports/multiqc/mosdepth/report_data"),
    params:
        extra="--verbose",
    log:
        "logs/multiqc/mosdepth.log",
    wrapper:
        "v4.5.0/bio/multiqc"

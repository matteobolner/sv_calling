rule read_depth_annotations:
    input:
        samples="config/samples.tsv",
        vars="data/sv_variants/smoove/vep/annotated.vcf.gz"
    output:
        merged_table="data/sv_variants/smoove/vep/read_depth_annotations.tsv"
    script:
        "../scripts/read_depth_annotations.py"


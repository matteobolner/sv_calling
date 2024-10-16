import pandas as pd


configfile: "config/config.yaml"


samples = pd.read_table(config["samples"]).set_index("sample", drop=False)


wildcard_constraints:
    sample="[^/]+",


def get_bam(wildcards):
    sample = samples.loc[wildcards.sample]
    return sample.bam_path


def get_bai(wildcards):
    sample = samples.loc[wildcards.sample]
    return sample.bam_path + ".bai"

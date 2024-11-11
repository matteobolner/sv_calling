import pandas as pd
import numpy as np
from vcforge import VCFClass

vars = VCFClass(
    vcf_path=snakemake.input.vars, sample_info=snakemake.input.samples, add_info=True
)
vars.get_var_stats()

all_formats = vars.variants["FORMAT"].str.split(":").explode().unique()

format_info = {}

for i in all_formats:
    format_info[i] = vars.vcf.get_header_type(i)


format_info = pd.DataFrame(format_info).transpose()


dp = vars.format("DP", allele=0)
ro = vars.format("RO", allele=0)
ao = vars.format("AO", allele=0)

dp.index = vars.var_ids
ro.index = vars.var_ids
ao.index = vars.var_ids

ro_by_dp = ro / dp
ro.columns = [i + "|REF" for i in ro.columns]
ao.columns = [i + "|ALT" for i in ao.columns]

ro_by_dp.columns = [i + "|ref_by_DP" for i in ro_by_dp.columns]

vars.variants = (
    vars.variants.merge(ro, left_index=True, right_index=True)
    .merge(ao, left_index=True, right_index=True)
    .merge(ro_by_dp, left_index=True, right_index=True)
)
vars.variants.index.name = "ID"
vars.variants.to_csv(snakemake.output.merged_table, index=True, sep="\t")


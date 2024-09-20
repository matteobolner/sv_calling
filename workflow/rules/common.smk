import pandas as pd
samples = pd.read_table(config["samples"]).set_index("sample", drop=False)

wildcard_constraints:
    sample="[^/]+",

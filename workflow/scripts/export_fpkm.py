import pandas as pd

fpkm_df = pd.read_csv(snakemake.input.fpkm, sep = "\t")
fpkm_df.to_excel(snakemake.output.xlsx)
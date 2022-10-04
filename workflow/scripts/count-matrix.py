import pandas as pd


def get_column(strandedness):
    if pd.isnull(strandedness) or strandedness == "none":
        return 8  # non stranded protocol
    elif strandedness == "yes":
        return 9  # 3rd column
    elif strandedness == "reverse":
        return 10  # 4th column, usually for Illumina truseq
    else:
        raise ValueError(
            (
                "'strandedness' column should be empty or have the "
                "value 'none', 'yes' or 'reverse', instead has the "
                "value {}"
            ).format(repr(strandedness))
        )


def get_fpkm_column(strandedness):
    if pd.isnull(strandedness) or strandedness == "none":
        return 11  # non stranded protocol
    elif strandedness == "yes":
        return 12  # 3rd column
    elif strandedness == "reverse":
        return 13  # 4th column, usually for Illumina truseq
    else:
        raise ValueError(
            (
                "'strandedness' column should be empty or have the "
                "value 'none', 'yes' or 'reverse', instead has the "
                "value {}"
            ).format(repr(strandedness))
        )


counts = [
    pd.read_table(
        f, index_col=0, usecols=[3, get_column(strandedness)], header=None, skiprows=1
    )
    for f, strandedness in zip(snakemake.input, snakemake.params.strand)
]

fpkm = [
    pd.read_table(
        f,
        index_col=0,
        usecols=[3, get_fpkm_column(strandedness)],
        header=None,
        skiprows=1,
    )
    for f, strandedness in zip(snakemake.input, snakemake.params.strand)
]

gene_names = pd.read_table(snakemake.input[0], index_col=0, usecols=[3, 6], header=0)

for t, sample in zip(counts, snakemake.params.samples):
    t.columns = [sample]

for t, sample in zip(fpkm, snakemake.params.samples):
    t.columns = [sample]

matrix = pd.concat(counts, axis=1)
fpkm_mat = pd.concat(fpkm, axis=1)


matrix.index.name = "gene"
fpkm_mat.index.name = "gene"
gene_names.index.name = "gene"
#  collapse technical replicates
matrix = matrix.groupby(matrix.columns, axis=1).sum()
fpkm_mat = fpkm_mat.groupby(fpkm_mat.columns, axis=1).sum()

fpkm_mat["gname"] = gene_names["name"]

matrix.to_csv(snakemake.output[0], sep="\t")
fpkm_mat.to_csv(snakemake.output[1], sep="\t")

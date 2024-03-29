% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/carnival_helpers.R
\name{rename_count_rownames}
\alias{rename_count_rownames}
\title{rename ENSG IDS from OTP using fpkm table}
\usage{
rename_count_rownames(
  count_table,
  fpkm_table,
  fpkm_col = "gene",
  tibble_col = NULL
)
}
\arguments{
\item{count_table}{Table with ENSG IDs with version as rownames}

\item{fpkm_table}{FPKM table from rna-star-deseq pipeline `fpkm/all.tsv` as tibble}

\item{fpkm_col}{column name of genes in fpkm table}

\item{tibble_col}{If count_table is tibble give tibble col}
}
\description{
rename ENSG IDS from OTP using fpkm table
}
\examples{
NULL
}

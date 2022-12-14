% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSEA_functions.R
\name{run_gsea}
\alias{run_gsea}
\title{Run KEGG gene set enrichment on DE_tb table of genes}
\usage{
run_gsea(DE_tb, input_type, p_valcut = 0.05, species)
}
\arguments{
\item{DE_tb}{Data.frame like object of col 1 gene symbols, col2 LFC values}

\item{input_type}{String denoting the type of input given based on \link[org.Mm.eg.db]{org.Mm.eg.db}}

\item{p_valcut}{Cutoff passed to gseKEGG}

\item{species}{Name of species analyzed in workflow}
}
\value{
\link[clusterProfiler]{gseKEGG}
}
\description{
Run KEGG gene set enrichment on DE_tb table of genes
}

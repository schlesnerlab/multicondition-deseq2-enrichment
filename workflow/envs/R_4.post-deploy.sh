#!/bin/bash
Rscript -e 'BiocManager::install("dorothea", upgrade = "never")'
Rscript -e 'BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Mmusculus.UCSC.mm10", "DO.db", "org.Mm.eg.db", "org.Hs.eg.db","ReactomePA", "ensembldb"))'
#Rscript -e 'devtools::install_github("saezlab/decoupleR", ref = "56dc1a3")'
Rscript -e 'devtools::install_github("saezlab/OmnipathR", ref = "c5f63b4")'
Rscript -e 'devtools::install_github("Christian-Heyer/mitch")'
Rscript -e 'devtools::install("workflow/scripts/RNAscripts", upgrade = "never", dependencies = F)'

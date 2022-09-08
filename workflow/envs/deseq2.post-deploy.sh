#!/bin/sh


#Rscript -e 'BiocManager::install(c("DO.db", "reactome.db", "ReactomePA", "org.Mm.eg.db", "org.Hs.eg.db", "GenomeInfoDb"))'

Rscript -e 'devtools::install("workflow/scripts/RNAscripts", upgrade = "never")'

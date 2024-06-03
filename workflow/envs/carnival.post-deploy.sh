#!/bin/sh
Rscript -e 'devtools::install_github("grimbough/biomaRt", upgrade = "never", dependencies = FALSE)'
Rscript -e 'BiocManager::install("dorothea", upgrade = "never")'
Rscript -e 'devtools::install_github("saezlab/CARNIVAL")' 
Rscript -e 'remotes::install_github("saezlab/omnipathr")'
Rscript -e 'remotes::install_github("saezlab/decoupleR")'
#, ref = "8b3f0dff")'
#Rscript -e 'BiocManager::install(c("DO.db", "GO.db", "org.Mm.eg.db", "org.Hs.eg.db"), upgrade = "never")'
Rscript -e 'devtools::install("workflow/scripts/RNAscripts", upgrade = "never", dependencies = FALSE)'
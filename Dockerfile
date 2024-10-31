FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="b5b91967977d9d3f6a9dbaa8c01e3f99ae6aff2e6f7d0553360c16ca731282d0"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/R_4.yaml
#   prefix: /conda-envs/4ff5cdad114a5133b20a8ce52f3ad433
#   name: R
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - bioconductor-viper
#     # - bioconductor-dorothea (NOTE: DOrothea is broken on bioconda)
#     - bioconductor-decoupler
#     - bioconductor-chipseeker
#     - bioconductor-deseq2
#     - bioconductor-omnipathr
#     - bioconductor-reactomepa
#     - bioconductor-genomicranges
#     - bioconductor-pcatools
#     - bioconductor-carnival
#     - bioconductor-txdb.mmusculus.UCSC.mm10.knownGene
#     - bioconductor-ensembldb
#     - bioconductor-progeny
#     - r-magrittr
#     - bioconductor-rtracklayer
#     - bioconductor-annotationhub
#     - bioconductor-org.mm.eg.db
#     - bioconductor-org.hs.eg.db
#     - bioconductor-mitch
#     - r-ggplot2
#     - r-tidyverse
#     - r-reshape2
#     - r-knitr
#     - r-data.table
#     - r-dt
#     - r-devtools
#     - bioconductor-complexheatmap
#     #- bioconductor-bsgenome.mmusculus.ucsc.mm10
#     #- bioconductor-bsgenome.hsapiens.ucsc.hg38
#     - bioconductor-clusterprofiler
#     - r-rmarkdown
#     - r-knitr
#     - r-readxl
#     - r-openxlsx
#     - r-rcpp
#     - r-pheatmap
#     - r-msigdbr
#     - r-babelgene
#     - r-ggsignif
#     - r-rvcheck
#     - bioconductor-enhancedvolcano
#     - r-ggnewscale
#     - r-svglite
#     - bioconductor-do.db
#     - bioconductor-apeglm
#     - r-ggupset
#     - r-furrr
RUN mkdir -p /conda-envs/4ff5cdad114a5133b20a8ce52f3ad433
COPY workflow/envs/R_4.yaml /conda-envs/4ff5cdad114a5133b20a8ce52f3ad433/environment.yaml
COPY workflow/envs/R_4.post-deploy.sh /conda-envs/4ff5cdad114a5133b20a8ce52f3ad433/post-deploy.sh

# Conda environment:
#   source: workflow/envs/deseq2.yaml
#   prefix: /conda-envs/42d453487007f54c94a3bbe3dda8980a
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - r-base
#     - bioconductor-deseq2
#     - bioconductor-annotate 
#     - bioconductor-apeglm 
#     - bioconductor-annotationdbi 
#     - bioconductor-biobase 
#     - bioconductor-biocgenerics
#     - bioconductor-biocparallel
#     - bioconductor-clusterprofiler
#     - bioconductor-complexheatmap
#     - bioconductor-delayedarray
#     - bioconductor-do.db
#     - bioconductor-edger
#     - bioconductor-enhancedvolcano
#     - bioconductor-genefilter
#     - bioconductor-geneplotter
#     - bioconductor-genomeinfodb
#     - bioconductor-genomeinfodbdata
#     - bioconductor-genomicranges
#     - bioconductor-goseq
#     - bioconductor-iranges
#     - bioconductor-org.mm.eg.db
#     - bioconductor-org.hs.eg.db
#     - bioconductor-pcatools
#     - bioconductor-reactome.db
#     - bioconductor-reactomepa
#     - bioconductor-s4vectors
#     - bioconductor-summarizedexperiment 
#     - bioconductor-tximport 
#     - bioconductor-xvector 
#     - bioconductor-zlibbioc 
#     - r-acepack 
#     - r-ashr
#     - r-backports 
#     - r-base64enc 
#     - r-bh 
#     - r-bit 
#     - r-bit64 
#     - r-bitops 
#     - r-blob 
#     - r-catools 
#     - r-checkmate 
#     - r-cluster 
#     - r-colorspace 
#     - r-devtools
#     - r-DT
#     - r-dbi 
#     - r-dichromat 
#     - r-digest 
#     - r-dplyr
#     - r-evaluate 
#     - r-foreign 
#     - r-formula 
#     - r-futile.logger 
#     - r-futile.options 
#     - r-gdata 
#     - r-getopt 
#     - r-ggplot2 
#     - r-ggupset
#     - r-gplots 
#     - r-gridextra 
#     - r-gtable 
#     - r-gtools 
#     - r-highr 
#     - r-hmisc 
#     - r-htmltable 
#     - r-htmltools 
#     - r-htmlwidgets 
#     - r-jsonlite 
#     - r-kernsmooth 
#     - r-knitr 
#     - r-labeling 
#     - r-lambda.r 
#     - r-lattice 
#     - r-latticeextra 
#     - r-lazyeval 
#     - r-locfit 
#     - r-magrittr 
#     - r-markdown 
#     - r-mass 
#     - r-matrix 
#     - r-matrixstats 
#     - r-memoise 
#     - r-mime 
#     - r-munsell 
#     - r-msigdbr
#     - r-nnet 
#     - r-patchwork
#     - r-pheatmap
#     - r-pkgconfig 
#     - r-plogr 
#     - r-plyr 
#     - r-rvcheck
#     - r-rcolorbrewer 
#     - r-rcpp 
#     - r-rcpparmadillo 
#     - r-rcurl 
#     - r-readr
#     - r-reshape2 
#     - r-rjson 
#     - r-rlang 
#     - r-rmarkdown
#     - r-rpart 
#     - r-rsqlite 
#     - r-scales 
#     - r-snow 
#     - r-stringi 
#     - r-stringr 
#     - r-survival 
#     - r-tibble 
#     - r-tidyverse
#     - r-upsetr
#     - r-viridis 
#     - r-viridislite 
#     - r-xml 
#     - r-xtable 
#     - r-yaml
RUN mkdir -p /conda-envs/42d453487007f54c94a3bbe3dda8980a
COPY workflow/envs/deseq2.yaml /conda-envs/42d453487007f54c94a3bbe3dda8980a/environment.yaml
COPY workflow/envs/deseq2.post-deploy.sh /conda-envs/42d453487007f54c94a3bbe3dda8980a/post-deploy.sh

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/4ff5cdad114a5133b20a8ce52f3ad433 --file /conda-envs/4ff5cdad114a5133b20a8ce52f3ad433/environment.yaml && \
    mamba env create --prefix /conda-envs/42d453487007f54c94a3bbe3dda8980a --file /conda-envs/42d453487007f54c94a3bbe3dda8980a/environment.yaml && \
    mamba clean --all -y

# Step 3 RUn post deploy scripts for each env
RUN mamba activate /conda-envs/4ff5cdad114a5133b20a8ce52f3ad433 && \
    /bin/bash /conda-envs/4ff5cdad114a5133b20a8ce52f3ad433/post-deploy.sh && \
    mamba deactivate && \
    mamba activate /conda-envs/42d453487007f54c94a3bbe3dda8980a && \
    /bin/bash /conda-envs/42d453487007f54c94a3bbe3dda8980a/post-deploy.sh && \
    mamba deactivate
RUN mamba activate 

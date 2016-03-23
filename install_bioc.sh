#!/bin/bash

Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("BiocInstaller"); biocLite("BSgenome.Hsapiens.UCSC.hg19"); install.packages("devtools"); biocLite("GenomicRanges"); devtools::install_github("mskilab/gUtils"); devtools::install_github("mskilab/skitools")'

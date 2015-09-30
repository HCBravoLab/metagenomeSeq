metagenomeSeq
=============

Statistical analysis for sparse high-throughput sequencing

[![Travis-CI Build Status](https://travis-ci.org/HCBravoLab/metagenomeSeq.svg?branch=master)](https://travis-ci.org/HCBravoLab/metagenomeSeq)

metagenomeSeq is designed to determine features (be it Operational Taxanomic Unit (OTU), species, etc.) 
that are differentially abundant between two or more groups of multiple samples. metagenomeSeq is designed 
to address the effects of both normalization and undersampling of microbial communities on disease 
association detection and the testing of feature correlations.

To install the latest release version of metagenomeSeq:
```S
source("http://bioconductor.org/biocLite.R")
biocLite("metagenomeSeq")
```

To install the latest development version of metagenomeSeq:
```S
install.packages("devtools")
library("devtools")
install_github("metagenomeSeq","nosson")
```

Author: [Joseph Nathaniel Paulson](http://www.cbcb.umd.edu/~jpaulson), Hisham Talukder, [Mihai Pop](http://www.cbcb.umd.edu/~mpop), [Hector Corrada Bravo](http://www.cbcb.umd.edu/~hcorrada)

Maintainer: Joseph N. Paulson : jpaulson at umiacs.umd.edu

Website: www.cbcb.umd.edu/software/metagenomeSeq

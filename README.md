metagenomeSeq
=============

Statistical analysis for sparse high-throughput sequencing

<a href="http://www.bioconductor.org/packages/devel/bioc/html/metagenomeSeq.html#since"><img border="0" src="http://www.bioconductor.org/shields/years-in-bioc/metagenomeSeq.svg" title="How long since the package was first in a released Bioconductor version (or is it in devel only)."></a> <a href="http://bioconductor.org/packages/stats/bioc/metagenomeSeq.html"><img border="0" src="http://www.bioconductor.org/shields/downloads/metagenomeSeq.svg" title="Percentile (top 5/20/50% or 'available') of downloads over last 6 full months. Comparison is done across all package categories (software, annotation, experiment)."></a> <a href="https://support.bioconductor.org/t/metagenomeSeq/"><img border="0" src="http://www.bioconductor.org/shields/posts/metagenomeSeq.svg" title="Support site activity, last 6 months: tagged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts."></a> <a href="http://www.bioconductor.org/packages/devel/bioc/html/metagenomeSeq.html#svn_source"><img border="0" src="http://www.bioconductor.org/shields/commits/bioc/metagenomeSeq.svg" title="average Subversion commits (to the devel branch) per month for the last 6 months"></a>

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
install_github("Bioconductor-mirror/metagenomeSeq")
```

Author: [Joseph Nathaniel Paulson](http://bcb.dfci.harvard.edu/~jpaulson), Hisham Talukder, [Mihai Pop](http://www.cbcb.umd.edu/~mpop), [Hector Corrada Bravo](http://www.cbcb.umd.edu/~hcorrada)

Maintainer: Joseph N. Paulson : jpaulson at jimmy.harvard.edu

Website: www.cbcb.umd.edu/software/metagenomeSeq

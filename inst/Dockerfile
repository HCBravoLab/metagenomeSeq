# Docker image to check metagenomeSeq with Bioc devel
FROM bioconductor/bioconductor_docker:devel

# Install all the latex stuff to build vignettes
RUN apt-get update \
	&& apt-get install -y --no-install-recommends apt-utils \
	&& apt-get install -y --no-install-recommends \
	texlive \
	texlive-latex-extra \
	texlive-fonts-extra \
	texlive-bibtex-extra \
	texlive-science \
	texi2html \
	texinfo \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

## Install BiocStyle
RUN R -e 'BiocManager::install("BiocStyle")'

# Install metagenomeSeq
RUN R -e 'BiocManager::install("metagenomeSeq")'

# Now dependencies to build vignettes
RUN R -e 'BiocManager::install(c("biomformat", "gss"))'

# build with command: docker build -t biocondcutor_docker_metagenomeseq:devel inst
# run rstudio with command: docker run -e PASSWORD=<password> -p 8787:8787 -v $PWD:/metagenomeSeq bioconductor_docker_metagenomeseq:devel

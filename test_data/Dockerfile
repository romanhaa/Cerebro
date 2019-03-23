FROM r-base:3.5.2

LABEL org.label-schema.name = "cerebro-example"
LABEL org.label-schema.description = "\
Environment for analysis of example data set for Cerebro."
LABEL org.label-schema.url = "https://github.com/romanhaa/Cerebro"
LABEL org.label-schema.schema-version = "1.0"
LABEL maintainer = "roman.hillje@googlemail.com"
LABEL org.label-schema.version = "1.0"
LABEL org.label-schema.build-date = "2018-01-11T00:00:00.0Z"

RUN apt update -y && apt upgrade -y
RUN apt install -y \
  python-dev \
  python-pip \
  python3-dev \
  python3-pip \
  python3-tk

# curl and RCurl (R package)
RUN apt install -y libcurl4-openssl-dev

# openssl (R package)
RUN apt install -y libssl-dev

# Cairo (R package)
RUN apt install -y libcairo2-dev libxt-dev

# XML and xml2 (R package)
RUN apt install -y libxml2-dev

# hdf5r (R package)
RUN apt install -y libhdf5-dev

# for UMAP in Seurat (R package)
# RUN pip install umap-learn
RUN pip3 install umap-learn
RUN pip3 install "numba<0.43" --force-reinstall

# for phate in Seurat (R package)
# not available for python2
RUN pip3 install phate

# use script here as soon as repo is public
# ADD https://github.com/romanhaa/cerebro-prepare/archive/master.tar.gz /cerebro/
# CMD Rscript install_packages.R
RUN Rscript -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install(c("destiny","phateR","reshape2","scales","sessioninfo","Seurat","tidyverse"))'
RUN Rscript -e 'BiocManager::install("devtools")'
RUN Rscript -e 'devtools::install_github("romanhaa/cerebroPrepare")'

CMD R

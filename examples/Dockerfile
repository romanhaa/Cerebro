FROM r-base:3.5.2

##----------------------------------------------------------------------------##
## LABELS
##----------------------------------------------------------------------------##

LABEL org.label-schema.name = "cerebro-example"
LABEL org.label-schema.description = "\
Environment for analysis of example data set for Cerebro."
LABEL org.label-schema.url = "https://github.com/romanhaa/Cerebro"
LABEL org.label-schema.schema-version = "1.0"
LABEL maintainer = "roman.hillje@googlemail.com"
LABEL org.label-schema.version = "1.0"
LABEL org.label-schema.build-date = "2019-04-29T00:00:00.0Z"

##----------------------------------------------------------------------------##
## SYSTEM LIBRARIES
##----------------------------------------------------------------------------##

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

##----------------------------------------------------------------------------##
## PYTHON PACKAGES
##----------------------------------------------------------------------------##

# for UMAP in Seurat (R package)
RUN pip3 install umap-learn
RUN pip3 install "numba<0.43" --force-reinstall

# for phate in Seurat (R package)
RUN pip3 install phate

##----------------------------------------------------------------------------##
## R PACKAGES
##----------------------------------------------------------------------------##

# base stuff + Seurat v3
RUN Rscript -e 'install.packages("BiocManager", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install(c("destiny","devtools","phateR","reshape2","scales","Seurat","sessioninfo","tidyverse"))'

# cerebroPrepare
RUN Rscript -e 'devtools::install_github("romanhaa/cerebroPrepare")'

# Seurat v2.3.4
RUN mkdir /other_R_packages
RUN Rscript -e 'devtools::install_version("Seurat", version = "2.3.4", repos = "http://cran.us.r-project.org", lib = "/other_R_packages")

CMD R

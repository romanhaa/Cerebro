FROM r-base:3.6.0

##----------------------------------------------------------------------------##
## LABELS
##----------------------------------------------------------------------------##

LABEL org.label-schema.name = "cerebro-example"
LABEL org.label-schema.description = "Environment for analysis of example data set for Cerebro."
LABEL org.label-schema.url = "https://github.com/romanhaa/Cerebro"
LABEL org.label-schema.schema-version = "1.0"
LABEL maintainer = "roman.hillje@googlemail.com"
LABEL org.label-schema.version = "1.0.0"
LABEL org.label-schema.build-date = "2019-04-29T00:00:00.0Z"

##----------------------------------------------------------------------------##
## SYSTEM LIBRARIES
##----------------------------------------------------------------------------##

RUN apt-get update --yes \
 && apt-get upgrade --yes \
 && apt-get install --yes \
  python-dev \
  python-pip \
  python3-dev \
  python3-pip \
  python3-tk \
  libcurl4-openssl-dev \
  libssl-dev \
  libcairo2-dev \
  libxt-dev \
  libxml2-dev \
  libhdf5-dev

##----------------------------------------------------------------------------##
## PYTHON PACKAGES
##----------------------------------------------------------------------------##

RUN pip3 install umap-learn phate

##----------------------------------------------------------------------------##
## R PACKAGES
##----------------------------------------------------------------------------##

# dependencies + Seurat v3
RUN Rscript -e 'install.packages("BiocManager", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install(c("destiny","devtools","phateR","reshape2","scales","Seurat","sessioninfo","tidyverse","monocle"))'

# cerebroPrepare
RUN Rscript -e 'devtools::install_github("romanhaa/cerebroPrepare")'

# Seurat v2.3.4
RUN mkdir /other_R_packages
RUN Rscript -e 'devtools::install_version("Seurat", version = "2.3.4", repos = "http://cran.us.r-project.org", lib = "/other_R_packages")'

CMD R

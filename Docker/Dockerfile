FROM r-base:3.6.1

##----------------------------------------------------------------------------##
## LABELS
##----------------------------------------------------------------------------##

LABEL org.label-schema.name = "cerebro"
LABEL org.label-schema.description = "Cerebro as Docker container."
LABEL org.label-schema.url = "https://github.com/romanhaa/Cerebro"
LABEL org.label-schema.schema-version = "1.1.0"
LABEL maintainer = "roman.hillje@googlemail.com"
LABEL org.label-schema.version = "1.1.0"
LABEL org.label-schema.build-date = "2019-09-23T00:00:00.0Z"

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

RUN pip3 install louvain python-igraph scanpy umap-learn

##----------------------------------------------------------------------------##
## R PACKAGES
##----------------------------------------------------------------------------##

# dependencies + Seurat v3 + monocle
RUN Rscript -e 'install.packages("BiocManager", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install(c("devtools","remotes"), update = TRUE, ask = FALSE)'
RUN Rscript -e 'BiocManager::install("satijalab/seurat@develop")'
RUN Rscript -e 'BiocManager::install("monocle")'

# cerebroApp
RUN Rscript -e 'devtools::install_github("romanhaa/cerebroApp@merge_packages")'

# Seurat v2.3.4
RUN mkdir /other_R_packages
RUN Rscript -e 'devtools::install_version("Seurat", version = "2.3.4", repos = "http://cran.us.r-project.org", lib = "/other_R_packages")'

##----------------------------------------------------------------------------##
## FINISHING
##----------------------------------------------------------------------------##

RUN mkdir /plots

EXPOSE 8080

ENV DOCKER=TRUE

CMD Rscript -e 'shiny::runApp(cerebroApp::launchCerebro(), port=8080, host="0.0.0.0", launch.browser=FALSE)'

# https://www.bjoern-hartmann.de/post/learn-how-to-dockerize-a-shinyapp-in-7-steps/
# docker run -p 8080:8080 romanhaa/cerebro:v1.1

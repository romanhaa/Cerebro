FROM r-base:3.6.2

##----------------------------------------------------------------------------##
## LABELS
##----------------------------------------------------------------------------##

LABEL org.label-schema.name = "cerebro"
LABEL org.label-schema.description = "Cerebro as Docker container."
LABEL org.label-schema.url = "https://github.com/romanhaa/Cerebro"
LABEL org.label-schema.schema-version = "1.2.1"
LABEL maintainer = "roman.hillje@googlemail.com"
LABEL org.label-schema.version = "1.2.1"
LABEL org.label-schema.build-date = "2020-04-04T00:00:00.0Z"

##----------------------------------------------------------------------------##
## SYSTEM LIBRARIES
##----------------------------------------------------------------------------##

RUN apt-get update --yes \
 && apt-get upgrade --yes \
 && apt-get install --yes \
  python3-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  libcairo2-dev \
  libxt-dev \
  libxml2-dev \
  libhdf5-dev

##----------------------------------------------------------------------------##
## R PACKAGES
##----------------------------------------------------------------------------##

# install some required packages
RUN Rscript -e 'install.packages("BiocManager", Ncpus = 2, repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install(c("devtools","remotes"), Ncpus = 2, update = TRUE, ask = FALSE)'

# cerebroApp
RUN Rscript -e 'devtools::install_github("romanhaa/cerebroApp", Ncpus = 2)'

##----------------------------------------------------------------------------##
## FINISHING
##----------------------------------------------------------------------------##

RUN mkdir /plots

EXPOSE 8080

ENV DOCKER=TRUE

CMD Rscript -e 'shiny::runApp(cerebroApp::launchCerebro(), port=8080, host="0.0.0.0", launch.browser=FALSE)'

# https://www.bjoern-hartmann.de/post/learn-how-to-dockerize-a-shinyapp-in-7-steps/
# docker run -p 8080:8080 romanhaa/cerebro:v1.1

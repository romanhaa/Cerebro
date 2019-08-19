FROM r-base:3.6.0

##----------------------------------------------------------------------------##
## LABELS
##----------------------------------------------------------------------------##

LABEL org.label-schema.name = "cerebro"
LABEL org.label-schema.description = "Cerebro as Docker container."
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
  libcurl4-openssl-dev libssl-dev

##----------------------------------------------------------------------------##
## R PACKAGES
##----------------------------------------------------------------------------##

RUN Rscript -e 'install.packages(c("BiocManager","remotes"))'
RUN Rscript -e 'BiocManager::install("romanhaa/cerebroApp")'

##----------------------------------------------------------------------------##
## FINISHING
##----------------------------------------------------------------------------##

RUN mkdir /plots

EXPOSE 8080

ENV DOCKER=TRUE

CMD Rscript -e 'shiny::runApp(cerebroApp::launchApp(), port=8080, host="0.0.0.0", launch.browser=FALSE)'

# https://www.bjoern-hartmann.de/post/learn-how-to-dockerize-a-shinyapp-in-7-steps/
# docker run -p 8080:8080 romanhaa/cerebro:1.0

FROM r-base:3.5.2

RUN apt-get update --yes && \
  apt-get upgrade --yes

RUN apt-get install --yes \
  libcurl4-openssl-dev \
  libssl-dev

RUN Rscript -e 'install.packages(c("BiocManager","remotes"))'
RUN Rscript -e 'BiocManager::install("romanhaa/cerebroApp", auth_token = "c21eae3bd91dd22b551915c330ed97f45352432e")'

EXPOSE 8080

CMD Rscript -e 'shiny::runApp(cerebroApp::launchApp(), port=8080, host="0.0.0.0", launch.browser=FALSE)'

# https://www.bjoern-hartmann.de/post/learn-how-to-dockerize-a-shinyapp-in-7-steps/
# docker run -p 8080:8080 romanhaa/cerebro:1.0
# add ubuntu 16.04 image to match the OS of the machine on which ensemblQueryR was developed
FROM ubuntu:16.04

RUN apt-get update
RUN apt-get install -y wget
RUN apt-get install -y build-essential

RUN apt-get -y install libxml2-dev
RUN export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/usr/lib/x86_64-linux-gnu/pkgconfig/

# RUN wget ftp://xmlsoft.org/libxml2/libxml2-2.9.2.tar.gz; tar -xvf libxml2-2.9.2.tar.gz; cd libxml2-2.9.2; ./configure --prefix=/usr --disable-static --with-history && make; make install

RUN apt-get -y install r-base-core
RUN apt-get -y install libssl-dev
RUN apt-get -y install libcurl4-openssl-dev
RUN apt-get -y install openssl
RUN apt-get -y install libfontconfig1-dev
RUN apt-get -y install libharfbuzz-dev
RUN apt-get -y install libfribidi-dev
RUN apt-get -y install libfreetype6-dev
RUN apt-get -y install libpng-dev
RUN apt-get -y install libtiff5-dev
RUN apt-get -y install libjpeg-dev

# add r-base 4.0.5 image to match the R version on which ensemblQueryR was developed
FROM r-base:4.0.5

# install ensemblQueryR dependencies at specific versions used in development
RUN R -e "install.packages('remotes')"
RUN R -e "library(remotes)"

RUN R -e "remotes::install_github('r-lib/httr')"
RUN R -e "remotes::install_github('r-lib/xml2')"

# RUN R -e "install.packages('dplyr')"
# RUN R -e "install.packages('jsonlite')"
# RUN R -e "install.packages('purrr')"
# RUN R -e "install.packages('tidyr')"
# RUN R -e "install.packages('vroom')"
# RUN R -e "install.packages('magrittr')"

# load ensemblQueryR dependencies
RUN R -e "library(httr)"
RUN R -e "library(xml2)"
RUN R -e "library(dplyr)"
RUN R -e "library(jsonlite)"
RUN R -e "library(parallel)"
RUN R -e "library(purrr)"
RUN R -e "library(tidyr)"
RUN R -e "library(vroom)"
RUN R -e "library(magrittr)"
RUN R -e "library(remotes)"

# install and load ensemblQueryR
# RUN R -e "remotes::install_github('ainefairbrother/ensemblQueryR')"
# RUN R -e "library(ensemblQueryR)"

# clear the R console
# RUN R -e "clear"

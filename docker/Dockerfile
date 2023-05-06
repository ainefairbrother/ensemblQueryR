# use rocker 4.0.0 and Ubuntu 20, as these resemble closely the development OS for ensemblQueryR
FROM rocker/r-ver:4.0.0

# system libraries of general use - informed by https://tinyurl.com/ypnkbrpb
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    pkg-config \
    libnlopt-dev \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libgsl-dev \
    libssl-dev \
    libssh2-1-dev \
    libssl1.0 \
    libxml2-dev \
    openssl
RUN apt-get update && apt-get install -y \
    libmpfr-dev
RUN apt-get install libcurl4-openssl-dev -y

# install and load ensemblQueryR dependencies - package versions are those used at the time of package development
RUN R -e "install.packages('remotes'); \
library(remotes); \
install.packages('httr'); \
install.packages('xml2'); \
install.packages('dplyr'); \
install.packages('jsonlite'); \
install.packages('purrr'); \
install.packages('tidyr'); \
install.packages('vroom'); \
install.packages('magrittr'); \
install.packages('parallel'); \
library(httr); \
library(xml2); \
library(dplyr); \
library(jsonlite); \
library(purrr); \
library(tidyr); \
library(vroom); \
library(magrittr); \
library(parallel)"

# install and load ensemblQueryR - install from clone
RUN apt-get install -y git
RUN git clone https://github.com/ainefairbrother/ensemblQueryR.git ./ensemblQueryR
RUN R -e "install.packages('devtools'); devtools::install('./ensemblQueryR'); library(ensemblQueryR)"

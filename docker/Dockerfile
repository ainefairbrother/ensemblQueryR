# add ubuntu 16.04 image to match the OS of the machine on which ensemblQueryR was developed
FROM ubuntu:16.04

# https://askubuntu.com/questions/1086870/r-script-returns-error-in-libraryhttr-there-is-no-package-called-httr
RUN apt-get update
RUN apt-get -y install r-base-core
RUN apt-get -y install libssl-dev
RUN apt-get -y install libcurl4-openssl-dev
RUN apt-get -y install xml2
RUN apt-get -y install openssl

# add r-base 4.0.5 image to match the R version on which ensemblQueryR was developed
FROM r-base:4.0.5

# install ensemblQueryR dependencies at specific versions used in development
RUN R -e "install.packages('devtools')"
RUN R -e "library(devtools)"
RUN R -e "devtools::install_github('r-lib/xml2')"
RUN R -e "devtools::install_github('r-lib/httr')"

# RUN R -e "install.packages('remotes')"
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

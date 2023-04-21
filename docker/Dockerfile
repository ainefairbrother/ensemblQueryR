# add ubuntu 16.04 image to match the OS of the machine on which ensemblQueryR was developed
FROM ubuntu:16.04

# https://askubuntu.com/questions/1086870/r-script-returns-error-in-libraryhttr-there-is-no-package-called-httr
RUN apt-get update
RUN apt-get -y install r-base-core libssl-dev libcurl4-openssl-dev
RUN apt-get -y install xml2 openssl

# add r-base 4.0.5 image to match the R version on which ensemblQueryR was developed
FROM r-base:4.0.5

# install ensemblQueryR dependencies at specific versions used in development
# RUN R -e "install.packages('devtools', dependencies= T)"
# RUN R -e "install.packages('remotes', dependencies= T)"
# RUN R -e "install.packages('dplyr', dependencies= T)"
# RUN R -e "install.packages('jsonlite', dependencies= T)"
# RUN R -e "install.packages('purrr', dependencies= T)"
# RUN R -e "install.packages('tidyr', dependencies= T)"
# RUN R -e "install.packages('vroom', dependencies= T)"
# RUN R -e "install.packages('magrittr', dependencies= T)"

RUN R -e "install.packages('xml2')"
RUN R -e "install.packages('httr')"

# load ensemblQueryR dependencies
RUN R -e "library(httr)"
RUN R -e "library(xml2)"
# RUN R -e "library(dplyr)"
# RUN R -e "library(jsonlite)"
# RUN R -e "library(parallel)"
# RUN R -e "library(purrr)"
# RUN R -e "library(tidyr)"
# RUN R -e "library(vroom)"
# RUN R -e "library(magrittr)"
# RUN R -e "library(remotes)"

# install and load ensemblQueryR
# RUN R -e "remotes::install_github('ainefairbrother/ensemblQueryR')"
# RUN R -e "library(ensemblQueryR)"

# clear the R console
# RUN R -e "clear"

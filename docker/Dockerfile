# add ubuntu 16.04 image to match the OS of the machine on which ensemblQueryR was developed
FROM ubuntu:16.04

# add r-base 4.0.5 image to match the R version on which ensemblQueryR was developed
FROM r-base:4.0.5

# sources required for libcurl4-gnutls-dev install
# RUN cp /etc/apt/sources.list /etc/apt/sources.list~
# RUN sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list
# RUN apt update

# libcurl4-gnutls-dev is a requirement of devtools, so install first
# RUN apt-get -y build-dep libcurl4-gnutls-dev
# RUN apt-get -y install libcurl4-gnutls-dev

# install ensemblQueryR dependencies at specific versions used in development
#RUN R -e "install.packages('devtools', dependencies= T)"

#RUN R -e "library(devtools)"

RUN R -e "install.packages('remotes', dependencies= T)"
RUN R -e "remotes::install_version(package = 'remotes', version = '2.4.2', dependencies= T)"
RUN R -e "remotes::install_version(package = 'dplyr', version = '1.1.1', dependencies= T)"
RUN R -e "remotes::install_version(package = 'httr', version = '1.4.5', dependencies= T)"
RUN R -e "remotes::install_version(package = 'jsonlite', version = '1.8.4', dependencies= T)"
RUN R -e "remotes::install_version(package = 'purrr', version = '1.0.1', dependencies= T)"
RUN R -e "remotes::install_version(package = 'tidyr', version = '1.3.0', dependencies= T)"
RUN R -e "remotes::install_version(package = 'vroom', version = '1.6.1', dependencies= T)"
RUN R -e "remotes::install_version(package = 'xml2', version = '1.3.3', dependencies= T)"
RUN R -e "remotes::install_version(package = 'magrittr', version = '2.0.3', dependencies= T)"

RUN R -e "library(devtools)"
RUN R -e "library(dplyr)"
RUN R -e "library(httr)"
RUN R -e "library(jsonlite)"
RUN R -e "library(parallel)"
RUN R -e "library(purrr)"
RUN R -e "library(tidyr)"
RUN R -e "library(vroom)"
RUN R -e "library(xml2)"
RUN R -e "library(magrittr)"
RUN R -e "library(remotes)"
RUN R -e "remotes::install_github("ainefairbrother/ensemblQueryR")"
RUN R -e "library(ensemblQueryR)"

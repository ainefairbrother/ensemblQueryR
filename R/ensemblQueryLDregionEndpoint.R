#' Function to query Ensembl LD data within a genomic window.
#'
#' @param chr String. Chromosome that the query region is located on.
#' @param start String. Base pair that the query region starts at.
#' @param end String. Base pair that the query region ends at.
#' @param pop String. Population for which to compute LD. Use `ensemblQueryGetPops()` to retrieve a list of all populations with LD data. Default is 1000GENOMES:phase_3:EUR.
#'
#' @return
#'
#' @import httr
#' @import xml2
#' @import jsonlite
#' @import dplyr
#' @import tidyr
#' @import vroom
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' ensemblQueryLDwithSNPregion(
#' chr="6",
#' start="25837556",
#' end="25843455",
#' pop="1000GENOMES:phase_3:EUR"
#' )
#'
ensemblQueryLDwithSNPregion = function(chr, start, end, pop="1000GENOMES:phase_3:EUR"){

  # # TEST
  # # load libs
  # require(httr)
  # require(xml2)
  # if( !("tidyverse" %in% (.packages())) ){
  #   require(jsonlite)
  # }
  # require(dplyr)
  # require(tidyr)
  # require(purrr)
  # require(vroom)
  # require(magrittr)

  # chr=6
  # start=25837556
  # end=25843455 #25843455
  # pop="1000GENOMES:phase_3:EUR"

  #------------------------------ check inputs -------------------------------

  stopifnot(is.character(chr) | is.numeric(chr))
  stopifnot(is.character(start) | is.numeric(start))
  stopifnot(is.character(end) | is.numeric(end))
  stopifnot(is.character(pop))

  #--------------------------------- run query -------------------------------

  server <- "https://rest.ensembl.org"
  ext <- paste0("/ld/human/region/",chr,":",start,"..",end,"/",pop,"?")

  r <- httr::GET(url=paste(server, ext, sep = ""), content_type("application/json"))

  #-------------------- check output and write out ---------------------------

  # stop_for_status(r)

  # error handling, if 400 error, set res.temp as NA
  if(r$status_code == 400){
    print("Error 400 thrown by httr::GET. This may not be a valid SNP rsID, check using dbSNP: https://www.ncbi.nlm.nih.gov/snp/.")
    res.temp = NA
  } else{
    # if no error, use this if you get a simple nested list back, otherwise inspect its structure
    res.temp = jsonlite::fromJSON(jsonlite::toJSON(content(r))) %>%
      data.frame()
  }

  # deal with null search result by testing for df length
  # then returning an empty df if no result
  if(is.data.frame(res.temp)){
    if(nrow(res.temp)==0){

      print(paste0("No LD data found for the genomic coordinates input: ", chr,":",start,"-",end))

      data.frame(rep(NA, 5), row.names = c("variation1", "variation2", "d_prime", "population_name", "r2")) %>%
        t() %>%
        `rownames<-`(NULL) %>%
        as.data.frame() %>%
        dplyr::rename(rsid1=variation1, rsid2=variation2) %>%
        dplyr::relocate(rsid1, rsid2, r2, d_prime, population_name) %>%
        return()
    } else{
      # if not 0-row (empty) df, then deal with it normally, format and prepare for return
      res.temp %>%
        data.frame() %>%
        dplyr::arrange(r2) %>%
        dplyr::rename(rsid1=variation1, rsid2=variation2) %>%
        dplyr::relocate(rsid1, rsid2, r2, d_prime, population_name) %>%
        return()
    }
    # deal with NA search result (result of 400 error) by testing if res.temp is NA
    # then returning an empty df of the same structure
  } else{
    if(is.na(res.temp)){

      print(paste0("No LD data found for the genomic coordinates input: ", chr,":",start,"-",end))

      data.frame(rep(NA, 5), row.names = c("variation1", "variation2", "d_prime", "population_name", "r2")) %>%
        t() %>%
        `rownames<-`(NULL) %>%
        as.data.frame() %>%
        dplyr::rename(rsid1=variation1, rsid2=variation2) %>%
        dplyr::relocate(rsid1, rsid2, r2, d_prime, population_name) %>%
        return()
    }
  }
}

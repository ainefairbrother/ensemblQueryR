#' Function to get list of populations that Ensembl has available to query LD.
#'
#' @return data.frame of populations.
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' ensemblQueryGetPops()
#'
ensemblQueryGetPops = function(){

  #--------------------------------- get pops --------------------------------

  server <- "https://rest.ensembl.org"
  ext <- "/info/variation/populations/homo_sapiens?filter=LD"

  httr::GET("http://cran.r-project.org/Rlogo.jpg")

  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))

  jsonlite::fromJSON(jsonlite::toJSON(content(r))) %>%
    data.frame() %>%
    return(.)

}

#' Function to query Ensembl LD data with a single rsID
#'
#' @param rsid String. Variant ID.
#' @param r2 Float. Measure of LD. If r-squared is provided only return pairs of variants whose r-squared value is equal to or greater than the value provided.
#' @param d.prime Float. Measure of LD. If D' is provided only return pairs of variants whose D' value is equal to or greater than the value provided.
#' @param window.size Integer. Window size in kb. The maximum allowed value for the window size is 500 kb. LD is computed for the given variant and all variants that are located within the specified window.
#' @param pop String. Population for which to compute LD. Use `ensemblQueryGetPops()` to retrieve a list of all populations with LD data. Default is 1000GENOMES:phase_3:EUR.
#'
#' @return data.frame with 5 columns.
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
#' ensemblQueryLDwithSNPwindow(rsid="rs3851179", r2=0.8, d.prime=0.8, window.size=500, pop="1000GENOMES:phase_3:EUR")
#'
ensemblQueryLDwithSNPwindow = function(rsid, r2=0.8, d.prime=0.8, window.size=500, pop="1000GENOMES:phase_3:EUR"){

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
  #
  # rsid="rs1210721423"
  # r2=0.8
  # d.prime=0.8
  # window.size=500
  # pop="1000GENOMES:phase_3:EUR"

  #--------------------------------- run query -------------------------------

  server <- "https://rest.ensembl.org"
  ext <- paste0("/ld/human/",rsid,"/",pop,"?d_prime=",d.prime,";window_size=",window.size,";r2=",r2)

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
      data.frame(rep(NA, 5), row.names = c("variation2", "population_name", "variation1", "d_prime", "r2")) %>%
        t() %>%
        `rownames<-`(NULL) %>%
        as.data.frame() %>%
        dplyr::rename(query=variation1, snp_in_ld=variation2) %>%
        dplyr::relocate(query, snp_in_ld, r2, d_prime, population_name) %>%
        dplyr::mutate(query=rsid) %>%
        return()
    } else{
      # if not 0-row (empty) df, then deal with it normally, format and prepare for return
      res.temp %>%
        data.frame() %>%
        dplyr::arrange(r2) %>%
        dplyr::rename(query=variation1, snp_in_ld=variation2) %>%
        dplyr::relocate(query, snp_in_ld, r2, d_prime, population_name) %>%
        return()
    }
    # deal with NA search result (result of 400 error) by testing if res.temp is NA
    # then returning an empty df of the same structure
  } else{
    if(is.na(res.temp)){
      data.frame(rep(NA, 5), row.names = c("variation2", "population_name", "variation1", "d_prime", "r2")) %>%
        t() %>%
        `rownames<-`(NULL) %>%
        as.data.frame() %>%
        dplyr::rename(query=variation1, snp_in_ld=variation2) %>%
        dplyr::relocate(query, snp_in_ld, r2, d_prime, population_name) %>%
        dplyr::mutate(query=rsid) %>%
        return()
    }
  }
}

#' Function to apply `ensemblQueryLDwithSNPwindow` to a list of SNPs of length<1000
#' Outputs a table containing all SNPs in LD with all query SNPs in rsid.list.
#' Outputs NA row if no SNPs in LD found.
#'
#' @param rsid.list Character vector of variant IDs.
#' @param r2 Float. Measure of LD. If r-squared is provided only return pairs of variants whose r-squared value is equal to or greater than the value provided.
#' @param d.prime Float. Measure of LD. If D' is provided only return pairs of variants whose D' value is equal to or greater than the value provided.
#' @param window.size Integer. Window size in kb. The maximum allowed value for the window size is 500 kb. LD is computed for the given variant and all variants that are located within the specified window.
#' @param pop String. Population for which to compute LD. Use `ensemblQueryGetPops()` to retrieve a list of all populations with LD data. Default is 1000GENOMES:phase_3:EUR.
#'
#' @return
#'
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' ensemblQueryLDwithSNPwindowList(rsid.list=c("rs7153434","rs1963154","rs12672022","rs3852802","rs12324408","rs56346870"),
#'                           r2=0.8,
#'                           d.prime=0.8,
#'                           window.size=500,
#'                           pop="1000GENOMES:phase_3:EUR")
#'
ensemblQueryLDwithSNPwindowList = function(rsid.list, r2=0.8, d.prime=0.8, window.size=500, pop="1000GENOMES:phase_3:EUR"){

  # max query length ensembl REST API will accept
  max.query.len=1000

  if(length(rsid.list)<=max.query.len){
    lapply(X=rsid.list, FUN=ensemblQueryR::ensemblQueryLDwithSNPwindow, r2=r2, d.prime=d.prime, window.size=window.size, pop=pop) %>%
      do.call("rbind", .) %>%
      return()
  } else{
    print(paste("Query length=", length(rsid.list), "which is longer than the max. query length ensembl REST API accepts. Please use ensemblQueryLDwithLargeSNPdf() instead."))
  }

}

# employs ensemblQueryLDwithSNPwindowList to make this work with an entire results dataframe,
# splits a df by phenotype, queries each one (splitting if rsid lists are >1000), then returns a df for each phenotype
# final output is a df of all input rsid, with all SNPs in LD with those, with accompanying phenotype labels
# in.table needs minimum cols: phenotype, rsID
#'
#' @param in.table data.frame with minimum columns `rsid`.
#' @param r2 Float. Measure of LD. If r-squared is provided only return pairs of variants whose r-squared value is equal to or greater than the value provided.
#' @param d.prime Float. Measure of LD. If D' is provided only return pairs of variants whose D' value is equal to or greater than the value provided.
#' @param window.size Integer. Window size in kb. The maximum allowed value for the window size is 500 kb. LD is computed for the given variant and all variants that are located within the specified window.
#' @param pop String. Population for which to compute LD. Use `ensemblQueryGetPops()` to retrieve a list of all populations with LD data. Default is 1000GENOMES:phase_3:EUR.
#'
#' @return
#'
#' @import purrr
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @export
#'
#' @examples
#'require(magrittr)
#'data.frame(rsid=rep(c("rs7153434","rs1963154","rs12672022","rs3852802","rs12324408","rs56346870"), 500)) %>%
#'    ensemblQueryLDwithLargeSNPdf(in.table=.,
#'                                 r2=0.8,
#'                                 d.prime=0.8,
#'                                 window.size=500,
#'                                 pop="1000GENOMES:phase_3:EUR")
ensemblQueryLDwithSNPwindowDataframe = function(in.table, r2=0.8, d.prime=0.8, window.size=500, pop="1000GENOMES:phase_3:EUR"){

  # library(dplyr)
  # library(magrittr)
  # library(purrr)

  # max query length ensembl REST API will accept
  max.query.len=1000

  # test for length of query
  # if query is longer than 999, add id every 500 rows and split into multiple queries
  if(length(in.table$rsid)>=max.query.len){

    print(paste("Query consists of", length(in.table$rsid), "rsid, so splitting query and running in chunks..."))

    in.table = in.table %>%
      dplyr::mutate(id = rep(seq(n()), each = 500, length = n())) %>%
      dplyr::group_by(id) %>%
      dplyr::group_split()

    # define empty vector to hold query res
    out.list = vector("list", length(in.table))

    # the queries then have to be run sequentially, not in parallel/or using lapply as the REST API will only handle sequential queries
    for(i in 1:length(in.table)){

      print(paste("Running query number", i, "of", length(in.table)))

      # run query
      if( length(in.table[[i]] %>% dplyr::pull(rsid)) <= max.query.len ){

        out.list[[i]] = in.table[[i]] %>%
          dplyr::pull(rsid) %>%
          ensemblQueryR::ensemblQueryLDwithSNPwindowList(rsid.list=., r2=r2, d.prime=d.prime, window.size=window.size, pop=pop)

      }
    }

    # bind res together and return
    print("Returning results table... ")
    out.list %>%
      lapply(X=., FUN=function(x){x %>% dplyr::mutate_all(as.character)}) %>%
      purrr::discard(., ~nrow(.) == 0) %>%
      do.call("rbind", .) %>%
      dplyr::mutate(r2 = as.numeric(r2),
                    d_prime = as.numeric(d_prime)) %>%
      tibble::tibble() %>%
      return(.)

    # otherwise, run query as usual
  } else{

    print(paste("Query consists of", length(in.table$rsid), "rsid, so running a single query..."))

    in.table %>%
      dplyr::pull(rsid) %>%
      ensemblQueryR::ensemblQueryLDwithSNPwindowList(rsid.list=., r2=r2, d.prime=d.prime, window.size=window.size, pop=pop) %>%
      dplyr::mutate(r2 = as.numeric(r2),
                    d_prime = as.numeric(d_prime)) %>%
      tibble::tibble() %>%
      return(.)
  }}

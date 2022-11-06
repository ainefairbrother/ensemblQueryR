#' Function to query Ensembl LD data with a single rsID
#'
#' @param rsid
#' @param r2
#' @param d.prime
#' @param window.size
#' @param pop
#'
#' @return data.frame with 5 columns:
#' `query` (original query SNP),
#' `snp_in_ld` (all SNPs in LD within parameters specified),
#' `r2` (the r-squared value for query-snp_in_ld),
#' `d_prime` (the D' value for query-snp_in_ld),
#' `population_name` (the population LD was calculated from)
#' @export
#'
#' @examples
#' ensemblQueryLDwithSNP(rsid="rs3851179", r2=0.8, d.prime=0.8, window.size=500, pop="1000GENOMES:phase_3:EUR")
ensemblQueryLDwithSNP = function(rsid, r2=0.8, d.prime=0.8, window.size=500, pop="1000GENOMES:phase_3:EUR"){

  # load libs
  require(httr)
  require(xml2)
  if( !("tidyverse" %in% (.packages())) ){
    require(jsonlite)
  }
  require(dplyr)
  require(tidyr)
  require(purrr)
  require(vroom)
  require(magrittr)

  #--------------------------------- get pops --------------------------------

  # server <- "https://rest.ensembl.org"
  # ext <- "/info/variation/populations/homo_sapiens?filter=LD"
  #
  # r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  #
  # jsonlite::fromJSON(jsonlite::toJSON(content(r))) %>% data.frame()

  #--------------------------------- run query -------------------------------

  server <- "https://rest.ensembl.org"
  ext <- paste0("/ld/human/",rsid,"/",pop,"?d_prime=",d.prime,";window_size=",window.size,";r2=",r2)

  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))

  stop_for_status(r)

  # use this if you get a simple nested list back, otherwise inspect its structure
  res.temp = jsonlite::fromJSON(jsonlite::toJSON(content(r))) %>%
    data.frame()

  # deal with null search result by testing for df length
  # then returning an empty df if no result
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

    res.temp %>%
      data.frame() %>%
      dplyr::arrange(r2) %>%
      dplyr::rename(query=variation1, snp_in_ld=variation2) %>%
      dplyr::relocate(query, snp_in_ld, r2, d_prime, population_name) %>%
      return()

  }
}

#' Function to apply `ensemblQueryLDwithSNP` to a list of SNPs of length<1000
#' Outputs a table containing all SNPs in LD with all query SNPs in rsid.list.
#' Outputs NA row if no SNPs in LD found.
#'
#' @param rsid.list
#' @param r2
#' @param d.prime
#' @param window.size
#' @param pop
#'
#' @return
#' @export
#'
#' @examples
#'
#' ensemblQueryLDwithSNPlist(rsid.list=c("rs7153434","rs1963154","rs12672022","rs3852802","rs12324408","rs56346870"),
#'                           r2=0.8,
#'                           d.prime=0.8,
#'                           window.size=500,
#'                           pop="1000GENOMES:phase_3:EUR")
#'
ensemblQueryLDwithSNPlist = function(rsid.list, r2=0.8, d.prime=0.8, window.size=500, pop="1000GENOMES:phase_3:EUR"){

  library(magrittr)

  # max query length ensembl REST API will accept
  max.query.len=1000

  if(length(rsid.list)<=max.query.len){
    lapply(X=rsid.list, FUN=ensemblQueryR::ensemblQueryLDwithSNP, r2=r2, d.prime=d.prime, window.size=window.size, pop=pop) %>%
      do.call("rbind", .) %>%
      return()
  } else{
    print(paste("Query length=", length(rsid.list), "which is longer than the max. query length ensembl REST API accepts. Please use ensemblQueryLDwithLargeSNPdf() instead."))
  }

}

# employs ensemblQueryLDwithSNPlist to make this work with an entire results dataframe,
# splits a df by phenotype, queries each one (splitting if rsid lists are >1000), then returns a df for each phenotype
# final output is a df of all input rsid, with all SNPs in LD with those, with accompanying phenotype labels
# in.table needs minimum cols: phenotype, rsID
#'
#' @param in.table
#' @param r2
#' @param d.prime
#' @param window.size
#' @param pop
#'
#' @return
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
ensemblQueryLDwithLargeSNPdf = function(in.table, r2=0.8, d.prime=0.8, window.size=500, pop="1000GENOMES:phase_3:EUR"){

  library(dplyr)
  library(magrittr)
  library(purrr)

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
          ensemblQueryR::ensemblQueryLDwithSNPlist(rsid.list=., r2=r2, d.prime=d.prime, window.size=window.size, pop=pop)

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
      ensemblQueryR::ensemblQueryLDwithSNPlist(rsid.list=., r2=r2, d.prime=d.prime, window.size=window.size, pop=pop) %>%
      dplyr::mutate(r2 = as.numeric(r2),
                    d_prime = as.numeric(d_prime)) %>%
      tibble::tibble() %>%
      return(.)
  }}

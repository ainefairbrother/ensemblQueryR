#' Function to query Ensembl LD data with a pair of rsIDs.
#' This function will return r-squared and D' values for the rsID pair.
#'
#' @param rsid1 String. Variant ID 1.
#' @param rsid2 String. Variant ID 2.
#' @param pop String. Population for which to compute LD. Use `ensemblQueryGetPops()` to retrieve a list of all populations with LD data. Default is 1000GENOMES:phase_3:EUR.
#'
#' @return A dataframe.
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
#' ensemblQueryLDwithSNPpair(
#'   rsid1="rs6792369",
#'   rsid2="rs1042779",
#'   pop="1000GENOMES:phase_3:EUR"
#' )
#'
ensemblQueryLDwithSNPpair = function(rsid1, rsid2, pop="1000GENOMES:phase_3:EUR"){

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

  #------------------------------ check inputs -------------------------------

  stopifnot(is.character(rsid1))
  stopifnot(is.character(rsid2))
  stopifnot(is.character(pop))

  #--------------------------------- run query -------------------------------

  server <- "https://rest.ensembl.org"
  ext <- paste0("/ld/human/pairwise/",rsid1,"/",rsid2,"?population_name=",pop)

  r <- httr::GET(url=paste(server, ext, sep = ""), content_type("application/json"))

  # stop_for_status(r)

  #-------------------- check output and write out ---------------------------

  # error handling, if 400 error, set res.temp as NA
  if(r$status_code == 400){
    print(paste0("Error 400 thrown by httr::GET. One or both of rsid1 (",rsid1,") or rsid2 (", rsid2,")",
                 " may be invalid variant rsID(s). You can check using dbSNP: https://www.ncbi.nlm.nih.gov/snp/."))
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
      print(paste0("No data found for variants ", rsid1, " and ", rsid2))
      data.frame(rep(NA, 5), row.names = c("variation1", "variation2", "d_prime", "population_name", "r2")) %>%
        t() %>%
        `rownames<-`(NULL) %>%
        as.data.frame() %>%
        dplyr::rename(query1=variation1, query2=variation2) %>%
        dplyr::relocate(query1, query2, r2, d_prime, population_name) %>%
        tibble::tibble() %>%
        return()
    } else{
      # if not 0-row (empty) df, then deal with it normally, format and prepare for return
      res.temp %>%
        data.frame() %>%
        dplyr::arrange(r2) %>%
        dplyr::rename(query1=variation1, query2=variation2) %>%
        dplyr::relocate(query1, query2, r2, d_prime, population_name) %>%
        tibble::tibble() %>%
        return()
    }
    # deal with NA search result (result of 400 error) by testing if res.temp is NA
    # then returning an empty df of the same structure
  } else{
    if(is.na(res.temp)){
      data.frame(rep(NA, 5), row.names = c("variation1", "variation2", "d_prime", "population_name", "r2")) %>%
        t() %>%
        `rownames<-`(NULL) %>%
        as.data.frame() %>%
        dplyr::rename(query1=variation1, query2=variation2) %>%
        dplyr::relocate(query1, query2, r2, d_prime, population_name) %>%
        tibble::tibble() %>%
        return()
    }
  }
}

#' `ensemblQueryLDwithSNPpairDataframe` applies `ensemblQueryLDwithSNPpair` to a data.frame of rsID pairs
#'
#' @param in.table data.frame containing SNP pairs. Columns must include `rsid1` for the first member of the pair and `rsid2` for the second member of the pair.
#' @param pop String. Population for which to compute LD. Use `ensemblQueryGetPops()` to retrieve a list of all populations with LD data. Default is 1000GENOMES:phase_3:EUR.
#' @param cores Integer. A value between 1 and 10 is accepted, as this prevents the server returning overload-related errors.
#'
#' @return A dataframe.
#' @export
#'
#' @examples
#'ensemblQueryLDwithSNPpairDataframe(
#'  in.table=data.frame(rsid1=rep("rs6792369", 10), rsid2=rep("rs1042779", 10)),
#'  pop="1000GENOMES:phase_3:EUR")
#'
ensemblQueryLDwithSNPpairDataframe = function(in.table, pop="1000GENOMES:phase_3:EUR", cores=1){ #keep.original.table.row.n=FALSE

  #------------------------------ test -------------------------------
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
  # in.table=data.frame(rsid1=rep("rs6792369", 10), rsid2=rep("rs1042779", 10))
  # pop="1000GENOMES:phase_3:EUR"
  # keep.original.table.row.n=FALSE
  # cores=2

  #------------------------------ check inputs -------------------------------

  stopifnot(is.data.frame(in.table))
  stopifnot(is.character(pop))
  # stopifnot(is.logical(keep.original.table.row.n))
  stopifnot(is.numeric(cores))

  # check that the cores arg is set above 0 but not above the max. available
  # cores.available = parallel::detectCores()
  # stopifnot("Cores must be a value between 0 and 10"= (cores>0) & (cores<=10))

  #--------------------------------- main ------------------------------------

  if( is.data.frame(in.table)==TRUE ){
    if( (("rsid1" %in% colnames(in.table)) & ("rsid2" %in% colnames(in.table))) ){

      print(paste("Running ensemblQueryLDwithSNPpairDataframe to retrieve LD metrics for", nrow(in.table), "variant pairs..."))

      if(cores>1){
        print(paste(
          "Parallelising query using", cores, "cores"
        ))
      }

      res = parallel::mclapply(X=c(1:nrow(in.table)), mc.cores=cores, FUN=function(x){

        ensemblQueryLDwithSNPpair(rsid1=in.table$rsid1[x],
                                  rsid2=in.table$rsid2[x],
                                  pop=pop) %>%
          tidyr::unnest(cols = c(query1, query2, r2, d_prime, population_name)) %>%
          as.data.frame()

      }) %>%
        do.call("rbind", .) %>%
        tibble::tibble() %>%
        return()

      # # either filter null rows, or keep depending on arg - this can clean up rows where no data was found for the snp pair
      # if(keep.original.table.row.n==FALSE){
      #   res.original.len = nrow(res)
      #
      #   res = res %>%
      #     tidyr::drop_na()
      #
      #   res.filtered.len = nrow(res)
      #
      #   if(res.original.len!=res.filtered.len){
      #     n.filtered = res.original.len-res.filtered.len
      #     print(paste0(n.filtered, " rows filtered out due to no data for rsID pairs."))
      #   }
      #
      #   return(res)
      #
      # } else{
      #   return(res)
      # }

    } else{
      print("Error: columns rsid1 and rsid2 do not exist in in.table.")
      stop()
    }
  } else{
    print("Error: in.table is not a data.frame.")
    stop()
  }
}




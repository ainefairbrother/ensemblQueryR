#' Function to query Ensembl LD data with a single rsID
#'
#' @param rsid String. Variant ID.
#' @param r2 Float. Measure of LD. If r-squared is provided only return pairs of variants whose r-squared value is equal to or greater than the value provided.
#' @param d.prime Float. Measure of LD. If D' is provided only return pairs of variants whose D' value is equal to or greater than the value provided.
#' @param window.size Integer. Window size in kb. The maximum allowed value for the window size is 500 kb. LD is computed for the given variant and all variants that are located within the specified window.
#' @param pop String. Population for which to compute LD. Use `ensemblQueryGetPops()` to retrieve a list of all populations with LD data. Default is 1000GENOMES:phase_3:EUR.
#'
#'
#' @return A dataframe.
#'
#' @import httr
#' @import xml2
#' @importFrom jsonlite fromJSON toJSON
#' @import dplyr
#' @import tidyr
#' @import vroom
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' ensemblQueryLDwithSNPwindow(rsid="rs3851179", r2=0.8, d.prime=0.8,
#'                             window.size=500, pop="1000GENOMES:phase_3:EUR")
#'
ensemblQueryLDwithSNPwindow = function(rsid, r2=0.8, d.prime=0.8, window.size=500, pop="1000GENOMES:phase_3:EUR"){

  # TEST
  # load libs
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

  #------------------------------ check inputs -------------------------------

  stopifnot(is.character(rsid))
  stopifnot(is.character(r2) | is.numeric(r2))
  stopifnot(is.character(d.prime) | is.numeric(d.prime))
  stopifnot(is.character(window.size) | is.numeric(window.size))
  stopifnot(is.character(pop))

  #--------------------------------- run query -------------------------------

  server <- "https://rest.ensembl.org"
  ext <- paste0("/ld/human/",rsid,"/",pop,"?d_prime=",d.prime,";window_size=",window.size,";r2=",r2)

  r <- httr::GET(url=paste(server, ext, sep = ""), httr::content_type("application/json"))

  #-------------------- check output and write out ---------------------------

  # stop_for_status(r)

  # error handling, if 400 error, set res.temp as NA
  if(r$status_code == 400){
    message("Error 400 thrown by httr::GET. This may not be a valid SNP rsID, check using dbSNP: https://www.ncbi.nlm.nih.gov/snp/.")
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
        tibble::tibble() %>%
        return()
    } else{
      # if not 0-row (empty) df, then deal with it normally, format and prepare for return
      res.temp %>%
        data.frame() %>%
        dplyr::arrange(r2) %>%
        dplyr::rename(query=variation1, snp_in_ld=variation2) %>%
        dplyr::relocate(query, snp_in_ld, r2, d_prime, population_name) %>%
        tibble::tibble() %>%
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
        tibble::tibble() %>%
        return()
    }
  }
}

#' `ensemblQueryLDwithSNPwindowDataframe` applies `ensemblQueryLDwithSNPwindow` to a data.frame of rsIDs.
#'
#' @param in.table data.frame containing SNP pairs. Columns must include `rsid1` for the first member of the pair and `rsid2` for the second member of the pair.
#' @param r2 Float. Measure of LD. If r-squared is provided only return pairs of variants whose r-squared value is equal to or greater than the value provided.
#' @param d.prime Float. Measure of LD. If D' is provided only return pairs of variants whose D' value is equal to or greater than the value provided.
#' @param window.size Integer. Window size in kb. The maximum allowed value for the window size is 500 kb. LD is computed for the given variant and all variants that are located within the specified window.
#' @param pop String. Population for which to compute LD. Use `ensemblQueryGetPops()` to retrieve a list of all populations with LD data. Default is 1000GENOMES:phase_3:EUR.
#' @param cores Integer. A value between 1 and 10 is accepted, as this prevents the server returning overload-related errors.
#'
#' @return A dataframe.
#'
#' @import purrr
#' @import parallel
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' in.table = data.frame(rsid = rep(c("rs7153434","rs1963154","rs12672022",
#'                                    "rs3852802","rs12324408","rs56346870"), 10))
#'
#' ensemblQueryLDwithSNPwindowDataframe(in.table=in.table,
#'                                      r2=0.8,
#'                                      d.prime=0.8,
#'                                      window.size=500,
#'                                      pop="1000GENOMES:phase_3:EUR",
#'                                      cores=2)
#'
ensemblQueryLDwithSNPwindowDataframe = function(in.table, r2=0.8, d.prime=0.8, window.size=500, pop="1000GENOMES:phase_3:EUR", cores=1){

  #------------------------------ test -------------------------------
  # library(purrr)
  # library(parallel)
  # library(magrittr)
  #
  # in.table = data.frame(rsid = rep(c("rs7153434","rs1963154","rs12672022","rs3852802","rs12324408","rs56346870"), 10))
  # r2=0.8
  # d.prime=0.8
  # window.size=500
  # pop="1000GENOMES:phase_3:EUR"
  # cores=2

  #------------------------------ check inputs -------------------------------

  stopifnot(is.data.frame(in.table))
  stopifnot(is.numeric(r2) | is.character(r2))
  stopifnot(is.numeric(d.prime) | is.character(d.prime))
  stopifnot(is.numeric(window.size) | is.character(window.size))
  stopifnot(is.character(pop))
  stopifnot(is.numeric(cores))

  #--------------------------------- main ------------------------------------

  if( is.data.frame(in.table)==TRUE ){
    if( ("rsid" %in% colnames(in.table)) ){

      message(paste("Running ensemblQueryLDwithSNPwindowDataframe to retrieve LD metrics for" , nrow(in.table), "central variants, where: \n",
                    paste0("The window around each central variant is ", window.size, ","),
                    paste0("r-squared = ",r2, ","),
                    paste0("D' = ", d.prime, ".")
      ))

      # check system
      sys = Sys.info()['sysname'] %>% grepl("Windows",.)

      if(sys==FALSE){

        # else if cores set to more than 1, parallelise
        if(cores>1){
          message(paste(
            "Parallelising query using", cores, "cores."
          ))
        }

        parallel::mclapply(X=c(1:nrow(in.table)), mc.cores=cores, FUN=function(x){
          ensemblQueryLDwithSNPwindow(rsid=in.table$rsid[x],
                                      r2=r2,
                                      d.prime=d.prime,
                                      window.size=window.size,
                                      pop=pop) %>%
            tidyr::unnest(cols = c(query, snp_in_ld, r2, d_prime, population_name)) %>%
            as.data.frame()

        }) %>%
          do.call("rbind", .) %>%
          tibble::tibble() %>%
          return(.)

      }else{

        message("Windows OS detected. Cannot run parallel queries using parallel::mclapply. Using lapply instead.")

        lapply(X=c(1:nrow(in.table)), FUN=function(x){
          ensemblQueryLDwithSNPwindow(rsid=in.table$rsid[x],
                                      r2=r2,
                                      d.prime=d.prime,
                                      window.size=window.size,
                                      pop=pop) %>%
            tidyr::unnest(cols = c(query, snp_in_ld, r2, d_prime, population_name)) %>%
            as.data.frame()

        }) %>%
          do.call("rbind", .) %>%
          tibble::tibble() %>%
          return(.)

      }

    } else{
      message("Error: column rsid does not exist in in.table.")
      stop()
    }
  } else{
    message("Error: in.table is not a data.frame.")
    stop()
  }
}








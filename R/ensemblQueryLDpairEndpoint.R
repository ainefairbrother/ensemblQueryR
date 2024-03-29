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
#' @importFrom jsonlite fromJSON toJSON
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
  #
  # rsid1="rs6792369"
  # rsid2="rs1042779"
  # pop="1000GENOMES:phase_3:EUR"

  #------------------------------ check inputs -------------------------------

  stopifnot(is.character(rsid1))
  stopifnot(is.character(rsid2))
  stopifnot(is.character(pop))

  #--------------------------------- run query -------------------------------

  server <- "https://rest.ensembl.org"
  ext <- paste0("/ld/human/pairwise/",rsid1,"/",rsid2,"?population_name=",pop)

  r <- httr::GET(url=paste(server, ext, sep = ""), httr::content_type("application/json"))

  # stop_for_status(r)

  #-------------------- check output and write out ---------------------------

  # error handling, if 400 error, set res.temp as NA
  if(r$status_code == 400){
    message(paste0("Error 400 thrown by httr::GET. One or both of rsid1 (",rsid1,") or rsid2 (", rsid2,")",
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
      message(paste0("No data found for variants ", rsid1, " and ", rsid2))
      data.frame(rep(NA, 5), row.names = c("variation1", "variation2", "d_prime", "population_name", "r2")) %>%
        t() %>%
        `rownames<-`(NULL) %>%
        as.data.frame() %>%
        dplyr::rename(query1=variation1, query2=variation2) %>%
        tidyr::unnest(cols=c(query1, query2, r2, d_prime, population_name)) %>%
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
        tidyr::unnest(cols=c(query1, query2, r2, d_prime, population_name)) %>%
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
        tidyr::unnest(cols=c(query1, query2, r2, d_prime, population_name)) %>%
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
#'  in.table=data.frame(rsid1=rep("rs6792369", 10),
#'                      rsid2=rep("rs1042779", 10)),
#'                      pop="1000GENOMES:phase_3:EUR"
#')
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
  # cores=1

  #------------------------------ check inputs -------------------------------

  stopifnot(is.data.frame(in.table))
  stopifnot(is.character(pop))
  stopifnot(is.numeric(cores))

  #--------------------------------- main ------------------------------------

  if( is.data.frame(in.table)==TRUE ){
    if( (("rsid1" %in% colnames(in.table)) & ("rsid2" %in% colnames(in.table))) ){

      message(paste("Running ensemblQueryLDwithSNPpairDataframe to retrieve LD metrics for", nrow(in.table), "variant pairs..."))

      # check system
      sys = Sys.info()['sysname'] %>% grepl("Windows",.)

      if(sys==FALSE){

        if(cores>1){
          message(paste(
            "Parallelising query using", cores, "cores"
          ))
        }

        parallel::mclapply(X=c(1:nrow(in.table)), mc.cores=cores, FUN=function(x){

          ensemblQueryLDwithSNPpair(rsid1=in.table$rsid1[x],
                                    rsid2=in.table$rsid2[x],
                                    pop=pop) %>%
            tidyr::unnest(cols = c(query1, query2, r2, d_prime, population_name)) %>%
            as.data.frame()

        }) %>%
          do.call("rbind", .) %>%
          tibble::tibble() %>%
          return()

      } else{
        if(sys==TRUE){

          message("Windows OS detected. Cannot run parallel queries using parallel::mclapply. Using lapply instead.")

          lapply(X=c(1:nrow(in.table)), FUN=function(x){

            ensemblQueryLDwithSNPpair(rsid1=in.table$rsid1[x],
                                      rsid2=in.table$rsid2[x],
                                      pop=pop) %>%
              tidyr::unnest(cols = c(query1, query2, r2, d_prime, population_name)) %>%
              as.data.frame()

          }) %>%
            do.call("rbind", .) %>%
            tibble::tibble() %>%
            return()
        }
      }

    } else{
      message("Error: columns rsid1 and rsid2 do not exist in in.table.")
      stop()
    }
  } else{
    message("Error: in.table is not a data.frame.")
    stop()
  }
}




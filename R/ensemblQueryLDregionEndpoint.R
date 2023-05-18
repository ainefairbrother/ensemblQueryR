#' Function to query Ensembl LD data within a genomic window to get all variant pairs in the specified region and associated LD metrics.
#'
#' @param chr String. Chromosome that the query region is located on.
#' @param start String. Base pair that the query region starts at.
#' @param end String. Base pair that the query region ends at.
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
  #
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

  r <- httr::GET(url=paste(server, ext, sep = ""), httr::content_type("application/json"))

  # stop_for_status(r)

  #-------------------- check output and write out ---------------------------

  # error handling, if 400 error, set res.temp as NA
  if(r$status_code == 400){
    message(paste(
      "Error 400 thrown by httr::GET. Please check your input data for validity.",
      "CHR=",chr,"START=",start,"END=",end, "",
      chr, "may be an invalid chromosome identifier - this should be 1,2,3 ... 22, X or Y.",
      start, "may be an invalid genomic position - this should be a number lower than the maximum number of bases in the query chromosome.",
      end, "may be an invalid genomic position - this should be a number lower than the maximum number of bases in the query chromosome."
    ))
    res.temp = NA
  } else{
    # if no error, use this if you get a simple nested list back, otherwise inspect its structure
    res.temp = jsonlite::fromJSON(jsonlite::toJSON(httr::content(r))) %>%
      data.frame()
  }

  # deal with null search result by testing for df length
  # then returning an empty df if no result
  if(is.data.frame(res.temp)){
    if(nrow(res.temp)==0){

      message(paste0("No LD data found for the genomic coordinates input: ", chr,":",start,"-",end))

      data.frame(rep(NA, 5), row.names = c("variation1", "variation2", "d_prime", "population_name", "r2")) %>%
        t() %>%
        `rownames<-`(NULL) %>%
        as.data.frame() %>%
        dplyr::rename(rsid1=variation1, rsid2=variation2) %>%
        dplyr::mutate(query_chr=chr, query_start=start, query_end=end) %>%
        dplyr::relocate(query_chr, query_start, query_end, rsid1, rsid2, r2, d_prime, population_name) %>%
        tibble::tibble() %>%
        return()
    } else{
      # if not 0-row (empty) df, then deal with it normally, format and prepare for return
      res.temp %>%
        data.frame() %>%
        dplyr::arrange(r2) %>%
        dplyr::rename(rsid1=variation1, rsid2=variation2) %>%
        dplyr::mutate(query_chr=chr, query_start=start, query_end=end) %>%
        dplyr::relocate(query_chr, query_start, query_end, rsid1, rsid2, r2, d_prime, population_name) %>%
        tibble::tibble() %>%
        return()
    }
    # deal with NA search result (result of 400 error) by testing if res.temp is NA
    # then returning an empty df of the same structure
  } else{
    if(is.na(res.temp)){

      message(paste0("No LD data found for the genomic coordinates input: ", chr,":",start,"-",end))

      data.frame(rep(NA, 5), row.names = c("variation1", "variation2", "d_prime", "population_name", "r2")) %>%
        t() %>%
        `rownames<-`(NULL) %>%
        as.data.frame() %>%
        dplyr::rename(rsid1=variation1, rsid2=variation2) %>%
        dplyr::mutate(query_chr=chr, query_start=start, query_end=end) %>%
        dplyr::relocate(query_chr, query_start, query_end, rsid1, rsid2, r2, d_prime, population_name) %>%
        tibble::tibble() %>%
        return()
    }
  }
}

#' `ensemblQueryLDwithSNPregionDataframe` applies `ensemblQueryLDwithSNPregion` to a data.frame of genomic coordinates and returns all variant pairs present in each specified genomic region and their associated LD metrics.
#'
#' @param in.table Dataframe containing genomic coordinates. Columns must include `chr` (the chromosome), `start` (the starting genomic coordinate) and `end` (the ending genomic coordinate).
#' @param pop String. Population for which to compute LD. Use `ensemblQueryGetPops()` to retrieve a list of all populations with LD data. Default is 1000GENOMES:phase_3:EUR.
#' @param cores Integer. A value between 1 and 10 is accepted, as this prevents the server returning overload-related errors.
#'
#' @return A dataframe.
#' @export
#'
#' @examples
#' library(magrittr)
#'
#' data.frame(
#'   chr=rep(c("6"), 5),
#'   start=rep(c("25837556"), 5),
#'   end=rep(c("25943455"), 5)
#' ) %>%
#'   ensemblQueryLDwithSNPregionDataframe(
#'     in.table=.,
#'     pop="1000GENOMES:phase_3:EUR",
#'     cores = 1
#'   )
#'
ensemblQueryLDwithSNPregionDataframe = function(in.table, pop="1000GENOMES:phase_3:EUR", cores=1){

  # # Test
  # library(purrr)
  # library(parallel)
  # library(magrittr)

  #------------------------------ check inputs -------------------------------

  stopifnot(is.data.frame(in.table))
  stopifnot(is.character(pop))
  stopifnot(is.numeric(cores))

  #--------------------------------- main ------------------------------------

  if( is.data.frame(in.table)==TRUE ){
    if( (("chr" %in% colnames(in.table)) & ("start" %in% colnames(in.table)) & ("end" %in% colnames(in.table))) ){

      message(paste("Running ensemblQueryLDwithSNPregionDataframe to retrieve LD metrics for variant pairs in", nrow(in.table), "genomic regions"))

      # check system
      sys = Sys.info()['sysname'] %>% grepl("Windows",.)

      if(sys==FALSE){

        if(cores>1){
          message(paste(
            "Parallelising query using", cores, "cores."
          ))
        }

        parallel::mclapply(X=c(1:nrow(in.table)), mc.cores=cores, FUN=function(x){

          ensemblQueryLDwithSNPregion(chr=in.table$chr[x],
                                      start=in.table$start[x],
                                      end=in.table$end[x],
                                      pop=pop) %>%
            tidyr::unnest(cols = c(query_chr, query_start, query_end, rsid1, rsid2, r2, d_prime, population_name)) %>%
            dplyr::relocate(query_chr, query_start, query_end, rsid1, rsid2, r2, d_prime, population_name) %>%
            as.data.frame()

        }) %>%
          do.call("rbind", .) %>%
          tibble::tibble() %>%
          return(.)
      } else{
        lapply(X=c(1:nrow(in.table)), FUN=function(x){

          message("Windows OS detected. Cannot run parallel queries using parallel::mclapply. Using lapply instead.")

          ensemblQueryLDwithSNPregion(chr=in.table$chr[x],
                                      start=in.table$start[x],
                                      end=in.table$end[x],
                                      pop=pop) %>%
            tidyr::unnest(cols = c(query_chr, query_start, query_end, rsid1, rsid2, r2, d_prime, population_name)) %>%
            dplyr::relocate(query_chr, query_start, query_end, rsid1, rsid2, r2, d_prime, population_name) %>%
            as.data.frame()

        }) %>%
          do.call("rbind", .) %>%
          tibble::tibble() %>%
          return(.)
      }

    } else{
      message("Error: one or more of 'chr', 'start' or 'end' column(s) does not exist in in.table.")
      stop()
    }
  } else{
    message("Error: in.table is not a data.frame.")
    stop()
  }
}

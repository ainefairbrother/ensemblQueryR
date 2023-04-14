# Computes and returns LD values between the given variants.

ensemblQueryLDwithSNPpair = function(rsid1, rsid2, pop="1000GENOMES:phase_3:EUR"){

  # TEST
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

  # rsid1="rs6792369"
  # rsid2="rs1042779"
  # pop="1000GENOMES:phase_3:EUR"

  #--------------------------------- run query -------------------------------

  server <- "https://rest.ensembl.org"
  ext <- paste0("/ld/human/pairwise/",rsid1,"/",rsid2,"?population_name=",pop)

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
      data.frame(rep(NA, 5), row.names = c("variation1", "variation2", "d_prime", "population_name", "r2")) %>%
        t() %>%
        `rownames<-`(NULL) %>%
        as.data.frame() %>%
        dplyr::rename(query1=variation1, query2=variation2) %>%
        dplyr::relocate(query1, query2, r2, d_prime, population_name) %>%
        return()
    } else{
      # if not 0-row (empty) df, then deal with it normally, format and prepare for return
      res.temp %>%
        data.frame() %>%
        dplyr::arrange(r2) %>%
        dplyr::rename(query1=variation1, query2=variation2) %>%
        dplyr::relocate(query1, query2, r2, d_prime, population_name) %>%
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
        return()
    }
  }
}

test.df = data.frame(rsid1=rep("rs6792369", 10), rsid2=rep("rs1042779", 10))

ensemblQueryLDwithLargeSNPdf = function(in.table, r2=0.8, d.prime=0.8, window.size=500, pop="1000GENOMES:phase_3:EUR"){

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


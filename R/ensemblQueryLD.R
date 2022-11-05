



# main function that queries with one rsid, outputs all SNPs in LD with the query
query.ensemble.ld.with.snp = function(rsid, r2=0.8){

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

  #--------------------------------- get pops --------------------------------

  # server <- "https://rest.ensembl.org"
  # ext <- "/info/variation/populations/homo_sapiens?filter=LD"
  #
  # r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  #
  # jsonlite::fromJSON(jsonlite::toJSON(content(r))) %>% data.frame()

  #--------------------------------- run query -------------------------------

  server <- "https://rest.ensembl.org"
  ext <- paste0("/ld/human/",rsid,"/1000GENOMES:phase_3:EUR?d_prime=0.8;window_size=500;r2=",r2)

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
      dplyr::rename(query=variation1, snp.in.ld=variation2) %>%
      dplyr::relocate(query, snp.in.ld, r2, d_prime, population_name) %>%
      dplyr::mutate(query=rsid) %>%
      return()

  } else{

    res.temp %>%
      data.frame() %>%
      dplyr::arrange(r2) %>%
      dplyr::rename(query=variation1, snp.in.ld=variation2) %>%
      dplyr::relocate(query, snp.in.ld, r2, d_prime, population_name) %>%
      return()

  }
}

# function that runs lapply, applying query.ensemble.ld.with.snp to a list of SNPs,
# outputs a table containing all SNPs in LD with all query SNPs in rsid.list, or NA row if no SNPs in LD found
apply.query.ensemble.ld.with.snp.to.list = function(rsid.list, r2=0.8){

  lapply(X=rsid.list, FUN=query.ensemble.ld.with.snp, r2=r2) %>%
    do.call("rbind", .) %>%
    return()

}

# employs apply.query.ensemble.ld.with.snp.to.list to make this work with an entire results dataframe,
# splits a df by phenotype, queries each one (splitting if rsid lists are >1000), then returns a df for each phenotype
# final output is a df of all input rsid, with all SNPs in LD with those, with accompanying phenotype labels
# in.table needs minimum cols: phenotype, rsID
apply.query.ensemble.ld.with.snp.to.df = function(in.table){

  # # TEST
  # in.table = meta %>%
  #   dplyr::filter(diagnosis=="Control") %>%
  #   dplyr::distinct()

  exp.res.split.by.pheno = in.table %>%
    dplyr::select(diagnosis, phenotype, rsId) %>%
    tidyr::drop_na(rsId) %>%
    split(., f=.$phenotype)

  out.list = vector("list", length(exp.res.split.by.pheno))
  names(out.list) = names(exp.res.split.by.pheno)
  max.query.len = 1000

  for(i in 1:length(exp.res.split.by.pheno)){

    print(paste("Running query number", i, "of", length(exp.res.split.by.pheno)))

    # run the query normally for shorter queries
    if( length(exp.res.split.by.pheno[[i]] %>% dplyr::pull(rsId)) <= max.query.len ){

      out.list[[i]] = exp.res.split.by.pheno[[i]] %>%
        dplyr::pull(rsId) %>%
        apply.query.ensemble.ld.with.snp.to.list(rsid.list=.)

      # split the query if longer than 1000, store the split queries in temp.out.list, then stick back together and store in out.list
      # because the ensembl API doesn't accept n>1000 queries
    } else{

      temp = exp.res.split.by.pheno[[i]] %>%
        dplyr::pull(rsId) %>%
        split(., ceiling(seq_along(.)/max.query.len))

      temp.out.list = vector("list", length(temp))

      print(paste("Query too long, splitting into", length(temp), "chunks"))

      for(j in 1:length(temp)){

        print(paste("Running sub-query number", j))

        temp.out.list[[j]] = temp[[j]] %>%
          apply.query.ensemble.ld.with.snp.to.list(rsid.list=.)
      }

      out.list[[i]] = temp.out.list %>%
        purrr::discard(., ~nrow(.) == 0) %>%
        do.call("rbind", .)

    }
  }

  print("Returning results table... ")
  out.list %>%
    lapply(X=., FUN=function(x){x %>% dplyr::mutate_all(as.character)}) %>%
    dplyr::bind_rows(., .id = "phenotype") %>%
    return(.)

}




#### -------------- implement functions -------------- ####

# start_time = Sys.time()
#
#
# vroom::vroom(
#   "/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/PP_PD_mega_analysis/log10/aggregated_tables/matrixeqtl_mega_res_aggregated_FDR=0.05_filter_peaksnp_varid_justPEAKSNP_annot_991_map_alan_aminah_res.csv", guess_max = Inf, show_col_types = FALSE) %>%
#   dplyr::filter(diagnosis=="Control") %>%
#   apply.query.ensemble.ld.with.snp.to.df(in.table=.) %>%
#   vroom::vroom_write(x=., file="/home/abrowne/projects/amppd_analysis/data/ensemble_LD_out/matrixeqtl_mega_res_aggregated_FDR=0.05_filter_peaksnp_varid_justPEAKSNP_Control.csv", delim=",")
#
#
# end_time = Sys.time()
# time_taken = end_time - start_time
# print(paste(time_taken, "taken to process", num.snps, "rsIDs"))


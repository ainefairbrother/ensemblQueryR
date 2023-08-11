#' Function to get list of populations that Ensembl has available to query LD.
#'
#' @return data.frame of populations.
#'
#' @import httr
#' @import xml2
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' ensemblQueryGetPops()
#'
ensemblQueryGetPops = function(){

  #--------------------------------- get pops --------------------------------

  server = "https://rest.ensembl.org"
  ext = "/info/variation/populations/homo_sapiens?filter=LD"

  r = httr::GET(paste(server, ext, sep = ""), httr::content_type("application/json"))

  jsonlite::fromJSON(jsonlite::toJSON(httr::content(r))) %>%
    data.frame() %>%
    return()

}

#' Function to check whether the Ensembl server is up by pinging it.
#'
#' @return Integer. Where 1 is indicative of a successful ping.
#'
#' @import httr
#' @import xml2
#' @importFrom jsonlite fromJSON toJSON
#'
#' @export
#'
#' @examples
#' pingEnsembl()
#'
pingEnsembl = function(){

  server = "https://rest.ensembl.org"
  ext = "/info/ping?"
  r = httr::GET(paste(server, ext, sep = ""), content_type("application/json"))
  response = jsonlite::fromJSON(jsonlite::toJSON(content(r)))$ping

  if(response==1){
    message("Server OK.")
  } else{
    message("Server may be experiencing issues.")
  }

  return(response)
}

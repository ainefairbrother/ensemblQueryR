#' Function to get list of populations that Ensembl has available to query LD.
#'
#' @return data.frame of populations.
#'
#' @import httr
#' @import xml2
#' @import jsonlite
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
#' @import jsonlite
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
    print("Server OK.")
  }

  return(response)
}

# estimateQueriesPerHour = function(cores, n.queries){
#
#   # based on a single-core run of 54000 queries (max ensembl API query rate per hour) to ensemblQueryLDwithSNPpairDataframe,
#   # the time taken to run 54000 queries was estimated to be 1.937771.
#   # Based on this, this function takes the cores and number of queries input by the user and will output the predicted queries per hour that your run will likely spawn
#   # this assumes a linear effect of additional cores
#   cores=20
#   n.queries=1000
#   api.limit.hourly = 54000
#   time_to_run_54000_in_hours = (116.266269091765/60)
#   per_hour_query_rate = 54000/time_to_run_54000_in_hours
#   per_minute_query_rate = per_hour_query_rate/60
#   per_second_query_rate = per_minute_query_rate/60
#   time_for_one_query_seconds = (time_to_run_54000_in_hours/54000)*60*60
#
#
#
#   if(time_for_user_queries < 15){
#     print(paste(
#       "Warning: your query of size",n.queries,"using",cores,"cores may exceed the Ensembl REST API hourly query limit. Consider using fewer cores or splitting your query into smaller chunks."
#     ))
#   }
#
#
#   predicted_requests_per_hour = n.queries/(time_to_run_54000_in_hours/cores)
#
#
#
# }



# -- 1. test LDlinkR --------------------------------------------

# clear env
rm(list = ls())

# set n queries
n=10000

start_time = Sys.time()

# load libs
library(LDlinkR)
library(magrittr)

t = list(
  c(rep("rs6792369", n)),
  c(rep("rs1042779", n))
)

u = lapply(X=c(1:n),
           FUN=function(x){

             LDlinkR::LDpair(
               var1=t[[1]][x],
               var2=t[[2]][x],
               pop = "CEU",
               token = "5a141d21fa57",
               output = "table"
             ) %>%
               return()

           }) %>%
  dplyr::bind_rows()

end_time = Sys.time()
time_taken = end_time - start_time
print(paste("Time taken for LDlinkR to process", n,"rsID pairs:", time_taken))

# 43.39 minutes, then 'Bad Gateway (HTTP 502)' error

# -- 2. test ensemblQueryR --------------------------------------------

# clear env
rm(list = ls())

start_time = Sys.time()

# set n queries
n=10000

# load libs
library(ensemblQueryR)
library(magrittr)

t = data.frame(
  rsid1 = c(rep("rs6792369", n)),
  rsid2 = c(rep("rs1042779", n))
)

u = ensemblQueryLDwithSNPpairDataframe(
  in.table=t,
  pop="1000GENOMES:phase_3:EUR",
  keep.original.table.row.n=FALSE,
  parallelise=FALSE)

end_time = Sys.time()
time_taken = end_time - start_time
print(paste("Time taken for ensemblQueryR to process", n,"rsID pairs:", time_taken))

# 9.6 minutes, completed successfully











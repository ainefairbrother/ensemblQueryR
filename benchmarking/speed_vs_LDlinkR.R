
# -- 1. test LDlinkR --------------------------------------------

# clear env
rm(list = ls())

# set n queries
n=1000

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
time_taken = difftime(end_time, start_time, units='mins')
print(paste("Time taken for LDlinkR to process", n,"rsID pairs:", time_taken))

# save results
data.frame(function_tested="LDlinkR::LDpair",
           n_queries=n,
           time.min=time_taken,
           nrow.outtable=try(nrow(u))[1]) %>%
  vroom::vroom_write(x=., file=paste0("/home/abrowne/projects/ensemblQueryR/benchmarking/results/", "LDpair", ".", n, ".csv"))

# -- 2a. test ensemblQueryR::ensemblQueryLDwithSNPpair --------------------------------------------

# clear env
rm(list = ls())

start_time = Sys.time()

# set n queries
n=1000

# load libs
library(ensemblQueryR)
library(magrittr)

t = list(
  c(rep("rs6792369", n)),
  c(rep("rs1042779", n))
)

u = lapply(X=c(1:n),
           FUN=function(x){
             ensemblQueryLDwithSNPpair(
               rsid1=t[[1]][x],
               rsid2=t[[2]][x],
               pop="1000GENOMES:phase_3:EUR"
             ) %>%
               return()
           }) %>%
  dplyr::bind_rows()

end_time = Sys.time()
time_taken = difftime(end_time, start_time, units='mins')
print(paste("Time taken for ensemblQueryR to process", n,"rsID pairs:", time_taken))

# save results
data.frame(function_tested="ensemblQueryR::ensemblQueryLDwithSNPpair",
           n_queries=n,
           time.min=time_taken,
           nrow.outtable=try(nrow(u))[1]) %>%
  vroom::vroom_write(x=., file=paste0("/home/abrowne/projects/ensemblQueryR/benchmarking/results/", "ensemblQueryLDwithSNPpair", ".", n, ".csv"))

# -- 2b. test ensemblQueryR::ensemblQueryLDwithSNPpairDataframe --------------------------------------------

# clear env
rm(list = ls())

start_time = Sys.time()

# set n queries
n=1000

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
time_taken = difftime(end_time, start_time, units='mins')
print(paste("Time taken for ensemblQueryR to process", n,"rsID pairs:", time_taken))

# save results
data.frame(function_tested="ensemblQueryR::ensemblQueryLDwithSNPpairDataframe",
           n_queries=n,
           time.min=time_taken,
           nrow.outtable=try(nrow(u))[1]) %>%
  vroom::vroom_write(x=., file=paste0("/home/abrowne/projects/ensemblQueryR/benchmarking/results/", "ensemblQueryLDwithSNPpairDataframe", ".", n, ".csv"))

# -- 2c. test ensemblQueryR::ensemblQueryLDwithSNPpairDataframe with parallel --------------------------------------------

# clear env
rm(list = ls())

start_time = Sys.time()

# set n queries
n=1000

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
  parallelise=TRUE,
  n.cores=10)

end_time = Sys.time()
time_taken = difftime(end_time, start_time, units='mins')
print(paste("Time taken for ensemblQueryR to process", n,"rsID pairs:", time_taken))

# save results
data.frame(function_tested="ensemblQueryR::ensemblQueryLDwithSNPpairDataframe",
           n_queries=n,
           time.min=time_taken,
           nrow.outtable=try(nrow(u))[1]) %>%
  vroom::vroom_write(x=., file=paste0("/home/abrowne/projects/ensemblQueryR/benchmarking/results/", "ensemblQueryLDwithSNPpairDataframe.parallel10", ".", n, ".csv"))

# -- 3. Plot performance --------------------------------------------




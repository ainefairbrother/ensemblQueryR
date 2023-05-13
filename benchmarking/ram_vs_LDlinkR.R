# load libs
devtools::load_all("/home/abrowne/projects/ensemblQueryR/")
library(LDlinkR)
library(magrittr)
library(peakRAM)

# -- 1. Run for 100,1000,10000 queries --------------------------------------------

for(i in c(1:10)){
  for(n in c(100,1000,10000)){

    print(paste(
      "Processing", i, n
    ))

    # -- 1. test LDlinkR --------------------------------------------

    # clear env
    rm(t, u)

    t = list(
      c(rep("rs6792369", n)),
      c(rep("rs1042779", n))
    )

    u=try(peakRAM(
      lapply(X=c(1:n),
             FUN=function(x){
               print(x)
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
    ))

    # save results
    if( grepl("Error", u[1])==FALSE ){
      data.frame(function_tested="LDlinkR::LDpair",
                 n_queries=n,
                 time.sec=u$Elapsed_Time_sec,
                 peak.ram = u$Peak_RAM_Used_MiB,
                 total.ram = u$Total_RAM_Used_MiB) %>%
        vroom::vroom_write(x=., file=paste0("/home/abrowne/projects/ensemblQueryR/benchmarking/results/ram/", "LDpair", ".", n,".",i,".csv"))
    } else{
      print(u)
    }

    # -- 2a. test ensemblQueryR::ensemblQueryLDwithSNPpair --------------------------------------------

    # clear env
    rm(t, u)

    t = list(
      c(rep("rs6792369", n)),
      c(rep("rs1042779", n))
    )

    u=try(peakRAM(
      lapply(X=c(1:n),
             FUN=function(x){

               ensemblQueryLDwithSNPpair(
                 rsid1=t[[1]][x],
                 rsid2=t[[2]][x],
                 pop="1000GENOMES:phase_3:EUR"
               ) %>%
                 return()

             }) %>%
        dplyr::bind_rows()
    ))

    # save results
    if( grepl("Error", u[1])==FALSE ){
      data.frame(function_tested="ensemblQueryR::ensemblQueryLDwithSNPpair",
                 n_queries=n,
                 time.sec=u$Elapsed_Time_sec,
                 peak.ram = u$Peak_RAM_Used_MiB,
                 total.ram = u$Total_RAM_Used_MiB) %>%
        vroom::vroom_write(x=., file=paste0("/home/abrowne/projects/ensemblQueryR/benchmarking/results/ram/", "ensemblQueryLDwithSNPpair", ".", n,".",i,".csv"))
    }else{
      print(u)
    }

    # -- 2b. test ensemblQueryR::ensemblQueryLDwithSNPpairDataframe --------------------------------------------

    # clear env
    rm(t, u)

    t = data.frame(
      rsid1 = c(rep("rs6792369", n)),
      rsid2 = c(rep("rs1042779", n))
    )

    u=try(peakRAM(
      ensemblQueryLDwithSNPpairDataframe(
        in.table=t,
        pop="1000GENOMES:phase_3:EUR",
        keep.original.table.row.n=FALSE)
    ))

    # save results
    if( grepl("Error", u[1])==FALSE ){
      data.frame(function_tested="ensemblQueryR::ensemblQueryLDwithSNPpairDataframe",
                 n_queries=n,
                 time.sec=u$Elapsed_Time_sec,
                 peak.ram = u$Peak_RAM_Used_MiB,
                 total.ram = u$Total_RAM_Used_MiB) %>%
        vroom::vroom_write(x=., file=paste0("/home/abrowne/projects/ensemblQueryR/benchmarking/results/ram/", "ensemblQueryLDwithSNPpairDataframe", ".", n,".",i,".csv"))
    }else{
      print(u)
    }

  }
}

# -- 2. Plot performance --------------------------------------------

library(magrittr)
library(dplyr)
library(ggplot2)
library(ggsci)
library(patchwork)

# Set defaults for ggplots
theme_rhr <- theme_set(
  theme_bw(base_family = "Helvetica",
           base_size = 10) +
    theme(panel.grid.major.x = element_blank(),
          legend.position = "top",
          legend.background = element_blank(),
          panel.spacing = unit(0.1, "lines"))
)

ram.line = list.files("/home/abrowne/projects/ensemblQueryR/benchmarking/results/ram",
                      pattern=".csv",
                      full.names=T) %>%
  lapply(X=., FUN=function(x){
    vroom::vroom(x)
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(package =
                  case_when(
                    grepl("ensemblQueryR", function_tested) ~ "ensemblQueryR",
                    grepl("LDlinkR", function_tested) ~ "LDlinkR"
                  ),
                function_tested=gsub("ensemblQueryR::", "", function_tested),
                function_tested=gsub("LDlinkR::", "", function_tested),
                n_queries = paste0(n_queries, " queries")
  ) %>%
  dplyr::mutate(n_queries = factor(n_queries, levels=c("100 queries", "1000 queries", "5000 queries", "10000 queries")),
                queries_numeric = readr::parse_number(as.character(n_queries))) %>%
  dplyr::filter(queries_numeric!=5000) %>%

  ggplot(data=., aes(x=queries_numeric, y=peak.ram, colour=package)) +
  geom_line(aes(linetype=function_tested)) +
  geom_point() +
  scale_x_continuous(breaks=c(100,1000,10000)) +
  # scale_y_continuous(n.breaks=20) +
  xlab("Number of queries") +
  ylab("Maximum RAM usage (MiB)") +
  scale_colour_npg() +
  theme(legend.title = element_blank()) +
  guides(linetype = guide_legend(nrow = 3),
         colour = guide_legend(nrow = 2))

# -- 3. Generate speed + RAM figure --------------------------------------------

(
  (
    (ram.line | speed.line) + plot_layout(guides = "collect") & theme(legend.position="bottom")
  ) & plot_annotation(tag_levels = "a")
) %>%
  ggsave(
    path="/home/abrowne/projects/ensemblQueryR/img/",
    filename="00-ram_speed_LDlinkR.png",
    plot = .,
    device = "png",
    scale = 1,
    width = 7,
    height = 4,
    units = c("in"),
    dpi = 300
  )

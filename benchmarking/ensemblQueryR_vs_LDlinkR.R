# load libs
devtools::load_all("/home/abrowne/projects/ensemblQueryR/")
library(LDlinkR)
library(magrittr)
library(peakRAM)

# -- 1. Run benchmarking for 100,1000,10000 queries --------------------------------------------

for(i in c(1:10)){
  for(n in c(100,1000,10000)){

    print(paste(
      "Processing", i, n
    ))

    # # -- 1. test LDlinkR --------------------------------------------
    #
    # # clear env
    # rm(t, u)
    #
    # t = list(
    #   c(rep("rs6792369", n)),
    #   c(rep("rs1042779", n))
    # )
    #
    # u=try(peakRAM(
    #   lapply(X=c(1:n),
    #          FUN=function(x){
    #            print(x)
    #            LDlinkR::LDpair(
    #              var1=t[[1]][x],
    #              var2=t[[2]][x],
    #              pop = "CEU",
    #              token = "5a141d21fa57",
    #              output = "table"
    #            ) %>%
    #              return()
    #
    #          }) %>%
    #     dplyr::bind_rows()
    # ))
    #
    # # save results
    # if( grepl("Error", u[1])==FALSE ){
    #   data.frame(function_tested="LDlinkR::LDpair",
    #              n_queries=n,
    #              time.sec=u$Elapsed_Time_sec,
    #              peak.ram = u$Peak_RAM_Used_MiB,
    #              total.ram = u$Total_RAM_Used_MiB,
    #              error = "none") %>%
    #     vroom::vroom_write(x=., file=paste0("/home/abrowne/projects/ensemblQueryR/benchmarking/results/ram/", "LDpair", ".", n,".",i,".csv"))
    # } else{
    #   data.frame(function_tested="LDlinkR::LDpair",
    #              n_queries=n,
    #              time.sec=NA_integer_,
    #              peak.ram =NA_integer_,
    #              total.ram =NA_integer_,
    #              error = u[1]) %>%
    #     vroom::vroom_write(x=., file=paste0("/home/abrowne/projects/ensemblQueryR/benchmarking/results/ram/", "LDpair", ".", n,".",i,".csv"))
    # }
    #
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
                 total.ram = u$Total_RAM_Used_MiB,
                 error = "none") %>%
        vroom::vroom_write(x=., file=paste0("/home/abrowne/projects/ensemblQueryR/benchmarking/results/ram/", "ensemblQueryLDwithSNPpair", ".", n,".",i,".csv"))
    }else{
      data.frame(function_tested="ensemblQueryR::ensemblQueryLDwithSNPpair",
                 n_queries=n,
                 time.sec=NA_integer_,
                 peak.ram =NA_integer_,
                 total.ram =NA_integer_,
                 error = u[1]) %>%
        vroom::vroom_write(x=., file=paste0("/home/abrowne/projects/ensemblQueryR/benchmarking/results/ram/", "ensemblQueryLDwithSNPpair", ".", n,".",i,".csv"))
    }

  }
}

# -- 2. Plot RAM --------------------------------------------

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
                      pattern="0.\\d+.csv",
                      full.names=T) %>%
  lapply(X=., FUN=function(x){
    vroom::vroom(x) %>%
      dplyr::mutate(iter = stringr::str_match_all(string=x, pattern="\\w+.\\d+.(\\d).csv")[[1]][2])
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(package =
                  case_when(
                    grepl("ensemblQueryR", function_tested) ~ "ensemblQueryR",
                    grepl("LDlinkR", function_tested) ~ "LDlinkR"
                  ),
                function_tested=gsub("ensemblQueryR::", "", function_tested),
                function_tested=gsub("LDlinkR::", "", function_tested),
                n_queries = paste0(n_queries, " queries")) %>%
  dplyr::mutate(n_queries = factor(n_queries, levels=c("100 queries", "1000 queries", "5000 queries", "10000 queries")),
                queries_numeric = readr::parse_number(as.character(n_queries))) %>%
  dplyr::mutate(time.mins = time.sec/60) %>%
  dplyr::group_by(package, function_tested, n_queries, queries_numeric) %>%
  dplyr::summarise(sd.peak.ram = sd(peak.ram, na.rm=T),
                   sd.time.mins = sd(time.mins, na.rm=T),
                   mean.peak.ram = mean(peak.ram, na.rm=T),
                   mean.time.mins = mean(time.mins, na.rm=T)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!grepl("Dataframe", function_tested)) %>%

  ggplot(., aes(x=queries_numeric, y=mean.peak.ram, group=package, color=package)) +
  geom_line() +
  geom_point(size=1) +
  geom_errorbar(aes(ymin=mean.peak.ram-sd.peak.ram, ymax=mean.peak.ram+sd.peak.ram), width=0.5) +
  scale_x_continuous(breaks=c(100,1000,10000)) +
  scale_y_continuous(n.breaks=20) +
  xlab("Number of queries") +
  ylab("Maximum RAM usage (MiB)") +
  scale_colour_npg() +
  theme(legend.title = element_blank()) +
  guides(linetype = guide_legend(nrow = 3),
         colour = guide_legend(nrow = 2))

# -- 3. Plot speed --------------------------------------------

speed.line = list.files("/home/abrowne/projects/ensemblQueryR/benchmarking/results/ram",
                        pattern="0.\\d+.csv",
                        full.names=T) %>%
  # .[!(grepl("LDpair.10000", .))] %>%
  lapply(X=., FUN=function(x){
    vroom::vroom(x) %>%
      dplyr::mutate(iter = stringr::str_match_all(string=x, pattern="\\w+.\\d+.(\\d).csv")[[1]][2])
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(package =
                  case_when(
                    grepl("ensemblQueryR", function_tested) ~ "ensemblQueryR",
                    grepl("LDlinkR", function_tested) ~ "LDlinkR"
                  ),
                function_tested=gsub("ensemblQueryR::", "", function_tested),
                function_tested=gsub("LDlinkR::", "", function_tested),
                n_queries = paste0(n_queries, " queries")) %>%
  dplyr::mutate(n_queries = factor(n_queries, levels=c("100 queries", "1000 queries", "5000 queries", "10000 queries")),
                queries_numeric = readr::parse_number(as.character(n_queries))) %>%
  dplyr::mutate(time.mins = time.sec/60) %>%
  dplyr::group_by(package, function_tested, n_queries, queries_numeric) %>%
  dplyr::summarise(sd.peak.ram = sd(peak.ram, na.rm=T),
                sd.time.mins = sd(time.mins, na.rm=T),
                mean.peak.ram = mean(peak.ram, na.rm=T),
                mean.time.mins = mean(time.mins, na.rm=T)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!grepl("Dataframe", function_tested)) %>%

  {
    ggplot(., aes(x=queries_numeric, y=mean.time.mins, group=package, color=package)) +
      geom_line() +
      geom_point(size=1) +
      geom_errorbar(aes(ymin=mean.time.mins-sd.time.mins, ymax=mean.time.mins+sd.time.mins), width=0.5) +
      scale_x_continuous(breaks=c(100,1000,10000)) +
      scale_y_continuous(breaks = c(seq(0,max(.$mean.time.mins)+30,15))) +
      xlab("Number of queries") +
      ylab("Time taken (mins)") +
      scale_colour_npg() +
      theme(legend.title = element_blank()) +
      guides(linetype = guide_legend(nrow = 3),
             colour = guide_legend(nrow = 2))
  }

# -- 4. Generate speed + RAM figure --------------------------------------------

(
  (
    (ram.line | speed.line) + plot_layout(guides = "collect") & theme(legend.position="right")
  ) & plot_annotation(tag_levels = "a")
) %>%
  ggsave(
    path="/home/abrowne/projects/ensemblQueryR/img/",
    filename="00-ram_speed_LDlinkR.png",
    plot = .,
    device = "png",
    scale = 1,
    width = 7,
    height = 3.5,
    units = c("in"),
    dpi = 300
  )

# -- 5. Check how many 10K LDpair runs failed --------------------------------------------

list.files("/home/abrowne/projects/ensemblQueryR/benchmarking/results/ram",
           pattern="0.\\d+.csv",
           full.names=T) %>%
  # .[!(grepl("LDpair.10000", .))] %>%
  lapply(X=., FUN=function(x){
    vroom::vroom(x, show_col_types = FALSE) %>%
      dplyr::mutate(iter = stringr::str_match_all(string=x, pattern="\\w+.\\d+.(\\d).csv")[[1]][2])
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::filter(n_queries==10000) %>%
  dplyr::filter(function_tested=="LDlinkR::LDpair") %>%
  dplyr::group_by(error) %>%
  dplyr::summarise(n=n())

# -- 6. Calculate performance stats --------------------------------------------

performance = list.files("/home/abrowne/projects/ensemblQueryR/benchmarking/results/ram",
                        pattern="0.\\d+.csv",
                        full.names=T) %>%
  # .[!(grepl("LDpair.10000", .))] %>%
  lapply(X=., FUN=function(x){
    vroom::vroom(x) %>%
      dplyr::mutate(iter = stringr::str_match_all(string=x, pattern="\\w+.\\d+.(\\d).csv")[[1]][2])
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(package =
                  case_when(
                    grepl("ensemblQueryR", function_tested) ~ "ensemblQueryR",
                    grepl("LDlinkR", function_tested) ~ "LDlinkR"
                  ),
                function_tested=gsub("ensemblQueryR::", "", function_tested),
                function_tested=gsub("LDlinkR::", "", function_tested),
                n_queries = paste0(n_queries, " queries")) %>%
  dplyr::mutate(n_queries = factor(n_queries, levels=c("100 queries", "1000 queries", "5000 queries", "10000 queries")),
                queries_numeric = readr::parse_number(as.character(n_queries))) %>%
  dplyr::mutate(time.mins = time.sec/60) %>%
  dplyr::filter(!grepl("Dataframe", function_tested)) %>%
  dplyr::group_by(package, function_tested, n_queries, queries_numeric) %>%
  dplyr::summarise(mean.peak.ram = mean(peak.ram, na.rm=T),
                   mean.time.mins = mean(time.mins, na.rm=T)) %>%
  dplyr::ungroup() %>%
  dplyr::select(package, n_queries, queries_numeric, mean.peak.ram, mean.time.mins) %>%
  tidyr::pivot_wider(names_from="package", values_from=c("mean.peak.ram", "mean.time.mins")) %>%
  # dplyr::mutate(RAM = ((mean.peak.ram_ensemblQueryR-mean.peak.ram_LDlinkR)/mean.peak.ram_LDlinkR)*100, # as a fraction of the RAM used by LDlinkR
  #               speed = ((mean.time.mins_ensemblQueryR-mean.time.mins_LDlinkR)/mean.time.mins_LDlinkR)*100) # how many times faster is ensemblQueryR
  dplyr::mutate(RAM = mean.peak.ram_ensemblQueryR/mean.peak.ram_LDlinkR,
                speed = mean.time.mins_LDlinkR/mean.time.mins_ensemblQueryR)

performance

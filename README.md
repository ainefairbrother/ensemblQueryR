
# ensemblQueryR

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/562138040.svg)](https://zenodo.org/badge/latestdoi/562138040)

<!-- badges: end -->

The goal of ensemblQueryR is to seemlessly integrate querying of Ensembl databases into your R workflow. It does this by formatting and submitting user queries to the Ensembl API. In its current iteration, the package can retreive SNPs in LD with any number of query SNPs. 

## Installation

You can install the development version of ensemblQueryR like so:

``` r
remotes::install_github("ainefairbrother/ensemblQueryR")

```

## Example: for one query SNP  

``` r
# load libraries
library(ensemblQueryR)
library(magrittr)
```

``` r
# get all SNPs in LD with query SNP
ensemblQueryR::ensemblQueryLDwithSNP(rsid="rs3851179", 
                      r2=0.8, 
                      d.prime=0.8, 
                      window.size=500, 
                      pop="1000GENOMES:phase_3:EUR")
```

To get a list of possible human Ensembl populations to use in the `pop` argument, run the `ensemblQueryGetPops()` function.

``` r
ensemblQueryR::ensemblQueryGetPops()
```

## Example: for <1000 query SNPs

``` r
rsid.list <- c("rs7153434","rs1963154","rs12672022","rs3852802","rs12324408","rs56346870")

# run query on rsid.list
ensemblQueryR::ensemblQueryLDwithSNPlist(rsid.list, 
                          r2=0.8, 
                          d.prime=0.8, 
                          window.size=500, 
                          pop="1000GENOMES:phase_3:EUR")
``` 

## Example for >1000 query SNPs

There is a separate function for large queries (>1000 SNPs) because of Ensembl's API query size limit. This function takes a `data.frame` as an input, and gets all SNPs in LD with a column containing query SNPs called `rsid`. 

``` r
# example input data
in.table <- data.frame(rsid=rep(c("rs7153434","rs1963154","rs12672022","rs3852802","rs12324408","rs56346870"), 500))

# run query on in.table
ensemblQueryR::ensemblQueryLDwithLargeSNPdf(in.table=in.table,
                             r2=0.8,
                             d.prime=0.8,
                             window.size=500,
                             pop="1000GENOMES:phase_3:EUR")
```

## Disclaimer

Please note that this code is still under development and may contain bugs or errors. It is not recommended for use in production environments. Use at your own risk. I am working on improving the code, addressing any issues, and expanding the package's capabilities so please check back for updates.


<!-- badges: start -->

[![DOI](https://zenodo.org/badge/562138040.svg)](https://zenodo.org/badge/latestdoi/562138040)

<!-- badges: end -->

# ensemblQueryR

The goal of ensemblQueryR is to seemlessly integrate querying of Ensembl databases into your R workflow. It does this by formatting and submitting user queries to the Ensembl API. At present, the package contains functions for the three Ensembl Linkage Disequilibrium (LD) 'endpoints': 1. Query LD in a window around one SNP, 2. Query LD for a pair of query SNPs and 3. Query LD for SNPs at a specified genomic locus. 

## Installation

You can install the development version of ensemblQueryR like so:

``` r
# load remotes package
library(remotes)

remotes::install_github("ainefairbrother/ensemblQueryR")
```

## Setup 

All functions in this package take the `pop` argument which defines the population for which to retrieve LD statistics. So, to get a list of options for this argument, run the `ensemblQueryGetPops()` function.

``` r
# load libraries
library(ensemblQueryR)
library(magrittr)

ensemblQueryR::ensemblQueryGetPops()
```

## Functionality 1: querying LD for window around one query rsID

#### For 1 query rsID

For one query rsID, get all rsIDs in LD using `ensemblQueryLDwithSNPwindow`

``` r
ensemblQueryR::ensemblQueryLDwithSNPwindow(rsid="rs3851179", 
                      r2=0.8, 
                      d.prime=0.8, 
                      window.size=500, 
                      pop="1000GENOMES:phase_3:EUR")
```

#### For >1 and <1000 query rsIDs

For a vector of query rsIDs, get all rsIDs in LD if your query is <1000 rsIDs in length. This is due to Ensembl's 1000 query limit. See next example for queries >1000 rsIDs in length.

``` r
rsid.vec <- c("rs7153434","rs1963154","rs12672022","rs3852802","rs12324408","rs56346870")

# run query on rsid.vec
ensemblQueryR::ensemblQueryLDwithSNPwindowList(rsid.vec, 
                          r2=0.8, 
                          d.prime=0.8, 
                          window.size=500, 
                          pop="1000GENOMES:phase_3:EUR")
``` 

#### For >1000 query rsIDs

There is a separate function for large queries (>1000 rsIDs) because of Ensembl's API query size limit. This function takes a `data.frame` as an input, and gets all rsIDs in LD with a column of query rsIDs called `rsid`. 

``` r
# example input data
in.table <- data.frame(rsid=rep(c("rs7153434","rs1963154","rs12672022","rs3852802","rs12324408","rs56346870"), 500))

# run query on in.table
ensemblQueryR::ensemblQueryLDwithSNPwindowDataframe(
  in.table=in.table,
  r2=0.8,
  d.prime=0.8,
  window.size=500,
  pop="1000GENOMES:phase_3:EUR"
)
```

## Functionality 2: querying LD for a pair of query SNPs

The `ensemblQueryLDwithSNPpair` takes a single pair of query SNPs and provides a `data.frame` of LD statistics for these. 

``` r
ensemblQueryLDwithSNPpair(
  rsid1="rs6792369",
  rsid2="rs1042779",
  pop="1000GENOMES:phase_3:EUR"
)
```

The `ensemblQueryLDwithSNPpairDataframe` takes a `data.frame` with columns named `rsid1` and `rsid2` containing many pairs of query SNPs and provides a `data.frame` of LD statistics for these. 

``` r
ensemblQueryLDwithSNPpairDataframe(
  in.table=data.frame(rsid1=rep("rs6792369", 10), rsid2=rep("rs1042779", 10)),
  pop="1000GENOMES:phase_3:EUR",
  keep.original.table.row.n=F,
  parallelise=F
)
```

## Functionality 3: querying LD for a genomic region

The `ensemblQueryLDwithSNPregion` function takes a genomic region as input and returns all SNP pairs and their LD statistics within the defined region.

``` r
ensemblQueryLDwithSNPregion(
  chr="6",
  start="25837556",
  end="25843455",
  pop="1000GENOMES:phase_3:EUR"
)
```

## Disclaimer

Please note that this code is still under development and may contain bugs or errors. It is not recommended for use in production environments. Use at your own risk. I am working on improving the code, addressing any issues, and expanding the package's capabilities so please check back for updates.

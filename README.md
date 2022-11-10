
# ensemblQueryR

<!-- badges: start -->

<!-- badges: end -->

The goal of ensemblQueryR is to seemlessly integrate querying of Ensembl databases into your R workflow. It does this by formatting and submitting user queries to the Ensembl API. In its current iteration, the package can retreive SNPs in LD with any number of query SNPs. 

## Installation

You can install the development version of ensemblQueryR like so:

``` r
remotes::install_github("ainefairbrother/ensemblQueryR")
```

## Example

``` r
library(ensemblQueryR)
library(magrittr)

# for one query SNP: get all SNPs in LD with query SNP
ensemblQueryLDwithSNP(rsid="rs3851179", 
                      r2=0.8, 
                      d.prime=0.8, 
                      window.size=500, 
                      pop="1000GENOMES:phase_3:EUR")

# for fewer than 1000 query SNPs: get all SNPs in LD with a list of query SNPs
rsid.list = c("rs7153434","rs1963154","rs12672022","rs3852802","rs12324408","rs56346870")

# run query on rsid.list
ensemblQueryLDwithSNPlist(rsid.list, 
                          r2=0.8, 
                          d.prime=0.8, 
                          window.size=500, 
                          pop="1000GENOMES:phase_3:EUR")

# for more than 1000 query SNPs: get all SNPs in LD with a data.frame column of query SNPs
# example input data
in.table <- data.frame(rsid=rep(c("rs7153434","rs1963154","rs12672022","rs3852802","rs12324408","rs56346870"), 500))

# run query on in.table
ensemblQueryLDwithLargeSNPdf(in.table=in.table,
                             r2=0.8,
                             d.prime=0.8,
                             window.size=500,
                             pop="1000GENOMES:phase_3:EUR")

```


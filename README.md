
<!-- badges: start -->

[![DOI](https://zenodo.org/badge/562138040.svg)](https://zenodo.org/badge/latestdoi/562138040)
[![HitCount](https://hits.dwyl.com/dwyl/start-here.svg)](https://hits.dwyl.com/dwyl/start-here)

<!-- badges: end -->

# ensemblQueryR

The goal of ensemblQueryR is to seemlessly integrate querying of Ensembl databases into your R workflow. It does this by formatting and submitting user queries to the Ensembl API. At present, the package contains functions for the three Ensembl Linkage Disequilibrium (LD) 'endpoints': 1. Query LD in a window around one SNP, 2. Query LD for a pair of query SNPs and 3. Query LD for SNPs at a specified genomic locus. 

## Installation

You can install ensemblQueryR as below.

``` r
# load remotes package
library(remotes)

# to install the development version
remotes::install_github("ainefairbrother/ensemblQueryR")

# to install the stable CRAN release
install.packages("ensemblQueryR")
```

## Setup 

To check that the Ensembl server is up and running, the server can be pinged. 

``` r
library(ensemblQueryR)

ensemblQueryR::pingEnsembl()
```

All functions in this package take the `pop` argument which defines the population for which to retrieve LD metrics. To get a list of options for this argument, run the `ensemblQueryGetPops()` function.

``` r
ensemblQueryR::ensemblQueryGetPops()
```

## 1. Query LD metrics for a window around a variant  

Get all variants in LD with one query variant using `ensemblQueryLDwithSNPwindow`. This function constrains the query by taking a minimum r-squared cut-off (`r2`), D-prime (`d.prime`) and window size around the variant in kilobases (`window.size`).

``` r
ensemblQueryR::ensemblQueryLDwithSNPwindow(rsid="rs3851179", 
  r2=0.8, 
  d.prime=0.8, 
  window.size=500, 
  pop="1000GENOMES:phase_3:EUR")
```

For more than one query variant, the `ensemblQueryLDwithSNPwindowDataframe` function takes a `data.frame` as input, and gets all variants in LD with all query variants in the `rsid` column. It is possible to parallelise this operation by setting the number of cores above 1.

``` r
# example input data
in.table <- data.frame(rsid=rep(c("rs7153434","rs1963154","rs12672022","rs3852802","rs12324408","rs56346870"), 500))

# run query on in.table
ensemblQueryR::ensemblQueryLDwithSNPwindowDataframe(
  in.table=in.table,
  r2=0.8,
  d.prime=0.8,
  window.size=500,
  pop="1000GENOMES:phase_3:EUR",
  cores=1
)
```

## 2. Query LD metrics for a pair of variants  

The `ensemblQueryLDwithSNPpair` takes a single pair of query SNPs and returns a `data.frame` of LD metrics.

``` r
ensemblQueryR::ensemblQueryLDwithSNPpair(
  rsid1="rs6792369",
  rsid2="rs1042779",
  pop="1000GENOMES:phase_3:EUR"
)
```

The `ensemblQueryLDwithSNPpairDataframe` takes a `data.frame` with columns `rsid1` and `rsid2` and returns a `data.frame` of LD metrics for all variant pairs. It is possible to parallelise this operation by setting the number of cores above 1.

``` r
# example input data
in.table <- data.frame(rsid1=rep("rs6792369", 10), rsid2=rep("rs1042779", 10))

# run query on in.table
ensemblQueryR::ensemblQueryLDwithSNPpairDataframe(
  in.table=in.table,
  pop="1000GENOMES:phase_3:EUR",
  keep.original.table.row.n=F,
  cores=1
)
```

## 3. Query LD metrics for a genomic region  

The `ensemblQueryLDwithSNPregion` function takes genomic coordinates as input and returns all variant pairs and their LD metrics within the defined region.

``` r
ensemblQueryR::ensemblQueryLDwithSNPregion(
  chr="6",
  start="25837556",
  end="25843455",
  pop="1000GENOMES:phase_3:EUR"
)
```

The `ensemblQueryLDwithSNPregionDataframe` takes a `data.frame` with columns `chr`, `start` and `end` and returns a `data.frame` of LD metrics for all variant pairs contained within each genomic region (each row of `in.table`). It is possible to parallelise this operation by setting the number of cores above 1.

```r
# example input data
in.table = data.frame(chr=rep(c("6"), 10),
                       start=rep(c("25837556"), 10),
                       end=rep(c("25843455"), 10))
                       
# run query on in.table
ensemblQueryR::ensemblQueryLDwithSNPregionDataframe(
  in.table= ,
  pop="1000GENOMES:phase_3:EUR",
  cores = 2
)
```

## Docker

We have provided a [Docker](https://www.docker.com/) image, enabling this tool to be run regardless of your local operating system or R version. This can be found [here](https://hub.docker.com/r/ainefairbrotherbrowne/ensemblqueryr/tags). As long as you have Docker installed, the code below will allow you to pull this image, run a container and execute it. You will then be able to use `ensemblQueryR` as described above. A working installation of Docker is required. 

```bash
docker pull ainefairbrotherbrowne/ensemblqueryr:1.0; \
docker run -t -d --name ensemblqueryr ainefairbrotherbrowne/ensemblqueryr:1.0; \ 
docker exec -i -t ensemblqueryr R
```

Aditionally, to mount a volume - enabling you to load a file containing your variant IDs, for example - the following command can be used, replacing `path/to/vol` with the path to the directory you wish to mount.

```bash
docker pull ainefairbrotherbrowne/ensemblqueryr:1.0; \
docker run -t -d --name ensemblqueryr ainefairbrotherbrowne/ensemblqueryr:1.0 --volume path/to/vol; \ 
docker exec -i -t ensemblqueryr R
```

## Singularity 

For HPC use-cases where Docker usage becomes problematic owing to user privilege limitations, we have provided a [singularity](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html) image. This can be found [here](https://cloud.sylabs.io/library/ainefairbrother/ensemblqueryr/ensemblqueryr). The code below will allow you to pull this image, run a container and execute it. You will then be able to use `ensemblQueryR` as described above. A working installation of singularity is required. 

```bash
singularity pull --arch amd64 library://ainefairbrother/ensemblqueryr/ensemblqueryr:sha256.e387ea11ae4eaea8f94d81c625c2c1d5a22dd351858ebcd03910a7736d76ca30; \
singularity exec ensemblqueryr_sha256.e387ea11ae4eaea8f94d81c625c2c1d5a22dd351858ebcd03910a7736d76ca30.sif R
```

## Contribute to `ensemblQueryR`

We value contributions from the community to improve `ensemblQueryR`. Here's how you can do this:\
* Explore Issues: Look through our GitHub Issues to find tasks, bugs, or new features you're interested in.\
* Fork the Repository: Fork the `ensemblQueryR` repository to your GitHub account using the "Fork" button at the top.\
* Branch Out: Create a new branch in your forked repository for your changes. Use a descriptive name for clarity.\
* Make Changes: Make your code, documentation, or design updates in the branch you've created.\
* Commit and Push: Commit your changes and push them to your forked repository.\
* Open a Pull Request (PR): Submit a PR from your branch to the main repository. Provide a clear title and description of your changes.\
* Engage in Discussion: Participate in the conversation around your PR. Address feedback and make necessary adjustments.\
* Collaborate: Collaborate with our maintainers to refine your contribution. Once approved, your changes will be merged.\

Thank you for considering making a contribution to `ensemblQueryR`.

## Disclaimer

Please note that this code is still under development and may contain bugs or errors. It is not recommended for use in production environments. Use at your own risk. I am working on improving the code, addressing any issues, and expanding the package's capabilities so please check back for updates.

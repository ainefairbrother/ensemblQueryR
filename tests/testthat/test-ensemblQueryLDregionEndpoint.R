# prevents the testing hanging as per https://github.com/r-lib/devtools/issues/1526
Sys.setenv(R_TESTS="")

test_that("ensemblQueryLDwithSNPregion returns data.frame", {
  expect_s3_class(

    ensemblQueryLDwithSNPregion(
      chr="6",
      start="25837556",
      end="25843455",
      pop="1000GENOMES:phase_3:EUR"
    ),

    "data.frame")
})

test_that("ensemblQueryLDwithSNPregion returns data.frame", {
  expect_s3_class(

    ensemblQueryLDwithSNPregion(
      chr="!",
      start="hello",
      end="world",
      pop="1000GENOMES:phase_3:EUR"
    ),

    "data.frame")
})

test_that("ensemblQueryLDwithSNPregionDataframe returns data.frame", {
  expect_s3_class(

    ensemblQueryLDwithSNPregionDataframe(
      in.table=data.frame(
        chr=rep(c("6"), 10),
        start=rep(c("badrsnumber"), 10),
        end=rep(c("badrsnumber"), 10)
      ),
      pop="1000GENOMES:phase_3:EUR",
      cores = 2
    ),

    "data.frame")
})

test_that("ensemblQueryLDwithSNPregionDataframe returns data.frame", {
  expect_s3_class(

    ensemblQueryLDwithSNPregionDataframe(
      in.table=data.frame(
        chr=rep(c("6"), 10),
        start=rep(c("25837556"), 10),
        end=rep(c("25843455"), 10)
      ),
      pop="1000GENOMES:phase_3:EUR",
      cores = 2
    ),

    "data.frame")
})



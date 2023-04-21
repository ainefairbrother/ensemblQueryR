library(testthat)

test_that("ensemblQueryLDwithSNPpair returns data.frame", {
  expect_s3_class(

    ensemblQueryLDwithSNPpair(
      rsid1="rs6792369",
      rsid2="rs1042779",
      pop="1000GENOMES:phase_3:EUR"),

    "data.frame")
})

test_that("ensemblQueryLDwithSNPpair returns data.frame", {
  expect_s3_class(

    ensemblQueryLDwithSNPpair(
      rsid1="badrsnumber",
      rsid2="otherbadrsnumber",
      pop="1000GENOMES:phase_3:EUR"),

    "data.frame")
})

test_that("ensemblQueryLDwithSNPpairDataframe returns data.frame", {
  expect_s3_class(

    ensemblQueryLDwithSNPpairDataframe(
      in.table=data.frame(rsid1=rep("rs6792369", 10), rsid2=rep("rs1042779", 10)),
      pop="1000GENOMES:phase_3:EUR",
      keep.original.table.row.n=FALSE,
      parallelise=FALSE),

    "data.frame")
})

test_that("ensemblQueryLDwithSNPpairDataframe returns data.frame", {
  expect_s3_class(

    ensemblQueryLDwithSNPpairDataframe(
      in.table=data.frame(rsid1=rep("badrsnumber", 10), rsid2=rep("badrsnumber", 10)),
      pop="1000GENOMES:phase_3:EUR",
      keep.original.table.row.n=FALSE,
      parallelise=FALSE),

    "data.frame")
})

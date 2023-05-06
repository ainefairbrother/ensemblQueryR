# prevents the testing hanging as per https://github.com/r-lib/devtools/issues/1526
Sys.setenv(R_TESTS="")

test_that("ensemblQueryGetPops returns data.frame", {
  expect_s3_class(

    ensemblQueryGetPops(),

    "data.frame")
})

test_that("ensemblQueryLDwithSNPwindow returns data.frame", {

  expect_s3_class(

    ensemblQueryLDwithSNPwindow(
      rsid="rs3851179",
      r2=0.8,
      d.prime=0.8,
      window.size=500,
      pop="1000GENOMES:phase_3:EUR"),

    "data.frame")
})

test_that("ensemblQueryLDwithSNPwindow returns data.frame", {

  expect_s3_class(

    ensemblQueryLDwithSNPwindow(
      rsid="junk",
      r2=0.8,
      d.prime=0.8,
      window.size=500,
      pop="1000GENOMES:phase_3:EUR"),

    "data.frame")
})

test_that("ensemblQueryLDwithSNPwindowDataframe returns data.frame", {

  expect_s3_class(

    ensemblQueryLDwithSNPwindowDataframe(
      in.table=data.frame(rsid=rep(c("crap","crap","crap","crap","crap","crap"), 10)),
      r2=0.8,
      d.prime=0.8,
      window.size=500,
      pop="1000GENOMES:phase_3:EUR"),

    "data.frame")
})

test_that("ensemblQueryLDwithSNPwindowDataframe returns data.frame", {

  expect_s3_class(

    ensemblQueryLDwithSNPwindowDataframe(
        in.table=data.frame(rsid=rep(c("rs7153434","rs1963154","rs12672022","rs3852802","rs12324408","rs56346870"), 10)),
        r2=0.8,
        d.prime=0.8,
        window.size=500,
        pop="1000GENOMES:phase_3:EUR"),

    "data.frame")
})

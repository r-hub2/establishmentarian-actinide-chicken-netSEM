test_that("netSEMp2 returns correct class", {
  result <- netSEMp2(acrylic, exogenous = "IrradTot", endogenous = "YI")
  expect_s3_class(result, c("netSEMp2", "list"))
})

test_that("netSEMp2 returns expected list elements", {
  result <- netSEMp2(acrylic, exogenous = "IrradTot", endogenous = "YI")
  
  expect_named(result, c("data", "res.print", "res.full.model", "res.model"))
  expect_s3_class(result$data, "data.frame")
  expect_s3_class(result$res.print, "data.frame")
})

test_that("netSEMp2 res.print contains required columns", {
  result <- netSEMp2(acrylic, exogenous = "IrradTot", endogenous = "YI")
  
  required_cols <- c("endogenous", "Variable", "GModel", "GAdj-R2")
  expect_true(all(required_cols %in% names(result$res.print)))
})

test_that("netSEMp2 errors when exogenous/endogenous missing", {
  expect_error(netSEMp2(acrylic), "exogenous and endogenous variables need to be specified")
  expect_error(netSEMp2(acrylic, exogenous = "IrradTot"), 
               "exogenous and endogenous variables need to be specified")
})

test_that("netSEMp2 works with BIC criterion", {
  result <- netSEMp2(acrylic, exogenous = "IrradTot", endogenous = "YI", criterion = "BIC")
  expect_s3_class(result, c("netSEMp2", "list"))
})

test_that("netSEMp2 works with PVmodule dataset", {
  result <- netSEMp2(PVmodule, exogenous = "Time", endogenous = "Pmax")
  expect_s3_class(result, c("netSEMp2", "list"))
})
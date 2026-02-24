test_that("netSEMp1 returns correct class", {
  result <- netSEMp1(acrylic, exogenous = "IrradTot", endogenous = "YI")
  expect_s3_class(result, c("netSEMp1", "list"))
})

test_that("netSEMp1 returns expected list elements", {
  result <- netSEMp1(acrylic, exogenous = "IrradTot", endogenous = "YI")
  
  expect_named(result, c("data", "exogenous", "endogenous", "bestModels", "allModels"))
  expect_equal(result$exogenous, "IrradTot")
  expect_equal(result$endogenous, "YI")
  expect_s3_class(result$bestModels, "data.frame")
  expect_type(result$allModels, "list")
})

test_that("netSEMp1 bestModels contains required columns", {
  result <- netSEMp1(acrylic, exogenous = "IrradTot", endogenous = "YI")
  
  required_cols <- c("Resp", "Var", "best_model", "adj.r.squared")
  expect_true(all(required_cols %in% names(result$bestModels)))
})

test_that("netSEMp1 errors when exogenous/endogenous missing", {
  expect_error(netSEMp1(acrylic), "exogenous and endogenous variables need to be specified")
  expect_error(netSEMp1(acrylic, exogenous = "IrradTot"), 
               "exogenous and endogenous variables need to be specified")
  expect_error(netSEMp1(acrylic, endogenous = "YI"), 
               "exogenous and endogenous variables need to be specified")
})

test_that("netSEMp1 works with PVmodule dataset", {
  result <- netSEMp1(PVmodule, exogenous = "Time", endogenous = "Pmax")
  
  expect_s3_class(result, c("netSEMp1", "list"))
  expect_equal(result$exogenous, "Time")
  expect_equal(result$endogenous, "Pmax")
})

test_that("netSEMp1 preserves data dimensions", {
  result <- netSEMp1(acrylic, exogenous = "IrradTot", endogenous = "YI")
  
  expect_equal(nrow(result$data), nrow(acrylic))
  expect_equal(ncol(result$data), ncol(acrylic))
})
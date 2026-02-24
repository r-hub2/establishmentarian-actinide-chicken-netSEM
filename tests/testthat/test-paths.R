# Test paths function (if exported)

test_that("paths function works with netSEMp1 output", {
  skip_if_not(exists("paths", envir = asNamespace("netSEM")), 
              "paths function not exported")
  
  x <- netSEMp1(acrylic, exogenous = "IrradTot", endogenous = "YI")
  
  # Test that paths function runs without error
  result <- paths(x)
  expect_type(result, "list")
})
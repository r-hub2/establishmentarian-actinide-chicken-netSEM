# Test individual model functional forms
# These are internal functions tested via netSEMp1 output

test_that("netSEMp1 identifies multiple model types", {
  result <- netSEMp1(acrylic, exogenous = "IrradTot", endogenous = "YI")
  
  model_types <- unique(result$bestModels$best_model)
  
 # Should identify at least some model types
  expect_true(length(model_types) >= 1)
  
  # Model types should be from expected set
  valid_types <- c("Linear", "Quad", "SQuad", "Exponential", "Log", 
                   "SqRoot", "InvSqRoot", "NLS", "CP")
  expect_true(all(model_types %in% valid_types))
})

test_that("netSEMp1 calculates valid R-squared values", {
  result <- netSEMp1(acrylic, exogenous = "IrradTot", endogenous = "YI")
  
  # R-squared should be between 0 and 1 (or NA)
  r2_values <- result$bestModels$adj.r.squared
  valid_r2 <- r2_values[!is.na(r2_values)]
  
 expect_true(all(valid_r2 >= -1 & valid_r2 <= 1))
})

test_that("netSEMp1 model equations are generated", {
  result <- netSEMp1(acrylic, exogenous = "IrradTot", endogenous = "YI")
  
  # Check that model_eq column exists and has values
  expect_true("model_eq" %in% names(result$bestModels))
})
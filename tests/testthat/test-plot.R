# Setup - create model objects once for reuse

# plot.netSEMp1 tests
test_that("plot.netSEMp1 returns correct class", {
  x <- netSEMp1(acrylic, exogenous = "IrradTot", endogenous = "YI")
  result <- plot(x)
  expect_s3_class(result, c("grViz", "htmlwidget"))
})

test_that("plot.netSEMp1 works with cutoff parameter", {
  x <- netSEMp1(acrylic, exogenous = "IrradTot", endogenous = "YI")
  
  # Single cutoff

  result1 <- plot(x, cutoff = c(0.5))
  expect_s3_class(result1, c("grViz", "htmlwidget"))
  
  # Two cutoffs
  result2 <- plot(x, cutoff = c(0.3, 0.6))
  expect_s3_class(result2, c("grViz", "htmlwidget"))
  
  # Three cutoffs
  result3 <- plot(x, cutoff = c(0.3, 0.6, 0.8))
  expect_s3_class(result3, c("grViz", "htmlwidget"))
})

test_that("plot.netSEMp1 works with title parameter", {
  x <- netSEMp1(acrylic, exogenous = "IrradTot", endogenous = "YI")
  result <- plot(x, cutoff = c(0.3, 0.6, 0.8), title = "Test Title")
  expect_s3_class(result, c("grViz", "htmlwidget"))
})

test_that("plot.netSEMp1 works with style parameter", {
  x <- netSEMp1(acrylic, exogenous = "IrradTot", endogenous = "YI")
  
  result_style_true <- plot(x, cutoff = c(0.3, 0.6, 0.8), style = TRUE)
  result_style_false <- plot(x, cutoff = c(0.3, 0.6, 0.8), style = FALSE)
  
  expect_s3_class(result_style_true, c("grViz", "htmlwidget"))
  expect_s3_class(result_style_false, c("grViz", "htmlwidget"))
})

test_that("plot.netSEMp1 works with latent parameter", {
  x <- netSEMp1(acrylic, exogenous = "IrradTot", endogenous = "YI")
  result <- plot(x, cutoff = c(0.3, 0.6, 0.8), 
                 latent = c('IAD1' = 'FundAbsEdge'))
  expect_s3_class(result, c("grViz", "htmlwidget"))
})

# plot.netSEMp2 tests
test_that("plot.netSEMp2 returns correct class", {
  x <- netSEMp2(acrylic, exogenous = "IrradTot", endogenous = "YI")
  result <- plot(x)
  expect_s3_class(result, c("grViz", "htmlwidget"))
})

test_that("plot.netSEMp2 works with cutoff parameter", {
  x <- netSEMp2(acrylic, exogenous = "IrradTot", endogenous = "YI")
  
  result <- plot(x, cutoff = c(0.3, 0.6, 0.8))
  expect_s3_class(result, c("grViz", "htmlwidget"))
})

test_that("plot.netSEMp2 works with title parameter", {
  x <- netSEMp2(acrylic, exogenous = "IrradTot", endogenous = "YI")
  result <- plot(x, cutoff = c(0.3, 0.6, 0.8), title = "Test Title")
  expect_s3_class(result, c("grViz", "htmlwidget"))
})

test_that("plot.netSEMp2 works with style parameter", {
  x <- netSEMp2(acrylic, exogenous = "IrradTot", endogenous = "YI")
  
  result_style_false <- plot(x, cutoff = c(0.3, 0.6, 0.8), style = FALSE)
  expect_s3_class(result_style_false, c("grViz", "htmlwidget"))
})
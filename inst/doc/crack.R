## ----message=FALSE, eval=FALSE------------------------------------------------
# 
# ## Load the crack data set and preview column metadata
# data(crack)
# ?crack
# 
# ## Run netSEMp1 model
# ans1 <- netSEMp1(crack, exogenous = "uva360dose", endogenous = "dAvgNorm")
# ## Plot the network model for principle 1
# plot(ans1, cutoff = c(0.4, 0.5, 0.6))
# 
# ## Run netSEMp2 model
# ans2 <- netSEMp2(crack, exogenous = "uva360dose", endogenous = "dAvgNorm")
# ## Plot the network model for principle 2
# plot(ans2, cutoff = c(0.4, 0.5, 0.6))

## ----out.width="675px", dpi=1000, echo=FALSE, fig.cap="Crack netSEMp1 model"----
knitr::include_graphics("crack1.png")

## ----out.width="800px", dpi=1000, echo=FALSE, fig.cap="Crack netSEMp2 model"----
knitr::include_graphics("crack2.png")


## ----message=FALSE, eval=FALSE------------------------------------------------
# 
# ## Load the acrylic data set
# data("IVfeature")
# 
# ## Run netSEMp1 model
# ans1 <- netSEMp1(IVfeature, exogenous = "dy", endogenous = "Pmp")
# ## Plot the network model for principle 1
# plot(ans1, cutoff = c(0.25, 0.5, 0.8))
# 
# ## Run netSEMp2 model
# ans2 <- netSEMp2(IVfeature, exogenous = "dy", endogenous = "Pmp")
# ## Plot the network model for principle 2
# plot(ans2, cutoff = c(0.25, 0.5, 0.8))

## ----out.width="675px", echo=FALSE, fig.cap="IVfeature netSEMp1 model"--------
knitr::include_graphics("IVfeature1.png")

## ----out.width="800px", echo=FALSE, fig.cap="IVfeature netSEMp2 model"--------
knitr::include_graphics("IVfeature2.png")


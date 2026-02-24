## ----message=FALSE, eval=FALSE------------------------------------------------
# 
# ## Load the acrylic data set
# data(acrylic)
# ?acrylic
# 
# ## Run netSEMp1 model
# ans1 <- netSEMp1(acrylic, exogenous = "IrradTot", endogenous = "YI")
# ## Plot the network model for principle 1
# plot(ans1, cutoff = c(0.3,0.6,0.9))
# 
# ## Run netSEMp2 model
# ans2 <- netSEMp2(acrylic, exogenous = "IrradTot", endogenous = "YI")
# ## Plot the network model for principle 2
# plot(ans2, cutoff = c(0.3,0.6,0.9))

## ----out.width="675px", echo=FALSE, fig.cap="Acrylic netSEMp1 model"----------
knitr::include_graphics("acrylic1.png")

## ----out.width="800px", echo=FALSE, fig.cap="Acrylic netSEMp2 model"----------
knitr::include_graphics("acrylic2.png")


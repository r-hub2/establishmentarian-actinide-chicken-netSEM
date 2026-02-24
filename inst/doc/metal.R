## ----message=FALSE, eval=FALSE------------------------------------------------
# 
# ## Load the metal data set
# data(metal)
# ?metal
# 
# ## Run netSEMp1 model
# ans1 <- netSEMp1(metal, exogenous = "Z", endogenous = "Hardness")
# ## Plot the network model for principle 1
# plot(ans1, cutoff = c(0.3,0.6,0.9))
# 
# ## Run netSEMp2 model
# ans2 <- netSEMp2(metal, exogenous = "Z", endogenous = "Hardness", criterion = "AIC")
# ## Plot the network model for principle 2
# plot(ans2, cutoff = c(0.3,0.6,0.9))

## ----out.width="675px", echo=FALSE, fig.cap="Metal netSEMp1 model"------------
knitr::include_graphics("metal1.png")

## ----out.width="675px", echo=FALSE, fig.cap="Metal netSEMp2 model"------------
knitr::include_graphics("metal2.png")


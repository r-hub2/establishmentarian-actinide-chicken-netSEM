## ----message=FALSE, eval=FALSE------------------------------------------------
# 
# ## Load the backsheet data set
# data(backsheet)
# 
# ## Run netSEMp1 model
# ans1 <- netSEMp1(backsheet, exogenous = "Hours", endogenous = "YI")
# ## Plot the network model for principle 1
# plot(ans1, cutoff = c(0.3, 0.6, 0.9))
# 
# ## Run netSEMp2 model
# ans2 <- netSEMp2(backsheet, exogenous = "Hours", endogenous = "YI")
# ## Plot the network model for principle 2
# plot(ans2, cutoff = c(0.3, 0.6, 0.9))
# 
# 

## ----out.width="675px", echo=FALSE, fig.cap="Backsheet netSEMp1 model"--------
knitr::include_graphics("backsheet1.png")

## ----out.width="800px", echo=FALSE, fig.cap="Backsheet netSEMp2 model"--------
knitr::include_graphics("backsheet2.png")


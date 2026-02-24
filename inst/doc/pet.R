## ----message=FALSE, eval=FALSE------------------------------------------------
# 
# ## Load the PET data set
# data(pet)
# ?pet
# 
# ## Run netSEMp1 model
# ans1 <- netSEMp1(pet, exogenous = "Time", endogenous = "YI")
# ## Plot the network model for principle 1
# plot(ans1, cutoff = c(0.3,0.6,0.9))
# 
# ## Run netSEMp2 model
# ans2 <- netSEMp2(pet, exogenous = "Time", endogenous = "YI")
# ## Plot the network model for principle 2
# plot(ans2, cutoff = c(0.3,0.6,0.9))

## ----out.width="675px", echo=FALSE, fig.cap="PET netSEMp1 model"--------------
knitr::include_graphics("pet1.png")

## ----out.width="675px", echo=FALSE, fig.cap="PET netSEMp2 model"--------------
knitr::include_graphics("pet2.png")

## ----message=FALSE, eval=FALSE------------------------------------------------
# # Use pathwayRMSE to examine errors of network paths
# pathwayRMSE(ans)


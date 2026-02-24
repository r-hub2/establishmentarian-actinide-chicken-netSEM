
# netSEM

<!-- badges: start -->
<!-- badges: end -->

The R package 'netSEM' conducts a net-SEM statistical analysis (network structural equation modeling) on a data frame of coincident observations of multiple continuous variables. Principle 1 generates an inferential model through pairwise correlation of variables based on the Markovian Spirit. Principle 2 provides a predictive model through multiple regression with the model complexity and performance evaluated using either Akaike Information Criterion (AIC) or Bayesian Information Criterion (BIC) specified by the user.


## netSEM Usage

This is a simple example for generating degradation models using netSEM principle 1 and principle 2.

``` r
# Load in netSEM library
library(netSEM)
# Load in example data
data(acrylic)
# Perform principle 1 and principle 2
acrylic_p1 <- netSEMp1(acrylic)
acrylic_p2 <- netSEMp2(acrylic, criterion = "AIC") #AIC by default
# Plotting netSEM diagrams
plot(acrylic_p1, cutoff = c(0.3,0.6,0.8))
plot(acrylic_p2, cutoff = c(0.3,0.6,0.8))
```


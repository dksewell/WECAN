
# WECAN

## Overview

The Weighted Edge Clustering Adjusting for Noise (WECAN) model is
designed to estimate edge clusters in weighted networks. Additionally,
this algorithm can identify edges that belong to the “noisy cluster,”
representing anomalous edges.

## Installation

Package aLSEC needs to be installed first.

``` r
install.packages("devtools")

# To install the aLSEC
devtools::install_github('hanhtdpham/aLSEC')

# To install the WECAN
devtools::install_github("HaominLi7/WECAN")
```

## Example of running the algorithm

``` r
library(igraph)
library(WECAN)

# Simulate a network

set.seed(123)

network_data <- wecan_sim(distribution = "Normal", Lambda_sc = 5, 
                                     UVShape=150, UVRate=40, mag=0.4, Ymag=0.9, noise_mean = 0.1,
                                     phi_max = 0.2, noise_prop = 0.15, beta = c(2,7,4,6), conc = 8)[[1]]


############################################################
# Run the WECAN model 15 times, find the ICL for each run.

ICL_temp <- rep(NA, 15)
esti_results <- list()

for (st in 1:15){

    esti_results[[st]] <- wecan(network_data, K = 20, p =  3,
                                 noise_mean = 0.5)
 
    ICL_temp[st] <- ICL(esti_results[[st]])
}


# The one with highest ICL is the final estimation

 est_res <- esti_results[[order(ICL_temp, decreasing = TRUE)[1]]]
 

######################################################
# Alternatively, we can use BIC for model selection
 
BIC_temp <- rep(NA, 15)
esti_results <- list()

for (st in 1:15){

    esti_results[[st]] <- wecan(patients_nwk, K = 20, p =  3,
                                 noise_mean = 0.5)
 
    BIC_temp[st] <- BIC.wecan(esti_results[[st]])
}


# The one with lowest BIC is the final estimation

 est_res <- esti_results[[order(BIC_temp, decreasing = FALSE)[1]]]
 
 
###########################################################
############## Estimation Results #########################
 
esti_clusters <- est_res$estimates$Z_est
esti_clusters <- as.numeric(as.factor(esti_clusters))
table(esti_clusters)
```

Notes that cluster 1 is always the noisy cluster.

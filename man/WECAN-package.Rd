\name{WECAN-package}
\alias{WECAN-package}
\alias{WECAN}
\docType{package}
\title{
  WECAN: Weighted Edge Clustering Adjusting for Noise
}
\description{
  This package implements the WECAN model of Li and Sewell, which clusters the edges of a network.  It uses both edge weights and connectivity patterns to automatically detect superfluous edges, or noise edges, which do not reflect any underlying structure in the network.
}
\details{
  The primary function is the eponymous wecan().  There are also the BIC and ICL functions which can be applied to a wecan object.
}
\author{
Maintainer: Haomin Li <haominli1995@gmail.com>

Contributor: Daniel K. Sewell <daniel-sewell@uiowa.edu>
}
\references{
  Li, Haomin and Sewell, Daniel K. (2025+). Model-based edge clustering for weighted networks with a noise component.  Under review.
}
\examples{
  \dontrun{
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
  }
}

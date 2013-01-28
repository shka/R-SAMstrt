
library(samr)
samr.estimate.depth <- function (xx) {
  x <- xx[which(substr(rownames(xx), 1, 10) == 'RNA_SPIKE_'), ]
  cmeans <- colSums(x)/sum(x)
  depth <- cmeans/mean(cmeans)
  return(depth)
}
assignInNamespace("samr.estimate.depth", samr.estimate.depth, "samr")
resample_norank <- function(x, d, nresamp = 20) {
  ng <- nrow(x)
  ns <- ncol(x)
  dbar <- exp(mean(log(d)))
  xresamp <- array(0, dim = c(ng, ns, nresamp))
  for (k in 1:nresamp) {
    for (j in 1:ns) { xresamp[, j, k] <- rpois(n = ng, lambda = (dbar/d[j]) * x[, j]) + runif(ng) * 0.1 }
  }
  return(xresamp)
}
mean_xresamp0 <- function(xresamp, nresamp) {
  if(nresamp == 1) xresamp[, , 1]
  else mean_xresamp0(xresamp, nresamp-1) + xresamp[, , nresamp]
}
mean_xresamp <- function(xresamp, nresamp=20) { mean_xresamp0(xresamp, nresamp)/nresamp }
SAMstrt.normalization <- function(reads, nresamp=20, spikesum=2234.8) {
  nreads.org <- mean_xresamp(resample_norank(reads, samr.estimate.depth(reads), nresamp), nresamp)
  rownames(nreads.org) <- rownames(reads)
  colnames(nreads.org) <- colnames(reads)
  spikes <- colSums(nreads.org[which(substr(rownames(nreads.org), 1, 10) == 'RNA_SPIKE_'), ])
  nreads.org/rep(spikes, each=nrow(nreads.org))*spikesum
}

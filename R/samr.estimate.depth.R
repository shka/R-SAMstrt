samr.estimate.depth <-
function (xx) {
  x <- xx[which(substr(rownames(xx), 1, 10) == 'RNA_SPIKE_'), ]
  cmeans <- colSums(x)/sum(x)
  depth <- cmeans/mean(cmeans)
  return(depth)
}

.onLoad <- function(libname, pkgname) {
  if (pkgname == packageName()) {
    assignInNamespace("samr.estimate.depth", samr.estimate.depth, "samr")
    getAnywhere(samr.estimate.depth)
  }
}

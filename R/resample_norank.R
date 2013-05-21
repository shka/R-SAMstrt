resample_norank <-
function(x, d, nresamp = 20) {
  ng <- nrow(x)
  ns <- ncol(x)
  dbar <- exp(mean(log(d)))
  xresamp <- array(0, dim = c(ng, ns, nresamp))
  for (k in 1:nresamp) {
    for (j in 1:ns) { xresamp[, j, k] <- rpois(n = ng, lambda = (dbar/d[j]) * x[, j]) + runif(ng) * 0.1 }
  }
  return(xresamp)
}

\name{SAMstrt-package}
\alias{SAMstrt-package}
\alias{SAMstrt}
\docType{package}
\title{Statistical test for differential expression using sequencing-based methods with spike-in normalization.}
\description{
  The SAMstrt package provides an significance analysis of sequencing data
  with spike-in normalization. For the normalization method, same amount
  of spike-in molecules must be added into the target samples, for
  example by Islam2011. The statistical backgrounds and the benefits depend
  on SAMseq of the samr package (Li2011); however one exception is
  estimation of sequencing depth to make spike-in read-counts equivalent
  over the samples. 
}
\details{
\tabular{ll}{
Package: \tab SAMstrt\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2013-01-28\cr
License: \tab LGPL-3\cr
}
  Loading of this package change the original normalization method of
  SAMseq to the spike-in based, however all usages are same with
  the original SAMseq. Function SAMstrt.normalization in this package
  provides the normalized values for visualization and the other
  purposes.
}
\author{Maintainer: Shintaro Katayama <shintaro.katayama@gmail.com>}
\references{
  Islam S, et al. Highly multiplexed and strand-specific single-cell RNA
  5'-end sequencing. Nat Protoc. 2012 May;7(5):813-28.
  
  Katayama S, et al. SAMstrt: Statistical test for differential
  expression in single-cell transcriptome with spike-in
  normalization. Submitted.

  Li J, Tibshirani R. Finding consistent patterns: A nonparametric
  approach for identifying differential expression in RNA-Seq data. Stat
  Methods Med Res. 2011 Nov 28.
}
\keyword{package}
\keyword{nonparametric}
\seealso{
  \code{\link[samr:SAMseq]{SAMseq}}, \code{\link[SAMstrt:SAMstrt.normalization]{SAMstrt.normalization}}
}
\examples{
set.seed(100)
mu <- matrix(100, 1000, 20)
mu[1:100, 11:20] <- 200
mu <- scale(mu, center=FALSE, scale=runif(20, 0.5, 1.5))
x <- matrix(rpois(length(mu), mu), 1000, 20)
rownames(x) <- sapply(1:nrow(x),
  function(i){ ifelse(i<4, sprintf("RNA_SPIKE_\%d",i),
                                sprintf("FEATURE_\%d",i-3)) })
y <- c(rep(1, 10), rep(2, 10))
samfit <- SAMseq(x, y, resp.type = "Two class unpaired")
print(samfit)
}

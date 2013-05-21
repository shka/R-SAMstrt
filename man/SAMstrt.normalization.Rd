\name{SAMstrt.normalization}
\alias{SAMstrt.normalization}
\title{Spike-in based normalization}
\description{
  This function normalizes given raw read-counts by "spike-in based
  normalization"; it corrects differences in counts of sequenced
  spike-in reads.
}
\usage{ SAMstrt.normalization(reads, nresamp = 20, spikesum = 2234.8) }
\arguments{
  \item{reads}{ Matrix of read-counts, p (number of features) by n
    (number of samples). The values must be integral numbers. }
  \item{nresamp}{ Number of resamples used to construct test statistic;
    normalized values are average of the resampling. Default is 20. }
  \item{spikesum}{ Expected total number of spike-in molecules in each
    sample; normalized values are estimated molecular counts by the
    expectation. Default is 2234.8. }
}
\details{
  This function performs Poisson resamplings to make spike-in
  read-counts equivalent over the samples, and returns average of the
  resampled values as normalized read-counts. In detail, sum of
  spike-in read-counts in sample $i$ is $d_i$ of the formula 2.4 in
  Li2011. Features, which name begins with "RNA_SPIKE_", are for
  the spike-ins.}
\value{ Matrix of normalized read-counts. }
\references{
  Katayama S, et al. SAMstrt: Statistical test for differential
  expression in single-cell transcriptome with spike-in
  normalization. Submitted.

  Li J, Tibshirani R. Finding consistent patterns: A nonparametric
  approach for identifying differential expression in RNA-Seq data. Stat
  Methods Med Res. 2011 Nov 28.
}
\author{Shintaro Katayama}
\examples{
set.seed(100)
mu <- matrix(100, 1000, 20)
mu[1:100, 11:20] <- 200
mu <- scale(mu, center=FALSE, scale=runif(20, 0.5, 1.5))
x <- matrix(rpois(length(mu), mu), 1000, 20)
rownames(x) <- sapply(1:nrow(x),
  function(i){ ifelse(i<4, sprintf("RNA_SPIKE_\%d",i),
				sprintf("FEATURE_\%d",i-3)) })
x.norm <- SAMstrt.normalization(x)
}
\keyword{nonparametric}

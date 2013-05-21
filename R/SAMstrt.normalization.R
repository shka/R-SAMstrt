SAMstrt.normalization <-
function(reads, nresamp=20, spikesum=2234.8) {
  nreads.org <- mean_xresamp(resample_norank(reads, samr.estimate.depth(reads), nresamp), nresamp)
  rownames(nreads.org) <- rownames(reads)
  colnames(nreads.org) <- colnames(reads)
  spikes <- colSums(nreads.org[which(substr(rownames(nreads.org), 1, 10) == 'RNA_SPIKE_'), ])
  nreads.org/rep(spikes, each=nrow(nreads.org))*spikesum
}

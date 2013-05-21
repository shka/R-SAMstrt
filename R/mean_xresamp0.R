mean_xresamp0 <-
function(xresamp, nresamp) {
  if(nresamp == 1) xresamp[, , 1]
  else mean_xresamp0(xresamp, nresamp-1) + xresamp[, , nresamp]
}

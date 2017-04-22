#' Scan through suggested pairs of neurons
#'
#' @description implementation of the algorithm for clustering neurons by synapse location from Schlegel et al. (2016). Assumes neurons are scaled to microns.
#'
#' @param someneuronlist a neuronlist or neuron object
#' @param pairs A data frame / matrix of two columns named skid.rgight and skid.left perhaps generated using the deformetricar package to estimate neuron cognates.
#' @param reference A reference object to plot. E.g. CNS cortex.
#' @param ... additional arguments passed to methods.
#'
#' @export
scan4matching <- function(someneuronlist, pairs, reference, ...){
  open3d(userMatrix = structure(c(0.999989449977875, -0.00419542612507939, -0.00183705519884825,
                                  0, -0.00426193978637457, -0.999274253845215, -0.0378473773598671,
                                  0, -0.00167691614478827, 0.037854727357626, -0.999281764030457,
                                  0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.5050682, windowRect = c(1580L, 136L, 2768L, 1051L))
  k = c()
  m = c()
  d = c()
  progress = 'return'
  for (row in 1:nrow(pairs)){
    rgl::points3d(xyzmatrix(reference), col = 'grey')
    rgl::plot3d(someneuronlist[c(as.character(c(pairs$skid.right)[row]), as.character(c(pairs$skid.left)[row]))] )
    print(paste(row, "/", nrow(pairs), sep = ""))
    while (progress == 'return'){
      progress = readline(prompt="Don't keep (d), maybe keep (m), or keep (k)? ")
      if (progress == 'k'){
        k = rbind(k, pairs[row,])
      } else if (progress == 'm'){
        m = rbind(m, pairs[row,])
      } else if (progress == 'd'){
        d = rbind(d, pairs[row,])
      } else{
        progress = 'return'
      }
    }
    progress = 'return'
    rgl::clear3d()
  }
  l = list(k, m, d)
  names(l) <- c("keep", "maybe", "loose")
  return(l)
}




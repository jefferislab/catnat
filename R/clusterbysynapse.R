#' Cluster neurons by pre- and postsynapse positons
#'
#' @description implementation of the algorithm for clustering neurons by synapse location from Schlegel et al. (2016). Assumes neurons are scaled to microns.
#'
#' @param someneuronlist a neuronlist or neuron object
#' @param sigma determines what distances between two synapses are considered close (defaults to 1 um)
#' @param omega synapse cluster radius. Defaults to sigma.
#' @param symmetric whether to return a symmetric martrix (average of scores between two neurons in both directions)
#' @param direction cluster neurons by presynapses (0), postysnapses (1) or both (2)
#' @param ... additional arguments passed to methods.
#'
#' @return A matrix of similarity scores between inputted neurons, based on synapse positions.
#' @export
#' @seealso \code{\link{seesplit3d}} \code{\link{get.synapses}}
clusterbysynapses <- function(someneuronlist, sigma = 1, omega = sigma, symmetric = T, direction = c(0,1,2)){
  m = matrix(nrow = length(someneuronlist), ncol = length(someneuronlist))
  colnames(m) = rownames(m) = names(someneuronlist)
  for (neuron in 1:length(someneuronlist)){
    print(paste(neuron,"/",length(someneuronlist)))
    g = as.data.frame(get.synapses(someneuronlist[neuron], "BOTH"))
    g[g[,4] == 0,4] <- 2
    if (direction > 0){g = g[g[,4] == direction,]}
    for (neuron2 in 1:length(someneuronlist)){
      t = as.data.frame(get.synapses(someneuronlist[neuron2], "BOTH"))
      t[t[,4] == 0,4] <- 2
      if (direction > 0){t = t[t[,4] == direction,]}
      scores = matrix(ncol = 2, nrow = nrow(g))
      for (syn in 1:nrow(g)){
        gg = subset(g, g$prepost == g[syn,"prepost"])[,-4]
        tt = subset(t, t$prepost == g[syn,"prepost"])[,-4]
        if(plyr::empty(tt)){ scores[syn,g[syn,"prepost"]] = 0; break}
        n = nabor::knn(tt, nat::xyzmatrix(g[syn,]), k =1)
        close = t[n$nn.idx,]
        gn = nabor::knn(nat::xyzmatrix(g[syn,]), gg, k =1)
        gn = sum(gn$nn.dists<omega) - 1
        tn = nabor::knn(nat::xyzmatrix(close), tt, k =1)
        tn = sum(tn$nn.dists<omega) - 1
        multiplier = exp(-abs(gn-tn)/(tn+gn))
        if(is.infinite(multiplier)|is.na(multiplier)|is.nan(multiplier)){ multiplier = 0}
        score = exp((-n$nn.dist^2)/2*(sigma)^2)*multiplier
        scores[syn,g[syn,"prepost"]] = score
      }
      m[neuron,neuron2] = ifelse(direction[1] >0, max(colSums(scores, na.rm = T)),min(colSums(scores, na.rm = T)))
    }
  }
  if (symmetric == T){
    #pmean <- function(x,y) (x+y)/2 take the mean, or take the lowest score as below.
    m[] <- pmin(m, matrix(m, nrow(m), byrow=TRUE))
  }
  return(m)
}


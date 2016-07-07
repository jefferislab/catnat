#' Cluster neurons by pre- and postsynapse positons
#'
#' @description implementation of the algorithm for clustering neurons by synapse location from Schlegel et al. (2016). Assumes neurons are scaled to microns.
#'
#' @param someneuronlist a neuronlist or neuron object
#' @param sigma determines what distances between two synapses are considered close (defaults to 2 um)
#' @param omega synapse cluster radius. Defaults to sigma.
#' @param symmetric whether to return a symmetric martrix (average of scores between two neurons in both directions)
#' @param ... additional arguments passed to methods.
#'
#' @details From Schneider-Mizell et al. (2016): "We use flow centrality for four purposes. First, to split an arbor into axon and dendrite at the maximum centrifugal SFC, which is a preliminary step for computing the segregation index, for expressing all kinds of connectivity edges (e.g. axo-axonic, dendro-dendritic) in the wiring diagram, or for rendering the arbor in 3d with differently colored regions. Second, to quantitatively estimate the cable distance between the axon terminals and dendritic arbor by measuring the amount of cable with the maximum centrifugal SFC value. Third, to measure the cable length of the main den- dritic shafts using centripetal SFC, which applies only to insect neurons with at least one output syn- apse in their dendritic arbor. And fourth, to weigh the color of each skeleton node in a 3d view, providing a characteristic signature of the arbor that enables subjective evaluation of its identity."
#'
#' @return A matrix of similarity scores between inputted neurons, based on synapse positions.
#'
#' @references Schlegel, P., Texada, M. J., Miroschnikow, A., Peters, M., Schneider-Mizell, C. M., Lacin, H., … Pankratz, M. J. (2016). Synaptic Transmission Parallels Neuromodulation in a Central Food-Intake Circuit. bioRxiv, 1–26. Retrieved from http://www.biorxiv.org/content/early/2016/04/07/044990.abstract
#'
#' @export
#' @rdname cluster.by.synapses
#' @seealso \code{\link{plot3d.split}} \code{\link{get.synapses}}
cluster.by.synapses <- function(someneuronlist, sigma = 2, omega = sigma, symmetric = T, ...){
  m = matrix(nrow = length(someneuronlist), ncol = length(someneuronlist))
  colnames(m) = rownames(m) = names(someneuronlist)
  for (neuron in 1:length(someneuronlist)){
    g = get.synapses(someneuronlist[neuron], "BOTH")
    for (neuron2 in 1:length(someneuronlist)){
      t = get.synapses(someneuronlist[neuron2], "BOTH")
      scores = c()
      for (syn in 1:nrow(g)){
        gg = subset(g, g$prepost == g[syn,"prepost"])[,-4]
        tt = subset(t, t$prepost == g[syn,"prepost"])[,-4]
        if(empty(tt)){ score = 0; break}
        n = nabor::knn(tt, nat::xyzmatrix(g[syn,]), k =1)
        close = t[n$nn.idx,]
        gn = nabor::knn(nat::xyzmatrix(g[syn,]), gg, k =1)
        gn = sum(gn$nn.dists<omega) - 1
        tn = nabor::knn(nat::xyzmatrix(close), tt, k =1)
        tn = sum(tn$nn.dists<omega) - 1
        multiplier = exp(abs(gn-tn)/(tn+gn))
        if(is.infinite(multiplier)|is.na(multiplier)|is.nan(multiplier)){ multiplier = 0}
        score = exp((-n$nn.dist^2)/2*(sigma)^2)*multiplier
        scores = c(scores,score)
      }
      m[neuron,neuron2] = mean(scores)
    }
  }
  if (symmetric == T){
    pmean <- function(x,y) (x+y)/2
    m[] <- pmean(m, matrix(m, nrow(s), byrow=TRUE))
  }
  return(m)
}

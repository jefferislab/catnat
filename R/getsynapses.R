#' Get 3D coordinates of synapse positions
#'
#' @param someneuronlist a neuronlist or neuron object
#' @param target whether post or presynapse are ot be returned, or both
#' @param polypre whether to consider the number of presynapses as a multiple of the numbers of connections each makes
#' @param ... additional arguments passed to methods.
#'
#' @return Anatomically accurate synapse position (i.e. not just connector positions) as 3D coordinates
#' @export
#' @rdname get.synapses
#' @seealso \code{\link{cluster.by.synapses}} \code{\link{flow.centrality}}
get.synapses <-function(someneuronlist, target = c("BOTH", "PRE", "POST"), polypre = T, ...) UseMethod("get.synapses")

#' @export
#' @rdname get.synapses
get.synapses.neuron <- function (someneuronlist, target = target, polypre = polypre){
  if (target%in%c("POST","BOTH")) {
    syns.in = neuron$connectors[neuron$connectors[,3]==1,][,1]
    point.no = rownames(neuron$d)[match(syns.in,neuron$d[,"PointNo"])]
    points = nat::xyzmatrix(neuron$d[point.no,])
    points = cbind(points, prepost = 1)
  }
  if (target%in%c("PRE","BOTH")) {
    pres = neuron$connectors[neuron$connectors[,3]==0,][,2]
    pre.cons = catmaid_get_connectors(pres)$connector_id
    syns.out = neuron$connectors[,1][match(pre.cons, neuron$connectors[,2])]
    point.no = rownames(neuron$d)[match(syns.out,neuron$d[,"PointNo"])]
    points.out = nat::xyzmatrix(neuron$d[point.no,])
    points.out = cbind(points.out, prepost = 0)
    if(target == "BOTH"){
      points = rbind(points, points.out)
    }else{ points = points.out}
  }
  points
}

#' @export
#' @rdname get.synapses
get.synapses.neuronlist <- function (someneuronlist, target = c("BOTH", "PRE", "POST"), polypre = T){
  points = nat::nlapply(someneuronlist, get.synapses.neuron, target = target, polypre = polypre)
  points = do.call(rbind, points)
  points
}

#' Prune nodes from a neuron, keeping the root node
#'
#' @description Remove points from a neuron, keeping the root node intact
#'
#' @param neuron a neuron object
#' @param bad_vertex_labels Nodes ids for removal
#' @param ... additional arguments passed to methods.
#' @return A pruned neuron object
#' @export
#' @rdname prune_distal
prune_distal <- function(neuron, bad_vertex_labels, invert=FALSE, ...) {
  # remove points **keeping same origin**
  origin=n$d$PointNo[rootpoints(neuron)]
  if(origin%in%bad_vertex_labels)
    n1=subset(neuron, !PointNo%in%bad_vertex_labels)
  else
    n1=subset(neuron, !PointNo%in%bad_vertex_labels, origin=origin)
  prune_vertices(n1, unlist(as.seglist(n1, all=F)), invert =! invert)
}

#' Find how many downstream partners tracers have connected up
#'
#' @description Remove points from a neuron, keeping the root node intact
#'
#' @param someneuron a neuronlist object. Should be the downstream partners of (a) neuron(s) of interest
#' @param names Full names of tracers to consider
#' @param ... additional arguments passed to methods.
#' @return The percentage of downstream singlet nodes named tracers have contributed to neuron in
#' @export
#' @rdname downstream.deletion.test
downstream.deletion.test <- function(someneuronlist,names = c("Alex Bates", "Ruairi Roberts"), ...){
  neurons.treenodes=nlapply(names(someneuronlist), catmaid_get_treenode_table)
  user_ids = subset(ul, full_name%in%names)$id
  neurons.treenodes %>%
    bind_rows %>%
    filter(user_id%in%user_ids) ->
    neurons.treenodes.names
  neurons.treenodes %>%
    bind_rows %>%
    filter(user_id!=66) ->
    neurons.treenodes.others
  # now we remove all the nodes that named tracers traced from those neurons
  someneuronlist.pruned=nlapply(someneuronlist, prune_distal, neurons.treenodes.names$id)
  treenode_ids=unlist(sapply(someneuronlist.pruned, function(x) x$d$PointNo))
  someneuronlist.names=nlapply(someneuronlist, subset, !PointNo%in%treenode_ids, OmitFailures = T)
  # and then we get the tree node ids that are left
  treenodes_names=unlist(sapply(someneuronlist.names, function(x) x$d$PointNo))
  treenodes_nonames=unlist(sapply(someneuronlist.pruned, function(x) x$d$PointNo))
  subset(cneurons, post_node_id %in% treenodes_nonames) -> synapses_nonames
  100 - nrow(synapses_nonames) / nrow(cneurons) *100
}


#' Prune nodes from a catmaid neuron, keeping the synapses
#'
#' @description Prune nodes from a catmaid neuron, keeping the synapses
#'
#' @param neuron a neuron object
#' @param target Nodes ids for removal
#' @param maxdist The threshold distance for keeping points
#' @param keep Whether to keep points in x that are near or far from the target
#' @param return.indices Whether to return the indices that pass the test rather than the 3D object/points (default FALSE)
#' @param ... additional arguments passed to methods.
#' @return A pruned neuron object
#' @export
#' @aliases prune
#' @importFrom nat prune
prune.catmaidneuron<- function (neuron,target,maxdist, keep = c("near", "far"),
                                return.indices = FALSE,...){
  pruned = prune(neuron,target,maxdist=maxdist, keep = keep,
                 return.indices = return.indices)
  pruned$connectors = neuron$connectors[neuron$connectors$treenode_id%in%pruned$d$PointNo,]
  pruned
}


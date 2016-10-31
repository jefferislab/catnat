#' Segregate neurite compartments from a neuron/neuronlist
#'
#' @description Fragment a skeleton into compartments that have been assigned using flow.centrality(). Note that this will break down a neuron regardless of the segregation score calculated by the flow.centrality algorithm.
#'
#' @param someneuronlist A neuronlist that has been processed by flow.centrality()
#' @param fragment The type of neurite fragment to retrieve. 'Nulls' refers to areas of zero flow. See Schnieder-Mizell et al. (2016)
#' @param ... additional arguments passed to methods.
#'
#' @return Neurites as a neuronlist object, complete with  synaptic information relevant to that fragment.
#' @export
#' @rdname neurites
#' @seealso \code{\link{get_connected_skeletons}} \code{\link{skeleton_connectivity_matrix}} \code{\link{flow.centrality}}
neurites <-function(someneuronlist, fragment = c("axons","dendrites","primary dendrite","primary neurite","nulls"), ...) UseMethod("neurites")

#' @export
#' @rdname neurites
neurites.neuron <- function(neuron, fragment, ...){
  if (is.null(neuron$d$flow.cent)) {
    warning("No flow centrality calculated, dropping neuron")
    break
  }
  dendrites.v = subset(rownames(neuron$d), neuron$d$compartment ==  "dendrite")
  axon.v = subset(rownames(neuron$d), neuron$d$compartment == "axon")
  nulls.v = subset(rownames(neuron$d), neuron$d$compartment == "null")
  p.d.v = subset(rownames(neuron$d), neuron$d$compartment == "primary dendrite")
  p.n.v = subset(rownames(neuron$d), neuron$d$compartment == "primary neurite")
  if (fragment == "dendrites"){
    tree = nat::prune_vertices(neuron, verticestoprune = as.integer(c(axon.v,nulls.v, p.d.v, p.n.v)))
    tree$connectors = neuron$connectors[neuron$connectors$treenode_id%in%tree$d$PointNo,]
  }
  else if (fragment == "axons"){
    tree = nat::prune_vertices(neuron, verticestoprune = as.integer(c(nulls.v, dendrites.v, p.d.v, p.n.v)))
    tree$connectors = neuron$connectors[neuron$connectors$treenode_id%in%tree$d$PointNo,]
  }
  else if (fragment == "nulls"){
    tree = nat::prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, dendrites.v, p.d.v, p.n.v)))
  }
  else if (fragment == "primary dendrite"){
    tree = prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, dendrites.v, nulls.v, p.n.v)))
    tree$connectors = neuron$connectors[neuron$connectors$treenode_id%in%tree$d$PointNo,]
  }
  else if (fragment == "primary neurite"){
    tree = prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, dendrites.v, nulls.v, p.d.v)))
    tree$connectors = neuron$connectors[neuron$connectors$treenode_id%in%tree$d$PointNo,]
  }
  tree
}

#' @export
#' @rdname neurites
neurites.neuronlist <- function(someneuronlist, fragment, ... ){
  trees = nat::nlapply(someneuronlist, neurites, fragment, OmitFailures = T)
  trees[,"name"] = unlist(lapply(trees[,"name"], function(x) paste(x,'#',fragment,sep='')))
  names(trees) = unlist(lapply(names(trees), function(x) paste(x,'#',fragment,sep='')))
  trees
}

#' Segregate neurite compartments from a neuron/neuronlist based on synapse clustering
#'
#' @description Fragment a skeleton into compartments that have been assigned cluster_synapses_within_skeleton() or flow.centrality(). Note that this will break down a neuron regardless of the segregation score calculated across the arbour.
#'
#' @param someneuronlist A neuronlist that has been processed by flow.centrality()
#' @param arbourlist A neuronlist produced from arbour.custers
#' @param ... additional arguments passed to methods.
#'
#' @return Segmented arbours as a neuronlist object, complete with  synaptic information relevant to that fragment. If a neuronlist is given, subsequently plotting the neuronlist will not reveal the fragments. However, plotting individual neurons double indexed in the lists will. arbours() returns a neurinlist, where each entry is an arbour fragment with a unique skid and name.
#' @seealso \code{\link{neurites}} \code{\link{cluster_synapses_within_skeleton}} \code{\link{seebroken3d}} \code{\link{flow.centrality}} \code{\link{seesplit3d}}
#' @export
#' @rdname arbour.clusters
arbour.clusters <-function(someneuronlist, ...) UseMethod("arbour.clusters")

#' @export
#' @rdname arbour.clusters
arbour.clusters.neuron <- function(someneuron, ...){
  if (!is.null(someneuron$d$cluster)){
    clusters = unique(someneuron$d$cluster)
    entropies = someneuron$d$cluster.entropy[!duplicated(someneuron$d$cluster)]
  }else if (!is.null(neuron$d$flow.cent)) {
    clusters = unique(someneuron$d$compartment)
    entropies = someneuron$d$cluster.entropy[!duplicated(someneuron$d$flow.cent)]
  }else{
    warning("No clustering calculated, dropping neuron")
    break
  }
  arbours = neuronlist()
  for (c in 1:length(clusters)){
    cluster = clusters[[c]]
    cluster.v = subset(rownames(someneuron$d), someneuron$d$cluster == cluster)
    cluster.g = nat::prune_vertices(someneuron, verticestoprune = rownames(someneuron$d)[!rownames(someneuron$d)%in%cluster.v])
    cluster.g$connectors = someneuron$connectors[someneuron$connectors$treenode_id%in%cluster.g$d$PointNo,]
    cluster.g$d$value = entropies[c]
    if (!is.null(cluster.g$cluster.segregation.index)){cluster.g$cluster.segregation.index = someneuron$cluster.segregation.index}
    if (!is.null(cluster.g$AD.segregation.index)){cluster.g$cluster.segregation.index = someneuron$AD.segregation.index}
    arbours = c(arbours, nat::as.neuronlist(cluster.g))
  }
  # Modify name
  names(arbours) = clusters
  arbours
}

#' @export
#' @rdname arbour.clusters
arbour.clusters.neuronlist <- function(someneuronlist, ...){
  arbours = nat::nlapply(someneuronlist, arbour.clusters.neuron, OmitFailures = T)
  # work out names
  # For some reason, plotting synapses WithConnectors does not work with these
  arbours
}

#' @export
#' @rdname arbour.clusters
arbours.neuron <- function(arbourcluster, ...){
  skeletons <- neuronlist()
  x = neuronlist()
  for (n in length(arbourcluster)){
    oldname = arbourcluster[,"name"][n]
    skid = name = arbourcluster[,"skid"][n]
    bits = nlapply(arbourcluster[[n]], as.neuron)
    newnames = unlist(lapply(names(bits), function(x) paste(oldname,'#',x,sep='')))
    names(bits) = unlist(lapply(names(bits), function(x) paste(skid,'#',x,sep='')))
    dfn = data.frame(pid = 1, skid = skid, name = newnames)
    rownames(dfn) = names(bits)
    attr(bits,'df') = dfn
    skeletons <- c(skeletons, bits)
  }
  skeletons
}

#' @export
#' @rdname arbour.clusters
arbours <- function(arbourcluster, ...){
  if (length(arbourcluster)>1){
    skeletons <- nlapply(arbourcluster, arbours, OmitFailures = T)
  }else{skeletons = arbours.neuron(arbourcluster)}
  skeletons
}




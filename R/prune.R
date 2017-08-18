#' Prune nodes from a neuron, keeping the root node
#'
#' @description Remove points from a neuron, keeping the root node intact
#'
#' @param neuron a neuron object
#' @param bad_vertex_labels Nodes ids for removal
#' @param invert whether to keep vertices rather than dropping them (default FALSE)
#' @param ... additional arguments passed to methods.
#' @return A pruned neuron object
#' @export
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
#' @param someneuronlist a neuronlist object. Should be the downstream partners of (a) neuron(s) of interest
#' @param names Full names of tracers to consider
#' @param ... additional arguments passed to methods.
#' @return The percentage of downstream singlet nodes named tracers have contributed to neuron in
#' @export
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
#' @param x a neuron object
#' @param target Nodes ids for removal
#' @param maxdist The threshold distance for keeping points
#' @param keep Whether to keep points in x that are near or far from the target
#' @param return.indices Whether to return the indices that pass the test rather
#'   than the 3D object/points (default FALSE)
#' @param ... additional arguments passed to methods (i.e.
#'   \code{\link[nat]{prune}}).
#' @return A pruned neuron object
#' @export
#' @aliases prune
#' @importFrom nat prune
#' @seealso \code{\link[nat]{prune}}
prune.catmaidneuron<- function (x,target,maxdist, keep = c("near", "far"),
                                return.indices = FALSE,...){
  class(x) = c("neuron")
  pruned = nat::prune(x,target=target, maxdist=maxdist, keep = keep,
                 return.indices = return.indices, ...)
  pruned$connectors = x$connectors[x$connectors$treenode_id%in%pruned$d$PointNo,]
  relevant.points = subset(x$d, PointNo%in%pruned$d$PointNo)
  y = pruned
  y$d = relevant.points[match(pruned$d$PointNo,relevant.points$PointNo),]
  y$d$Parent = pruned$d$Parent
  y
}

#' Prune vertices from a CATMAID neuron, keeping the synapses
#'
#' @description Prune nodes from a catmaid neuron, keeping the synapses
#'
#' @param x a CATMAID neuron object
#' @inheritParams nat::prune_vertices
#' @param ... additional arguments passed to methods
#' @return A pruned neuron object
#' @export
#' @rdname prune_vertices.catmaidneuron
prune_vertices.catmaidneuron<- function (x,verticestoprune, invert = FALSE,...){
  class(x) = c("neuron")
  pruned = nat::prune_vertices(x,verticestoprune,invert = invert)
  pruned$connectors = x$connectors[x$connectors$treenode_id%in%pruned$d$PointNo,]
  relevant.points = subset(x$d, PointNo%in%pruned$d$PointNo)
  y = pruned
  y$d = relevant.points[match(pruned$d$PointNo,relevant.points$PointNo),]
  y$d$Parent = pruned$d$Parent
  pruned
}

#' Prune a neuron interactively
#'
#' @description Remove points from a neuron, keeping the root node intact
#'
#' @param x a neuron/neuronlist object
#' @param ... additional arguments passed to methods
#' @return A pruned neuron/neuronlist object
#' @export
#' @rdname prune_online
prune_online <-function(x, ...) UseMethod("prune_online")

#' @export
#' @rdname prune_online
prune_online.neuron <- function(x, ...){
  continue = "no"
  while(!continue%in%c("y","yes")){
    selected = select.points(nat::xyzmatrix(x), plot3d = x)
    neuron = nat::prune(x, target = selected, keep = "near", maxdist = 0)
    rgl::plot3d(neuron, col ="black")
    continue = readline("Finished with this neuron? yes/no ")
  }
  neuron
}

#' @export
#' @rdname prune_online
prune_online.neuronlist <- function(x, ...){
  nat::nlapply(x,prune_online.neuron)
}

#' Manually assign the dendrite and axon to a neuron
#'
#' @description Manually assign the dendrite and axon to neurons / a neuron
#'
#' @param x a neuron/neuronlist object
#' @param ... additional arguments passed to methods
#' @return The neuron/neuronlist object with axon/dendrite info assigned in SWC format to neuron$d
#' @export
#' @rdname manually_assign_axon_dendrite
manually_assign_axon_dendrite <-function(x, ...) UseMethod("manually_assign_axon_dendrite")

#' @export
#' @rdname manually_assign_axon_dendrite
manually_assign_axon_dendrite.neuron <- function(x){
  happy = "no"
  skid = x$skid
  x$d$Label = 0
  while(!happy%in%c("y","yes")){
    rgl::clear3d()
    message("Please choose dendrites for your neuron")
    dend = prune_online.neuron(x)
    if(!"catmaidneuron"%in%class(x)){
      x$d$Label[x$d$X%in%dend$d$X&x$d$Y%in%dend$d$Y] = 3
    }else{
      x$d$Label[x$d$PointNo%in%dend$d$PointNo] = 3
    }
    rgl::clear3d()
    message("Please choose axon for your neuron")
    axon = prune_online.neuron(x)
    if(!"catmaidneuron"%in%class(x)){
      x$d$Label[x$d$X%in%axon$d$X&x$d$Y%in%axon$d$Y] = 2
    }else{
      x$d$Label[x$d$PointNo%in%axon$d$PointNo] = 2
    }
    x$d$Label[nat::rootpoints(x)] = 1
    rgl::clear3d()
    rgl::plot3d(dend,col="blue")
    rgl::plot3d(axon, col = "orange")
    rgl::plot3d(x, col = "purple")
    happy = readline("Happy with this division? yes/no  ")
  }
  x$skid = skid
  x
}

#' @export
#' @rdname manually_assign_axon_dendrite
manually_assign_axon_dendrite.neuronlist<-function(x){
  nat::nlapply(x, manually_assign_axon_dendrite.neuron)
}

#' Give connector data in a CATMAID neuron the same attributes as node data
#'
#' @description Give connector data in a CATMAID neuron the same attributes as
#'   node data. I.e. adding Label information to indicate compartments such as
#'   axon and dendrite
#'
#' @param x a neuron/neuronlist object that has primary neurites marked (Label =
#'   7) and soma as the root
#' @param ... Additional arguments passed to nlapply
#' @export
#' @rdname assign.connector.info
assign.connector.info <-function(x, ...) UseMethod("assign.connector.info")

assign.connector.info.neuron<-function(x){
  relevant.points = subset(x$d, PointNo%in%x$connectors$treenode_id)
  x$connectors = cbind(x$connectors,relevant.points[match(x$connectors$treenode_id,relevant.points$PointNo),colnames(relevant.points)[!colnames(relevant.points)%in%c("PointNo", "Label", "X", "Y", "Z", "W", "Parent")]])
  x
}
assign.connector.info.neuronlist<-function(x, ...){
  nlapply(x,assign.connector.info.neuron, ...)
}

#' Prune neuron within a mesh3d volume
#'
#' @description Give connector data in a CATMAID neuron the same attributes as
#'   node data. I.e. adding Label information to indicate compartments such as
#'   axon and dendrite
#'
#' @param x a neuron/neuronlist object that has primary neurites marked (Label =
#'   7) and soma as the root
#' @param brain The \code{\link[nat]{hxsurf}} object containing the neuropil of
#'   interest, e.g. \code{\link[nat.flybrains]{FCWBNP.surf}}
#' @param neuropil Character vector specifying the neuropil
#' @param invert Logical when \code{TRUE} indicating that points outside the
#'   surface should be pruned.
#' @inheritParams nat::prune
#' @export
prune.in.volume<- function(x, brain, neuropil = "LH_R", maxdist = 0, invert = FALSE){
  if(invert){
    keep = "far"
  }else{
    keep = "near"
  }
  nat::prune(x, nat::xyzmatrix(x)[nat::pointsinside(nat::xyzmatrix(pns),surf=nat::as.mesh3d(subset(brain,"LH_R"))),],maxdist=maxdist,keep = keep)
}


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
#' @param return.indices Whether to return the indices that pass the test rather than the 3D object/points (default FALSE)
#' @param verticestoprune	an integer vector describing which vertices to remove
#' @param invert	whether to keep vertices rather than dropping them (default FALSE)
#' @param ... additional arguments passed to methods
#' @return A pruned neuron object
#' @export
#' @aliases prune
#' @importFrom nat prune
prune.catmaidneuron<- function (x,target,maxdist, keep = c("near", "far"),
                                return.indices = FALSE,...){
  class(x) = c("neuron")
  pruned = prune(x,target,maxdist=maxdist, keep = keep,
                 return.indices = return.indices)
  pruned$connectors = x$connectors[x$connectors$treenode_id%in%pruned$d$PointNo,]
  pruned$skid = x$skid
  pruned
}


#' @aliases prune_vertices
#' @importFrom nat prune_vertices
prune_vertices.catmaidneuron<- function (x,verticestoprune, invert = FALSE,...){
  class(x) = c("neuron")
  pruned = prune_vertices(x,verticestoprune,invert = invert)
  pruned$connectors = x$connectors[x$connectors$treenode_id%in%pruned$d$PointNo,]
  pruned$skid = x$skid
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
    selected = catnat:::select.points(nat::xyzmatrix(x), plot3d = x)
    neuron = nat::prune(x, target = selected, keep = "near", maxdist = 0)
    rgl::plot3d(neuron, col ="black")
    continue = readline("Finished with this neuron? yes/no ")
  }
  neuron$skid = x$skid
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


#' Functions to assign and visualise microtubule rich and twig portions of a neuron
#'
#' @description Manually assign the dendrite and axon to neurons / a neuron
#'
#' @param x a neuron/neuronlist object
#' @param microtubules whether to return the microtubule containing arbour (TRUE) or twigs (FALSE)
#' @return The neuron/neuronlist object with axon/dendrite info assigned in SWC format to neuron$d
#' @export
#' @rdname microtubules
mark.microtubules <- function(x){
  if(is.null(x$d$microtubules)){
    if(is.null(x$tags$`microtubules end`)){
      message("No microtubular endings marked in CATMAID neuron")
      break
    }
    root = nat::rootpoints(x)
    microtubule.endings.pointno = x$tags$`microtubules end`
    microtubule.endings = as.numeric(rownames(subset(x$d,PointNo%in%microtubule.endings.pointno)))
    p = unique(unlist(igraph::shortest_paths(igraph::as.directed(as.ngraph(x)), from = root, to = microtubule.endings)))
    x$d$microtubules = FALSE
    x$d[p,]$microtubules = TRUE
  }
  x
}

#' @export
#' @rdname microtubules
prune_microtubules <- function(x, microtubules = TRUE){
  if(is.null(x$d$microtubules)){
    x = mark.microtubules(x)
  }
  mt = as.numeric(rownames(subset(x$d,microtubules==TRUE)))
  nat::prune_vertices(x, verticestoprune = mt, invert = microtubules)
}

#' @export
#' @rdname microtubules
visualise.microtubules <-function(x, soma = TRUE, WithConnectors = FALSE,...){
  mt = prune_microtubules(x,microtubules = TRUE)
  twigs = prune_microtubules(x,microtubules = FALSE)
  rgl::plot3d(mt, col = "darkred", WithNodes = FALSE, soma = soma, WithConnectors = WithConnectors)
  rgl::plot3d(twigs, col = "chartreuse4", WithNodes = FALSE, WithConnectors = WithConnectors)
}





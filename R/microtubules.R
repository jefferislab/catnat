# Functions for working with microtubules

#' Functions to assign and visualise microtubule rich and twig portions of a
#' neuron
#'
#' @description Manually assign the dendrite and axon to neurons / a neuron
#'
#' @param x a neuron/neuronlist object
#' @param microtubules whether to return the microtubule containing arbour
#'   (TRUE) or twigs (FALSE)
#' @param skid skeleton ID of CATMAID neuron for checking whether there are
#'   presynapses marked as being on a microtubule-lacking twig
#' @param ... Additional arguments passed to nlapply
#' @return The neuron/neuronlist object with axon/dendrite info assigned in SWC
#'   format to neuron$d
#' @export
#' @rdname microtubules
mark.microtubules <-function(x, ...) UseMethod("mark.microtubules")

#' @export
#' @rdname microtubules
mark.microtubules.neuron <- function(x, ...){
  if(is.null(x$d$microtubules)){
    if(is.null(x$tags$`microtubules end`)){
      x$d$microtubules = NA
      x$connectors$microtubules = NA
      warning("No microtubular endings marked in CATMAID neuron")
    }else{
      root = nat::rootpoints(x)
      leaves = nat::endpoints(x)
      microtubule.endings.pointno = x$tags$`microtubules end`
      microtubule.endings = as.numeric(rownames(subset(x$d,PointNo%in%microtubule.endings.pointno)))
      splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))
      i = igraph::shortest_paths(igraph::as.directed(nat::as.ngraph(x)), from = root, to = leaves, mode = "out")$vpath
      p = c()
      for(ii in 1:length(i)){
        pos = which(i[ii][[1]]%in%microtubule.endings)
        iii = splitAt(i[ii][[1]], pos)
        p = c(p,unlist(iii[1]))
      }
      p = unique(p)
      x$d$microtubules = FALSE
      x$d[p,]$microtubules = TRUE
      relevant.points = subset(x$d, PointNo%in%x$connectors$treenode_id)
      x$connectors$microtubules = relevant.points[match(x$connectors$treenode_id,relevant.points$PointNo),]$microtubules
    }
  }
  class(x) = c("catmaidneuron","neuron")
  x
}

#' @export
#' @rdname microtubules
mark.microtubules.neuronlist <- function(x, ...){
  nat::nlapply(x, mark.microtubules.neuron, ...)
}

#' @export
#' @rdname microtubules
prune_microtubules <-function(x, ...) UseMethod("prune_microtubules")

#' @export
#' @rdname microtubules
prune_microtubules.neuron <- function(x, microtubules = TRUE, ...){
  if(is.null(x$d$microtubules)){
    x = mark.microtubules.neuron(x)
  }
  mt = as.numeric(rownames(subset(x$d,microtubules==TRUE)))
  x = prune_vertices.catmaidneuron(x, verticestoprune = mt, invert = microtubules)
  x
}

#' @export
#' @rdname microtubules
prune_microtubules.neuronlist <- function(x, microtubules = TRUE, ...){
  nat::nlapply(x,prune_microtubules.neuron, microtubules = microtubules, ...)
}

#' @export
#' @rdname microtubules
#' @param soma for \code{visualise.microtubules} whether to show the soma
#' @param WithConnectors for \code{visualise.microtubules} whether to show
#'   connectors
visualise.microtubules <-function(x, soma = TRUE, WithConnectors = FALSE,...){
  mt = prune_microtubules(x,microtubules = TRUE)
  twigs = prune_microtubules(x,microtubules = FALSE)
  rgl::plot3d(mt, col = "darkred", WithNodes = FALSE, soma = soma, WithConnectors = WithConnectors)
  rgl::plot3d(twigs, col = "chartreuse4", WithNodes = FALSE, WithConnectors = WithConnectors)
}

#' @export
#' @rdname microtubules
microtubules.errors<-function(skid){
  x = read.neuron.catmaid(skid)
  x = mark.microtubules(x)
  df = subset(x$connectors,prepost==0&microtubules==FALSE)
  if(nrow(df)>0){
    message("Presynapses have been marked on microtubule-lacking twigs")
    catmaid_urls(df)
  }else{
    message("There are no presynapses marked without microtubule")
  }
}


#' Assign Strahler stream order to neurites
#'
#' @description Assign Strahler stream order to to neurons / a neuron
#'
#' @param x a neuron/neuronlist object
#' @param ... Additional arguments passed to nlapply
#' @export
#' @rdname assign_strahler
assign_strahler <-function(x, ...) UseMethod("assign_strahler")

#' @export
#' @rdname assign_strahler
assign_strahler.neuron<-function(x, ...){
  if(ifelse(!is.null(x$nTrees),x$nTrees!=1,FALSE)){
    warning("Neuron has multiple trees, calculating Strahler order for each subtree separately")
    x$d$strahler_order = 1
    for(tree in 1:x$nTrees){
      v = unique(unlist(x$SubTrees[tree]))
      if(length(v)<2){
        x$d[x$d$PointNo%in%v,]$strahler_order = 1
      }else{
        neuron = tryCatch(nat::prune_vertices(x,verticestoprune = v, invert = TRUE), error=function(e) NULL)
        if(sum(branchpoints(x)%in%v)==0){
          x$d[x$d$PointNo%in%v,]$strahler_order = 1
        }else if (!is.null(neuron)){
          s = nat::strahler_order(neuron)
          x$d[x$d$PointNo%in%v,]$strahler_order = s$points
        }
      }
    }
  }else{
    s = nat::strahler_order(x)
    x$d$strahler_order = s$points
  }
  if("catmaidneuron"%in%class(x)){
    relevant.points = subset(x$d, PointNo%in%x$connectors$treenode_id)
    x$connectors$strahler_order = relevant.points[match(x$connectors$treenode_id,relevant.points$PointNo),]$strahler_order
  }
  x
}

#' @export
#' @rdname assign_strahler
assign_strahler.neuronlist<-function(x, ...){
  nlapply(x, assign_strahler.neuron, ...)
}


#' Calculate geodesic distance from nodes to a neuron's axon-dendrite
#' branchpoint
#'
#' @description alculate geodesic distance from nodes to a neuron's primary,
#'   axon-dendrite branchpoint
#'
#' @param x a neuron/neuronlist object that has primary neurites marked (Label =
#'   7) and soma as the root
#' @param graph.distance whether to calculate the graph distance (defualt)
#'   between nodes and the primary branchpoint, or the cable length
#' @param ... Additional arguments passed to nlapply
#' @export
#' @rdname distance.from.first.branchpoint
distance.from.first.branchpoint <-function(x, ...) UseMethod("distance.from.first.branchpoint")

#' @export
#' @rdname distance.from.first.branchpoint
distance.from.first.branchpoint.neuron<-function(x, graph.distance = TRUE, ...){
  # Find axon-dendrite branch point
  pn = subset(x$d,Label==7)
  bp = nat::endpoints(pn)[!nat::endpoints(pn)%in%nat::rootpoints(pn)]
  bp = as.numeric(rownames(subset(x$d,PointNo==bp)))
  n=nat::as.ngraph(x)
  path = suppressWarnings(igraph::shortest_paths(n,from=bp, mode= "out")$vpath)
  x$d$geodesic.distance = sapply(path,length)
  if(!graph.distance){
    conns = as.numeric(rownames((subset(x$d,PointNo%in%x$connectors$treenode_id))))
    paths = suppressWarnings(igraph::shortest_paths(n,from=bp, to = conns, mode= "out")$vpath)
    real.lengths = c()
    for(p in paths){
      if(length(p)>2){
        rl = summary(nat::prune_vertices(x,verticestoprune = p,invert=TRUE))$cable.length
        real.lengths = c(real.lengths,rl)
      }else{
        real.lengths = c(real.lengths,0)
      }
    }
    x$d$geodesic.distance = NA
    x$d$geodesic.distance[conns] = real.lengths
  }
  relevant.points = subset(x$d, PointNo%in%x$connectors$treenode_id)
  x$connectors$geodesic.distance = relevant.points[match(x$connectors$treenode_id,relevant.points$PointNo),]$geodesic.distance
  x
  #real.distance = summary(nat::prune_vertices(x,verticestoprune = unlist(path),invert=TRUE))
}
#' @export
#' @rdname distance.from.first.branchpoint
distance.from.first.branchpoint.neuronlist<-function(x, graph.distance = TRUE, ...){
  nlapply(x,distance.from.first.branchpoint.neuron, graph.distance = graph.distance, ...)
}



#' Get a subtree in a neuron object as a separate neuron object
#'
#' @description Get a subtree in a neuron object as a separate neuron object
#'
#' @param neuron a neuron object from which to extract a subtree
#' @param subtree the number of the subtree to extract
#' @export
#' @rdname subtree
subtree<-function(neuron, subtree = 1){
  if(neuron$nTrees<1){
    warning("Neuron has no subtree")
  }else{
    v = unique(unlist(neuron$SubTrees[subtree]))
    if(length(v)>1){
      neuron = nat::prune_vertices(neuron, verticestoprune = v, invert = TRUE)
    }else{
      neuron = NULL
    }
  }
  neuron
}




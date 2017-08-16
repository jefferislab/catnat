# Functions for working with microtubules

#' Functions to assign and visualise microtubule rich and twig portions of a neuron
#'
#' @description Manually assign the dendrite and axon to neurons / a neuron
#'
#' @param x a neuron/neuronlist object
#' @param microtubules whether to return the microtubule containing arbour (TRUE) or twigs (FALSE)
#' @return The neuron/neuronlist object with axon/dendrite info assigned in SWC format to neuron$d
#' @export
#' @rdname microtubules
mark.microtubules <-function(x, ...) UseMethod("mark.microtubules")

#' @export
#' @rdname microtubules
mark.microtubules.neuron <- function(x){
  if(is.null(x$d$microtubules)){
    if(is.null(x$tags$`microtubules end`)){
      message("No microtubular endings marked in CATMAID neuron")
      break
    }
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
    #p = c()
    #for(m in microtubule.endings){
    #  p = c(p,unique(unlist(suppressWarnings(igraph::shortest_paths(igraph::as.directed(as.ngraph(x)), from = microtubule.endings[1], to = leaves, mode = "out")))))
    #}
    #p = unique(p)
    x$d$microtubules = FALSE
    x$d[p,]$microtubules = TRUE
    relevant.points = subset(x$d, PointNo%in%x$connectors$treenode_id)
    x$connectors$microtubules = relevant.points[match(x$connectors$treenode_id,relevant.points$PointNo),]$microtubules
  }
  #class(x) = c("catmaidneuron","list")
  x
}

#' @export
#' @rdname microtubules
mark.microtubules.neuronlist <- function(x){
  nat::nlapply(x, mark.microtubules.neuron)
}

#' @export
#' @rdname microtubules
prune_microtubules <-function(x, ...) UseMethod("prune_microtubules")

#' @export
#' @rdname microtubules
prune_microtubules.neuron <- function(x, microtubules = TRUE){
  if(is.null(x$d$microtubules)){
    x = mark.microtubules.neuron(x)
  }
  mt = as.numeric(rownames(subset(x$d,microtubules==TRUE)))
  x = nat::prune_vertices(x, verticestoprune = mt, invert = microtubules)
  #class(x) = c("catmaidneuron","list")
  x
}

#' @export
#' @rdname microtubules
prune_microtubules.neuronlist <- function(x, microtubules = TRUE){
  nat::nlapply(x,prune_microtubules.neuron, microtubules = microtubules)
}

#' @export
#' @rdname microtubules
visualise.microtubules <-function(x, soma = TRUE, WithConnectors = FALSE,...){
  mt = prune_microtubules(x,microtubules = TRUE)
  twigs = prune_microtubules(x,microtubules = FALSE)
  rgl::plot3d(mt, col = "darkred", WithNodes = FALSE, soma = soma, WithConnectors = WithConnectors)
  rgl::plot3d(twigs, col = "chartreuse4", WithNodes = FALSE, WithConnectors = WithConnectors)
}

# function needed to find presynapses marked on non microtubular fragments and so correct
microtubules.errors<-function(x){
  message("Neuron must be in FAFB13 space!")
  df = subset(x$connectors,prepost==0&microtubules==FALSE)
  catmaid_urls(df)
}

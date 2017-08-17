# Functions for working with microtubules

#' Functions to assign and visualise microtubule rich and twig portions of a neuron
#'
#' @description Manually assign the dendrite and axon to neurons / a neuron
#'
#' @param x a neuron/neuronlist object
#' @param microtubules whether to return the microtubule containing arbour (TRUE) or twigs (FALSE)
#' @param skid skeleton ID of CATMAID neuron for checking whether there are presynapses marked as being on a microtubule-lacking twig
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
#' @export
#' @rdname assign_strahler
assign_strahler <-function(x, ...) UseMethod("assign_strahler")

#' @export
#' @rdname assign_strahler
assign_strahler.neuron<-function(x){
  s = strahler_order(x)
  x$d$strahler_order = s$points
  if("catmaidneuron"%in%class(x)){
    relevant.points = subset(x$d, PointNo%in%x$connectors$treenode_id)
    x$connectors$strahler_order = relevant.points[match(x$connectors$treenode_id,relevant.points$PointNo),]$strahler_order
  }
  x
}

#' @export
#' @rdname assign_strahler
assign_strahler.neuronlist<-function(x){
  nlapply(x, assign_strahler.neuron)
}


#' Calculate geodesic distance from nodes to a neuron's axon-dendrite branchpoint
#'
#' @description alculate geodesic distance from nodes to a neuron's primary, axon-dendrite branchpoint
#'
#' @param x a neuron/neuronlist object that has primary neurites marked (Label = 7) and soma as the root
#' @param graph.distance whether to calculate the graph distance (defualt) between nodes and the primary branchpoint, or the cable length
#' @export
#' @rdname distance.from.first.branchpoint
distance.from.first.branchpoint <-function(x, ...) UseMethod("distance.from.first.branchpoint")

#' @export
#' @rdname distance.from.first.branchpoint
distance.from.first.branchpoint.neuron<-function(x, graph.distance = TRUE){
  # Find axon-dendrite branch point
  pn = subset(x$d,Label==7)
  bp = endpoints(pn)[!endpoints(pn)%in%rootpoints(pn)]
  bp = as.numeric(rownames(subset(x$d,PointNo==bp)))
  n=as.ngraph(x)
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
distance.from.first.branchpoint.neuronlist<-function(x, graph.distance = TRUE){
  nlapply(x,distance.from.first.branchpoint.neuron, graph.distance = graph.distance)
}

# Hidden
catmaid_urls <- function (df) {
  if (!is.data.frame(df)){
    stop("Please give me a data frame!")
  }
  base = "https://neuropil.janelia.org/tracing/fafb/v13"
  catmaid_url = paste0(base, "?pid=1")
  catmaid_url = paste0(catmaid_url, "&zp=", df[["z"]])
  catmaid_url = paste0(catmaid_url, "&yp=", df[["y"]])
  catmaid_url = paste0(catmaid_url, "&xp=", df[["x"]])
  catmaid_url = paste0(catmaid_url, "&tool=tracingtool")
  id = if (is.null(df$partner_skid)) {
    df$connector_id
  }else {
    df$partner_skid
  }
  catmaid_url = paste0(catmaid_url, "&active_skeleton_id=",
                       id)
  catmaid_url = paste0(catmaid_url, "&sid0=5&s0=0")
  catmaid_url
}

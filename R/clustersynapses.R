#' Cluster synapses within a neuron's skeleton
#'
#' @description implementation of the algorithm for clustering synapses from Schneider-Mizell et al. (2016). Note that the abrbour.cluster() function will retrieve these clusters as separate neuron objects in a neuronlist.
#'
#' @param x a neuronlist or neuron object
#' @param polyadic Whether to count presynapses as a single synapse or as the number of connections that prsynapse makes for the purpose of clustering. Defaults to true.
#' @param lambda A bandwidth parameter that effectively determines the size of clusters
#' @param order Helps determine cluster size. How many nodes to consider as being inside the neighbourhood at each step during gradient ascent.
#' @param e The entropy value calculated between post and pre synapses in a cluster abvoe which we assign that cluster as part of a dendrite, and below which we assign as an axonic segment. In between these values we assign the cluster as mixed.
#' @param ... additional arguments passed to methods.
#'
#' @details From Schneider-Mizell et al. (2016): "This approach involves convolving synapse locations with a Gaussian kernel to estimate the density of synapses in space. A cluster is then the set of synapses for which, starting at their location, gradient ascent reaches the same density peak. However, loca- tions on one neuron that are close in space can be very far apart along the neuron. Here, instead of considering the density of a neuron’s synapses in 3d space, we use a similar procedure to estimate the density of synapses at every point on the arbor (following the cable) and define synapse clusters in the same manner. The only parameter in both approaches is the width of the Gaussian kernel, a physically meaningful parameter."
#'
#' @references Schneider-Mizell, C. M., Gerhard, S., Longair, M., Kazimiers, T.,
#'   Li, F., Zwart, M. F., … Cardona, A. (2015). Quantitative neuroanatomy for
#'   connectomics in Drosophila. bioRxiv, 026617. http://doi.org/10.1101/026617
#'
#' @return the neuron or neuron list object inputted, with centipetal flow
#'   centrality information added to neuron$d, a segregation idnex score and
#'   estimation of neuronal type (intertneuron or PN) based on this score (>0.05
#'   = PN).
#' @export
#' @rdname cluster_synapses_within_skeleton
#' @seealso \code{\link{seebroken3d}} \code{\link{flow.centrality}}
cluster_synapses_within_skeleton <-function(neuron, polyadic = T, lambda = 30, order = 150, e = c(0.3,0.7),...) UseMethod("cluster_synapses_within_skeleton")

#' @export
#' @rdname cluster_synapses_within_skeleton
cluster_synapses_within_skeleton.neuron <- function(neuron, polyadic = T, lambda = 30, order = 150, e = c(0.3,0.7),...){
  require(igraph)
  el = neuron$d[neuron$d$Parent != -1, c("Parent", "PointNo")] # Get list of soma=leaf directed conenctions
  n = nat::ngraph(data.matrix(el[,2:1]), neuron$d$PointNo, directed = TRUE, xyz = nat::xyzmatrix(neuron$d),
                  diam = neuron$d$W) # Make ngraph object, but for centripetal, invert the el list
  # Get comprehensive paths list
  leaves = which(igraph::degree(n, v = igraph::V(n), mode = "in")==0, useNames = T)
  root= which(igraph::degree(n, v = igraph::V(n), mode = "out")==0, useNames = T)
  segs = neuron$SegList # Get the segment list
  nodes = neuron$d # get neuron's node data
  nodes[,"post"] <- 0 # raw synapse number at this location
  nodes[,"pre"] <- 0 # raw synapse number at this location
  nodes = nodes[unlist(c(root, lapply(segs, function (x) x[-1]))),]
  syns.in = neuron$connectors[neuron$connectors[,3]==1,][,1]
  if (polyadic == T){
    pres = neuron$connectors[neuron$connectors[,3]==0,][,2]
    pre.cons = catmaid_get_connectors(pres)$connector_id
    syns.out = neuron$connectors[,1][match(pre.cons, neuron$connectors[,2])]
  }else{
    syns.out = neuron$connectors[neuron$connectors[,3]==0,][,1]
  }
  # Rearrange so nodes are the indices and we can count no. of synapses to which nodes connect
  point.no.in = rownames(nodes)[match(syns.in,nodes[,"PointNo"])]
  nodes.in = rep(1,length(point.no.in))
  names(nodes.in) = point.no.in
  nodes.in = tapply(nodes.in, point.no.in, sum)
  point.no.out = rownames(nodes)[match(syns.out,nodes[,"PointNo"])]
  nodes.out = rep(1,length(point.no.out))
  names(nodes.out) = point.no.out
  nodes.out = tapply(nodes.out, point.no.out, sum)
  # Add more accurate synapose positions
  nodes[names(nodes.in),"post"] <- nodes.in
  nodes[names(nodes.out),"pre"] <- nodes.out
  # Claculate distance matrices
  dis.matrix.in = igraph::distances(as.undirected(n))[,as.integer(point.no.in)]
  rownames(dis.matrix.in) = rownames(nodes)
  colnames(dis.matrix.in) = point.no.in
  dis.matrix.out = igraph::distances(as.undirected(n))[,as.integer(point.no.out)]
  rownames(dis.matrix.out) = rownames(nodes)
  colnames(dis.matrix.out) = point.no.out
  synapse.density.in = lapply(rownames(nodes), function(x) sum(unlist(lapply(dis.matrix.in[x,], function(y) exp(-y^2/(2*lambda^2))))))
  synapse.density.out = lapply(rownames(nodes), function(x) sum(unlist(lapply(dis.matrix.out[x,], function(y) exp(-y^2/(2*lambda^2))))))
  nodes[,"post.density.score"] <- unlist(synapse.density.in)
  nodes[,"pre.density.score"] <- unlist(synapse.density.out)
  # Calculate for both signs of synapse combined
  dis.matrix = distances(as.undirected(n))[,as.integer(c(point.no.in,point.no.out))]
  rownames(dis.matrix) = rownames(nodes)
  colnames(dis.matrix) = c(point.no.in,point.no.out)
  synapse.density = lapply(rownames(nodes), function(x) sum(unlist(lapply(dis.matrix[x,], function(y) exp(-y^2/(2*lambda^2))))))
  nodes[,"syn.density.score"] <- unlist(synapse.density)
  #synapse.density.in = lapply(rownames(nodes), function(x) sum(unlist(lapply(point.no.in, function(y) exp(-distances(as.undirected(n),v=x,to=y)^2/(2*lambda^2))))))
  #synapse.density.out = lapply(rownames(nodes), function(x) lapply(point.no.in, function(y) exp(-distances(n,v=x,to=y)/(2*lambda^2))))
  ## Gradient ascent
  gradient.ascent <- function(x, ngraph, nodes, clusters, order){
    peak = F
    positions = c()
    while (peak == F){
      neighbours = unlist(igraph::neighborhood(igraph::as.undirected(ngraph), order = order, nodes = x))
      scores = nodes[as.character(neighbours),"syn.density.score"]
      if (match(max(scores),scores) != 1){
        positions = c(positions,unlist(igraph::shortest_paths(graph = as.undirected(n), from = x, to = neighbours[match(max(scores),scores)])))
        positions = positions[!positions%in%unlist(clusters)]
        x = neighbours[match(max(scores),scores)]
        if (x%in%unlist(clusters)){
          clusters = lapply(seq_along(clusters), function(c) if(x%in%clusters[[c]]){c(positions,clusters[[c]])}else{clusters[[c]]})
          peak = T
        }
      }else{
        positions = c(positions, x)
        if(!is.null(positions)){clusters[[length(clusters)+1]] <- unique(positions)}
        peak = T
      }
    }
    return(clusters) # the last position in the list will be the peak
  }
  # Set the scores of leaves to zero, so that they must ascend
  # nodes[as.character(leaves),"syn.density.score"] <- 0
  # To make this run faster, we can start from the lesser scored nodes
  density.scores = nodes[,"syn.density.score"]
  names(density.scores) = rownames(nodes)
  descending = sort(density.scores, decreasing = F)
  # Start the algorithm
  clusters = list()
  for (node in names(descending)){
    if(!node%in%unlist(clusters)){
      clusters = gradient.ascent(x = node, ngraph = n, nodes =nodes, clusters = clusters, order = order)
    }
  }
  sis = c()
  ns = c()
  for (c in 1:length(clusters)){
    cluster = clusters[[c]]
    # Claculate the entropy
    n.in = sum(point.no.in%in%cluster)
    n.out = sum(point.no.out%in%cluster)
    n.pi = n.in/(n.in+n.out)
    if(is.nan(n.pi)){n.pi = 0}
    n.si = -(n.pi*log(n.pi)+(1-n.pi)*log(1-n.pi))
    if(is.nan(n.si)){n.si = 0}
    sis = c(sis, n.si)
    ns =  c(ns,(n.in+n.out))
    nodes[as.character(cluster),"cluster.entropy"] = n.si
    if(n.pi > e[2]){type = "dendrite"
    }else if(n.pi < e[1]){ type = "axon"}
    if (n.pi < e[2] & n.pi > e[1]){type = "mixed"}
    nodes[as.character(cluster),"cluster"] = paste(type, c,sep = '-')
  }
  entropy.score = (1/(sum(ns)))*sum(sis*ns)
  both.comps = (sum(point.no.in%in%unlist(clusters)))/(sum(ns))
  control.score = -(both.comps*log(both.comps)+(1-both.comps)*log(1-both.comps))
  segregation.index = 1 - (entropy.score/control.score)
  if(is.na(segregation.index)) { segregation.index = 0 }
  neuron$cluster.segregation.index = segregation.index
  nodes = nodes[order(as.numeric(rownames(nodes))),]
  neuron$d = nodes
  neuron
}

#' @export
#' @rdname cluster_synapses_within_skeleton
cluster_synapses_within_skeleton.neuronlist <- function(someneuronlist, polyadic = T, lambda = 30, order = 150, e = c(0.3,0.7),...) {
  neurons = nat::nlapply(someneuronlist, cluster_synapses_within_skeleton, polyadic = polyadic, lambda = lambda, order = order, e = e, OmitFailures = T)
  neurons
}

#' Plot neurons split up by synapse clusters
#'
#' @param someneuronlist a neuronlist or neuron object that has been modified by flow.centrality
#' @param col colours of sections. Defaults to orange or axons, green for primary dendrite, blue for dendrites and pink for nodes with no flow.
#' @param WithConnectors whether ot plot the anatomical location of pre (red) and post (cyan) synapses.
#' @param soma whether to plot a soma, and what the radius should be
#' @param WithNodes whether to plot branch points
#' @param ... additional arguments passed to methods.
#'
#' @return Plots neuron(s) with arbours coloured by synapse cluster
#' @export
#' @rdname seebroken3d
#' @seealso \code{\link{cluster_synapses_within_skeleton}} \code{\link{seesplit3d}}
seebroken3d = function(neuron, WithConnectors = T, WithNodes = F, soma = 100){
  if(is.null(neuron$cluster.segregation.index)){
    warning("No synapse clustering calculated, dropping neuron")
    break
  }
  clusters = unique(neuron$d$cluster)
  count = length(clusters)
  col = rainbow(count)
  #dend.col = colorRampPalette(colors = c("skyblue", "darkblue"))(count)
  #mixed.col = colorRampPalette(colors = c("purple", "violetred"))(count)
  #axon.col = colorRampPalette(colors = c("orange", "red4"))(count)
  for (c in 1:length(clusters)){
    cluster = clusters[[c]]
    soma = F
    cluster.v = subset(rownames(neuron$d), neuron$d$cluster == cluster)
    cluster.g = nat::prune_vertices(neuron, verticestoprune = rownames(neuron$d)[!rownames(neuron$d)%in%cluster.v])
    #p.in = 100*cluster.g$d$pre/(cluster.g$d$pre+cluster.g$d$post)# colour by proportion of input synapses
    # col = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100)
    #rgl::plot3d(cluster.g, col = col[p.in], WithNodes = WithNodes, soma = soma)
    #if(grepl("dendrite",cluster)){col = dend.col} # Colour by compartment assignation
    #if(grepl("mixed",cluster)){col = mixed.col}
    #if(grepl("axon",cluster)){col = axon.col}
    if (neuron$StartPoint%in%cluster.v){soma = soma.size}
    rgl::plot3d(cluster.g, col = col[c], WithNodes = WithNodes, soma = soma)
  }
  if (WithConnectors == T){
    rgl::points3d(subset(xyzmatrix(neuron$d),neuron$d$post>0), col = 'cyan')
    rgl::points3d(subset(xyzmatrix(neuron$d),neuron$d$pre>0), col = 'red')
  }
}

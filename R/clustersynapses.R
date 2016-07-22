cluster.synapses <- function(neuron, polyadic = T, lambda = 30, order = 150){
  #neuron = read.neuron.catmaid("13581830")/1000
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
  dis.matrix.in = distances(as.undirected(n))[,as.integer(point.no.in)]
  rownames(dis.matrix.in) = rownames(nodes)
  colnames(dis.matrix.in) = point.no.in
  dis.matrix.out = distances(as.undirected(n))[,as.integer(point.no.out)]
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
      neighbours = unlist(neighborhood(as.undirected(ngraph), order = order, nodes = x))
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
    if(n.pi > 0.7){type = "dendrite"
    }else if(n.pi < 0.3){ type = "axon"}
    if (n.pi < 0.7 & n.pi > 0.3){type = "local"}
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




#' Plot neurons split up by synapse clusters
#'
#' @param someneuronlist a neuronlist or neuron object that has been modified by flow.centrality
#' @param col colours of sections. Defaults to orange or axons, green for primary dendrite, blue for dendrites and pink for nodes with no flow.
#' @param primary.dendrite Type of object for Deformetrica deformation. See Deformetrica's documentation. Default is appropriate for neuron skeletons.
#' @param WithConnectors whether ot plot the anatomical location of pre (red) and post (cyan) synapses.
#' @param soma whether to plot a soma, and what the radius should be
#' @param WithNodes whether to plot branch points
#' @param highflow wheather to plot the nodes of highest (with in one standard deviation less than maximum) flow centrality (pink points)
#' @param ... additional arguments passed to methods.
#'
#' @return Plots coloured neuron(s)
#' @export
#' @rdname seesplit3d
#' @seealso \code{\link{flow.centrality}} \code{\link{get.synapses}}
seebroken3d = function(neuron, WithConnectors = T, WithNodes = F, soma.size = 100){
  if(is.null(neuron$cluster.segregation.index)){
    warning("No synapse clustering calculated, dropping neuron")
    break
  }
  clusters = unique(neuron$d$cluster)
  count = length(clusters)
  col = rainbow(count)
  #dend.col = colorRampPalette(colors = c("skyblue", "darkblue"))(count)
  #local.col = colorRampPalette(colors = c("purple", "violetred"))(count)
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
    #if(grepl("local",cluster)){col = local.col}
    #if(grepl("axon",cluster)){col = axon.col}
    if (neuron$StartPoint%in%cluster.v){soma = soma.size}
    rgl::plot3d(cluster.g, col = col[c], WithNodes = WithNodes, soma = soma)
  }
  if (WithConnectors == T){
    rgl::points3d(subset(xyzmatrix(neuron$d),neuron$d$post>0), col = 'cyan')
    rgl::points3d(subset(xyzmatrix(neuron$d),neuron$d$pre>0), col = 'red')
  }
}



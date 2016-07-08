#' Cluster neurons by pre- and postsynapse positons
#'
#' @description implementation of the algorithm for clustering neurons by synapse location from Schlegel et al. (2016). Assumes neurons are scaled to microns.
#'
#' @param someneuronlist a neuronlist or neuron object
#' @param sigma determines what distances between two synapses are considered close (defaults to 2 um)
#' @param omega synapse cluster radius. Defaults to sigma.
#' @param symmetric whether to return a symmetric martrix (average of scores between two neurons in both directions)
#' @param ... additional arguments passed to methods.
#'
#' @details From Schneider-Mizell et al. (2016): "We use flow centrality for four purposes. First, to split an arbor into axon and dendrite at the maximum centrifugal SFC, which is a preliminary step for computing the segregation index, for expressing all kinds of connectivity edges (e.g. axo-axonic, dendro-dendritic) in the wiring diagram, or for rendering the arbor in 3d with differently colored regions. Second, to quantitatively estimate the cable distance between the axon terminals and dendritic arbor by measuring the amount of cable with the maximum centrifugal SFC value. Third, to measure the cable length of the main den- dritic shafts using centripetal SFC, which applies only to insect neurons with at least one output syn- apse in their dendritic arbor. And fourth, to weigh the color of each skeleton node in a 3d view, providing a characteristic signature of the arbor that enables subjective evaluation of its identity."
#'
#' @return A matrix of similarity scores between inputted neurons, based on synapse positions.
#' @export
#' @rdname cluster.by.synapses
#' @seealso \code{\link{plot3d.split}} \code{\link{get.synapses}}
cluster.by.synapses <- function(someneuronlist, sigma = 2, omega = sigma, symmetric = T, ...){
  m = matrix(nrow = length(someneuronlist), ncol = length(someneuronlist))
  colnames(m) = rownames(m) = names(someneuronlist)
  for (neuron in 1:length(someneuronlist)){
    g = get.synapses(someneuronlist[neuron], "BOTH")
    for (neuron2 in 1:length(someneuronlist)){
      t = get.synapses(someneuronlist[neuron2], "BOTH")
      scores = c()
      for (syn in 1:nrow(g)){
        gg = subset(g, g$prepost == g[syn,"prepost"])[,-4]
        tt = subset(t, t$prepost == g[syn,"prepost"])[,-4]
        if(empty(tt)){ score = 0; break}
        n = nabor::knn(tt, nat::xyzmatrix(g[syn,]), k =1)
        close = t[n$nn.idx,]
        gn = nabor::knn(nat::xyzmatrix(g[syn,]), gg, k =1)
        gn = sum(gn$nn.dists<omega) - 1
        tn = nabor::knn(nat::xyzmatrix(close), tt, k =1)
        tn = sum(tn$nn.dists<omega) - 1
        multiplier = exp(abs(gn-tn)/(tn+gn))
        if(is.infinite(multiplier)|is.na(multiplier)|is.nan(multiplier)){ multiplier = 0}
        score = exp((-n$nn.dist^2)/2*(sigma)^2)*multiplier
        scores = c(scores,score)
      }
      m[neuron,neuron2] = mean(scores)
    }
  }
  if (symmetric == T){
    pmean <- function(x,y) (x+y)/2
    m[] <- pmean(m, matrix(m, nrow(s), byrow=TRUE))
  }
  return(m)
}









flow.centrality.neuron <- function(neuron, mode = c("average","centrifugal","centripetal"), polypre = T, primary.dendrite = 0.85, ...){
  # prune Strahler first...and use segmentgraph?
  # Generate ngraph object
  el = neuron$d[neuron$d$Parent != -1, c("Parent", "PointNo")] # Get list of soma=leaf directed conenctions
  n = nat::ngraph(data.matrix(el[,2:1]), neuron$d$PointNo, directed = TRUE, xyz = nat::xyzmatrix(neuron$d),
                  diam = neuron$d$W) # Make ngraph object, but for centripetal, invert the el list
  # Get comprehensive paths list
  leaves = which(igraph::degree(n, v = V(n), mode = "in")==0, useNames = T)
  root= which(igraph::degree(n, v = V(n), mode = "out")==0, useNames = T)
  segs = neuron$SegList # Get the segment list
  nodes = neuron$d # get neuron's node data
  nodes[,"post"] <- 0 # raw synapse number at this location
  nodes[,"pre"] <- 0 # raw synapse number at this location
  nodes[,"up.syns.in"] <- 0 # Here, we are going to culmulatively count synapses
  nodes[,"up.syns.out"] <- 0 # Same but for outputs
  nodes[,"flow.cent"] <- 0 # We'll addd the score for each node here
  nodes[,"compartment"] <- 'dendrite'
  nodes = nodes[unlist(c(root, lapply(segs, function (x) x[-1]))),]
  syns.in = neuron$connectors[neuron$connectors[,3]==1,][,1]
  if (polypre == T){
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
  # Add running totals from the seglist and synapse positions to data matrix
  ins = c(0,lapply(segs, function(x) if(!is.null(x)){rev(cumsum(rev(unlist(lapply(x, function(x) ifelse(x%in%names(nodes.in),nodes[as.character(x),"post"],0))))))}))
  outs = c(0,lapply(segs, function(x) if(!is.null(x)){rev(cumsum(rev(unlist(lapply(x, function(x) ifelse(x%in%names(nodes.out),nodes[as.character(x),"pre"],0))))))}))
  nodes[,"up.syns.in"] = nodes[,"up.syns.in"] + c(0,unlist(lapply(ins, function(x) x[-1])))
  nodes[,"up.syns.out"] = nodes[,"up.syns.out"] + c(0,unlist(lapply(outs, function(x) x[-1])))
  # Add in the branch node point carried on from other upstream branches
  in.bps = unlist(lapply(ins, function(x) x[1]))[-1]
  out.bps = unlist(lapply(outs, function(x) x[1]))[-1]
  names(in.bps) = names(out.bps) = bps = c(root,unlist(lapply(segs, function(x) x[1]))[-1])
  in.bps[names(in.bps)%in%names(nodes.in)] = in.bps[names(in.bps)%in%names(nodes.in)] - nodes[names(in.bps)[names(in.bps)%in%names(nodes.in)],"post"]
  out.bps[names(out.bps)%in%names(nodes.out)] = out.bps[names(out.bps)%in%names(nodes.out)] - nodes[names(out.bps)[names(out.bps)%in%names(nodes.out)],"pre"]
  # Propagate scores from branch nodes
  for (i in 1:length(bps)){
    bp = bps[i]
    new.in = new.out = c(rep(0,nrow(nodes)))
    names(new.in) = names(new.out) = rownames(nodes)
    vertices = unlist(igraph::shortest_paths(n, bp, to = root)$vpath)[-1]
    new.in[as.character(vertices)] = new.in[as.character(vertices)] + in.bps[i]
    new.out[as.character(vertices)] = new.out[as.character(vertices)] + out.bps[i]
    nodes[,"up.syns.in"] = nodes[,"up.syns.in"]+ new.in
    nodes[,"up.syns.out"] = nodes[,"up.syns.out"] + new.out
  }
  # Calculate flow centrality
  in.total = nodes[1,"up.syns.in"] = length(point.no.in)
  out.total = nodes[1,"up.syns.out"] = length(point.no.out)
  if(mode[1] == "centrifugal"){nodes[,"flow.cent"] = (in.total - nodes[,"up.syns.in"])*nodes[,"up.syns.out"]}
  if(mode[1] == "centripetal"){nodes[,"flow.cent"] = (out.total - nodes[,"up.syns.out"])*nodes[,"up.syns.in"]}
  if(mode[1] == "average"){nodes[,"flow.cent"] = (in.total - nodes[,"up.syns.in"])*nodes[,"up.syns.out"] + (out.total - nodes[,"up.syns.out"])*nodes[,"up.syns.in"]}
  nodes = nodes[order(as.numeric(rownames(nodes))),]
  ais = which(apply(nodes, 1, function(x) x["flow.cent"] == max(nodes[,"flow.cent"])))
  if (length(ais)>0){
    runstosoma = unlist(lapply(ais, function(x) length(unlist(igraph::shortest_paths(n, x, to = root)$vpath))))
    ais = ais[match(min(runstosoma),runstosoma)]
  }
  downstream = suppressWarnings(unique(unlist(igraph::shortest_paths(n, ais, to = leaves, mode = "in")$vpath)))
  nodes[as.character(downstream),"compartment"] = "axon"
  # Work out the primary neurite and empty leaf compartments
  zeros = subset(rownames(nodes),nodes[,"flow.cent"]==0)
  igraph::V(n)$name = igraph::V(n)
  nn = igraph::delete_vertices(n, v = igraph::V(n)[!igraph::V(n)%in%zeros])
  p.n = suppressWarnings(unique(unlist(igraph::shortest_paths(nn, from =root))))
  nodes[p.n,"compartment"] = "primary neurite"
  nodes[zeros[!zeros%in%p.n],"compartment"] = "primary neurite"
  if(!is.null(primary.dendrite)){
    highs = subset(rownames(nodes),nodes[,"flow.cent"]>=primary.dendrite*max(nodes[,"flow.cent"]))
    nodes[as.character(highs),"compartment"] = "primary dendrite"
  }
  # Calculate segregation score
  dendrites = subset(nodes, nodes$compartment == "dendrite")
  dendrites.post = sum(subset(dendrites$post,dendrites$post>0))
  dendrites.pre = sum(subset(dendrites$pre,dendrites$pre>0))
  dendrites.both = dendrites.post + dendrites.pre
  dendrites.pi = dendrites.post/dendrites.both
  dendrites.si = -(dendrites.pi*log(dendrites.pi)+(1-dendrites.pi)*log(1-dendrites.pi))
  if(is.nan(dendrites.si)){dendrites.si = 0}
  axon = subset(nodes, nodes$compartment == "axon")
  axon.post = sum(subset(axon$post,axon$post>0))
  axon.pre = sum(subset(axon$pre,axon$pre>0))
  axon.both = axon.post + axon.pre
  axon.pi = axon.post/axon.both
  axon.si = -(axon.pi*log(axon.pi)+(1-axon.pi)*log(1-axon.pi))
  if(is.nan(axon.si)){axon.si = 0}
  entropy.score = (1/(dendrites.both+axon.both))*(axon.si*axon.both+dendrites.si*dendrites.both)
  both.comps = nrow(subset(nodes,nodes$post>0))/(dendrites.both+axon.both)
  control.score = -(both.comps*log(both.comps)+(1-both.comps)*log(1-both.comps))
  segregation.index = 1 - (entropy.score/control.score)
  # Add new data to object
  neuron$d = nodes
  neuron$segregation.index = segregation.index
  neuron$type = ifelse(segregation.index > 0.05, "interneuron", "PN")
  neuron
}

#' @export
#' @rdname flow.centrality
flow.centrality.neuronlist <- function(someneuronlist, mode = c("average","centrifugal","centripetal"), polypre = T, primary.dendrite = 0.85, ...){
  neurons = nat::nlapply(someneuronlist, flow.centrality, mode = mode, polypre = polypre, primary.dendrite = primary.dendrite)
  neurons
}

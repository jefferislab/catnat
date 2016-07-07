
library(nabor)

# Get connectors
get.connectors <-function(x, target = c("BOTH", "PRE", "POST"), polypre = T, ...) UseMethod("get.connectors")

get.connectors.neuron <- function (someneuronlist, target = c("BOTH", "PRE", "POST"), polypre = T){
    if (target%in%c("POST","BOTH")) { 
      syns.in = neuron$connectors[neuron$connectors[,3]==1,][,1]
      point.no = rownames(neuron$d)[match(syns.in,neuron$d[,"PointNo"])]
      points = nat::xyzmatrix(neuron$d[point.no,])
      points = cbind(points, prepost = 1)
    }
    if (target%in%c("PRE","BOTH")) { 
      pres = neuron$connectors[neuron$connectors[,3]==0,][,2]
      pre.cons = catmaid_get_connectors(pres)$connector_id
      syns.out = neuron$connectors[,1][match(pre.cons, neuron$connectors[,2])]
      point.no = rownames(neuron$d)[match(syns.out,neuron$d[,"PointNo"])]
      points.out = nat::xyzmatrix(neuron$d[point.no,])
      points.out = cbind(points.out, prepost = 0)
      if(target == "BOTH"){
        points = rbind(points, points.out)
      }else{ points = points.out}
    }
    points
}

get.connectors.neuronlist <- function (someneuronlist, target = c("BOTH", "PRE", "POST"), polypre = T){
  points = nat::nlapply(someneuronlist, get.connectors.neuron, target = target, polypre = polypre)
  points = do.call(rbind, points)
  points
}

# Cluster by synapse
clusterbysynapses <- function(someneuronlist, sigma = 2, omega = sigma, symmetric = T){
  m = matrix(nrow = length(someneuronlist), ncol = length(someneuronlist))
  colnames(m) = rownames(m) = names(someneuronlist)
  for (neuron in 1:length(someneuronlist)){
    g = get.connectors(someneuronlist[neuron], "BOTH")
    for (neuron2 in 1:length(someneuronlist)){
      t = get.connectors2(someneuronlist[neuron2], "BOTH")
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

# Calculate flow centrality
flow.centrality <- function(someneuronlist, mode = c("centripetal", "centrifugal", "both"), polypre = T, primary.dendrite = 0.8){
  meta = attr(someneuronlist,'df')
  nlist = neuronlist()
  for (someneuron in 1:length(someneuronlist)){
  # prune Strahler first...and use segmentgraph
    skid = names(someneuronlist)[someneuron]
    neuron = someneuronlist[[someneuron]]
    # Generate ngraph object
    el = neuron$d[neuron$d$Parent != -1, c("Parent", "PointNo")] # Get list of soma=leaf directed conenctions
    n = ngraph(data.matrix(el[,2:1]), neuron$d$PointNo, directed = TRUE, xyz = nat::xyzmatrix(neuron$d), 
           diam = neuron$d$W) # Make ngraph object, but for centripetal, invert the el list
    # Get comprehensive paths list
    leaves = which(degree(n, v = V(n), mode = "in")==0, useNames = T)
    root= which(degree(n, v = V(n), mode = "out")==0, useNames = T)
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
    if (polpre == T){
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
      vertices = unlist(shortest_paths(n, bp, to = root)$vpath)[-1]
      new.in[as.character(vertices)] = new.in[as.character(vertices)] + in.bps[i]
      new.out[as.character(vertices)] = new.out[as.character(vertices)] + out.bps[i]
      nodes[,"up.syns.in"] = nodes[,"up.syns.in"]+ new.in
      nodes[,"up.syns.out"] = nodes[,"up.syns.out"] + new.out
    }
    # Calculate flow centrality
    in.total = nodes[1,"up.syns.in"] = length(point.no.in)
    out.total = nodes[1,"up.syns.out"] = length(point.no.out)
    nodes[,"flow.cent"] = (in.total - nodes[,"up.syns.in"])*nodes[,"up.syns.out"] + (out.total - nodes[,"up.syns.out"])*nodes[,"up.syns.in"]
    nodes = nodes[order(as.numeric(rownames(nodes))),]
    ais = which(apply(nodes, 1, function(x) x["flow.cent"] == max(nodes[,"flow.cent"])))
    if (length(ais)>0){
      runstosoma = unlist(lapply(ais, function(x) length(unlist(shortest_paths(n, x, to = root)$vpath))))
      ais = ais[match(min(runstosoma),runstosoma)]
    }
    downstream = suppressWarnings(unique(unlist(shortest_paths(n, ais, to = leaves, mode = "in")$vpath)))
    nodes[as.character(downstream),"compartment"] = "axon"
    zeros = subset(rownames(nodes),nodes[,"flow.cent"]==0)
    nodes[zeros,"compartment"] = "null"
    if(!is.null(primary.dendrite)){
      highs = subset(rownames(nodes),nodes[,"flow.cent"]>=primary.dendrite*max(nodes[,"flow.cent"]))
      nodes[as.character(highs),"compartment"] = "primary dendrite"
    }
    # Calculate segregation score
    dendrites = subset(nodes, nodes$compartment == "dendrite")
    dendrites.post = nrow(subset(dendrites,dendrites$post>0))
    dendrites.pre = nrow(subset(dendrites,dendrites$pre>0))
    dendrites.both = dendrites.post + dendrites.pre
    dendrites.pi = dendrites.post/(dendrites.pre+dendrites.post)
    dendrites.si = -(dendrites.pi*log(dendrites.pi)+(1-dendrites.pi)*log(1-dendrites.pi))
    axon = subset(nodes, nodes$compartment == "axon")
    axon.post = nrow(subset(axon,axon$post>0))
    axon.pre = nrow(subset(axon,axon$pre>0))
    axon.both = axon.post + axon.pre
    axon.pi = axon.post/(axon.pre+axon.post)
    axon.si = -(axon.pi*log(axon.pi)+(1-axon.pi)*log(1-axon.pi))
    entropy.score = (1/(dendrites.both+axon.both))*(axon.si*axon.both+dendrites.si*dendrites.both)
    both.comps = nrow(subset(nodes,nodes$post>0))/(dendrites.both+axon.both)
    control.score = -(both.comps*log(both.comps)+(1-both.comps)*log(1-both.comps))
    segregation.index = 1 - (entropy.score/control.score)
    # Add new data to object
    neuron$d = nodes
    neuron$segregation.index = segregation.index
    neuron$type = ifelse(segregation.index > 0.05, "interneuron", "PN")
    nlist = c(nlist, as.neuronlist(neuron))
  }
  attr(nlist,'df') = meta
  return(nlist)
}

  
#
plot3d.split = function(someneuronlist, col = c("blue", "orange", "purple","green","pink"), WithConnectors = T, WithNodes = F, soma = T, highflow = F){
  for (neuron in 1:length(someneuronlist)){
    neuron = someneuronlist[[1]]
    if(is.null(neuron$d$flow.cent)){
      warning("No flow centrality calculated, dropping neuron")
      break
    }
    dendrites.v = subset(rownames(neuron$d), neuron$d$compartment == "dendrite")
    axon.v = subset(rownames(neuron$d), neuron$d$compartment == "axon")
    nulls.v = subset(rownames(neuron$d), neuron$d$compartment == "null")
    p.d.v = subset(rownames(neuron$d), neuron$d$compartment == "primary dendrite")
    dendrites = prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, nulls.v, p.d.v)))
    axon = prune_vertices(neuron, verticestoprune = as.integer(c(nulls.v, dendrites.v, p.d.v)))
    nulls = prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, dendrites.v, p.d.v)))
    p.d = prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, dendrites.v, nulls.v)))
    plot3d(dendrites, col = col[1], WithNodes = WithNodes)
    plot3d(axon, col = col[2], WithNodes = WithNodes)
    plot3d(nulls, col = col[3], WithNodes = WithNodes, soma = soma)
    plot3d(p.d, col = col[4], WithNodes = WithNodes)
    if (WithConnectors == T){
      points3d(subset(xyzmatrix(neuron$d),neuron$d$post>0), col = 'cyan')
      points3d(subset(xyzmatrix(neuron$d),neuron$d$pre>0), col = 'red')
    }
    if (highflow == T){
      highest = max(neuron$d[,"flow.cent"])
      s.d = sd(neuron$d[,"flow.cent"], na.rm = T)
      high = subset(neuron$d, neuron$d[,"flow.cent"] > (highest - s.d))
      points3d(nat::xyzmatrix(high), col = col[5])
    }
  }
}




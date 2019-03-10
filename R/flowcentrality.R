#' Determine dendritic/axonal by calculating flow centrality
#'
#' @description implementation of the algorithm for calculating  flow
#'   centralities from Schneider-Mizell et al. (2016). Note that the neurites() function will retrieve these clusters as separate neuron objects in a neuronlist.
#'
#' @param x a neuronlist or neuron object
#' @param mode type of flow centrality to calculate. There are three flavors:
#'   (1) centrifugal, which counts paths from proximal inputs to distal outputs;
#'   (2) centripetal, which counts paths from distal inputs to proximal outputs;
#'   and (3) the sum of both.
#' @param polypre whether to consider the number of presynapses as a multiple of
#'   the numbers of connections each makes
#' @param primary.dendrite whether to try to assign nodes to a 'primary
#'   dendrite'. Defaults to considering nodes of 0.9*maximal flow centrality.
#'   Assigning to NULL will prevent generating this compartment.
#' @param bending.flow we may need to add the 'bending flow' to all the branchpoints if looking at centripetal flow centrality
#' @param bending.flow the algorithm will assign two main neurite compartments, which as per SWC format will be indicates as either axon (Label =2)
#' or dendrite (Label = 3) in the returned objects, at neuron$d$Label.
#' This assignment can be based which compartment contains the most postsynapses ("postsynapses") or presynapses ("presynapses"),
#' or the Euclidean distance of its first branch point from the primary branch point (i.e. the first branch point from the soma) ("distance").
#' @param ... additional arguments passed to methods.
#'
#' @details From Schneider-Mizell et al. (2016): "We use flow centrality for
#'   four purposes. First, to split an arbor into axon and dendrite at the
#'   maximum centrifugal SFC, which is a preliminary step for computing the
#'   segregation index, for expressing all kinds of connectivity edges (e.g.
#'   axo-axonic, dendro-dendritic) in the wiring diagram, or for rendering the
#'   arbor in 3d with differently colored regions. Second, to quantitatively
#'   estimate the cable distance between the axon terminals and dendritic arbor
#'   by measuring the amount of cable with the maximum centrifugal SFC value.
#'   Third, to measure the cable length of the main dendritic shafts using
#'   centripetal SFC, which applies only to insect neurons with at least one
#'   output syn- apse in their dendritic arbor. And fourth, to weigh the color
#'   of each skeleton node in a 3d view, providing a characteristic signature of
#'   the arbor that enables subjective evaluation of its identity."
#'
#' @references Schneider-Mizell, C. M., Gerhard, S., Longair, M., Kazimiers, T.,
#'   Li, F., Zwart, M. F., … Cardona, A. (2015). Quantitative neuroanatomy for
#'   connectomics in Drosophila. bioRxiv, 026617. http://doi.org/10.1101/026617
#'
#' @return the neuron or neuron list object inputted, with centripetal flow
#'   centrality information added to neuron$d, a segregation index score and
#'   estimation of neuronal type (interneuron or PN) based on this score (>0.05
#'   = PN).
#' @export
#' @seealso \code{\link{seesplit3d}} \code{\link{get.synapses}} \code{\link{neurites}}
flow.centrality <-function(x, mode = c("sum","centrifugal","centripetal"), polypre = TRUE, primary.dendrite = 0.9, bending.flow = FALSE,split = c("postsynapses","presynapses","distance"),...) UseMethod("flow.centrality")

#' @export
#' @rdname flow.centrality
flow.centrality.neuron <- function(x, mode = c("sum","centrifugal","centripetal"),
                                   polypre = TRUE, primary.dendrite = 0.9,
                                   bending.flow = FALSE, split = c("postsynapses","presynapses","distance"), ...){
  # prune Strahler first...and use segmentgraph?
  split = match.arg(split)
  mode = match.arg(mode)
  # Generate ngraph object
  x$d$Label = 0
  el = x$d[x$d$Parent != -1, c("Parent", "PointNo")] # Get list of soma=leaf directed connections
  n = nat::ngraph(data.matrix(el[,2:1]), x$d$PointNo, directed = TRUE, xyz = nat::xyzmatrix(x$d),
                  diam = x$d$W) # Make ngraph object, but for centripetal, invert the el list
  # Get comprehensive paths list
  leaves = which(igraph::degree(n, v = igraph::V(n), mode = "in")==0, useNames = T)
  root= which(igraph::degree(n, v = igraph::V(n), mode = "out")==0, useNames = T)
  segs = x$SegList # Get the segment list
  nodes = x$d # get neuron's node data
  nodes[,"post"] <- 0 # raw synapse number at this location
  nodes[,"pre"] <- 0 # raw synapse number at this location
  nodes[,"up.syns.in"] <- 0 # Here, we are going to culmulatively count synapses
  nodes[,"up.syns.out"] <- 0 # Same but for outputs
  nodes[,"flow.cent"] <- 0 # We'll addd the score for each node here
  nodes[,"Label"] <- 3
  nodes = nodes[unlist(c(root, lapply(segs, function (x) x[-1]))),]
  syns.in = x$connectors[x$connectors[,3]==1,][,1]
  if (polypre == T){
    pres = x$connectors[x$connectors[,3]==0,][,2]
    pre.cons = catmaid::catmaid_get_connectors(pres)$connector_id
    pre.cons = c(pre.cons,pres[!pres%in%pre.cons])
    syns.out = x$connectors[,1][match(pre.cons, x$connectors[,2])]
  }else{
    syns.out = x$connectors[x$connectors[,3]==0,][,1]
  }
  # Rearrange so nodes are the indices and we can count no. of synapses to which nodes connect  -(count things multiples???)
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
  in.bps.child = tapply(in.bps,names(in.bps),function(x) ifelse(names(x)[1]%in%names(nodes.in),sum(x)-(nodes.in[names(x)[1]]*length(x)),sum(x))) # Add all the branch point values together for the immediate child fragments, and only coutn synapses on the BP itself ONCE
  out.bps.child = tapply(out.bps,names(out.bps),function(x) ifelse(names(x)[1]%in%names(nodes.out),sum(x)-(nodes.out[names(x)[1]]*length(x)),sum(x)))
  nodes[names(in.bps.child),"up.syns.in"] = in.bps.child
  nodes[names(out.bps.child),"up.syns.out"] = out.bps.child # However! If there is a synapse on the BP itself, it is counted twice...
  #Old way of preventing over counting synapses from branchpoints
  #in.bps[names(in.bps)%in%names(nodes.in)] = in.bps[names(in.bps)%in%names(nodes.in)] - nodes[names(in.bps)[names(in.bps)%in%names(nodes.in)],"post"]
  #out.bps[names(out.bps)%in%names(nodes.out)] = out.bps[names(out.bps)%in%names(nodes.out)] - nodes[names(out.bps)[names(out.bps)%in%names(nodes.out)],"pre"]
  # Propagate scores from branch nodes
  #plot to check: points3d(xyzmatrix(nodes),col=rbPal(max(nodes$up.syns.in+1))[nodes$up.syns.in+1])
  bps = as.numeric(names(in.bps.child))
  for (i in 1:length(bps)){
    bp = bps[i]
    new.in = new.out = c(rep(0,nrow(nodes)))
    names(new.in) = names(new.out) = rownames(nodes)
    vertices = unlist(igraph::shortest_paths(n, bp, to = root)$vpath)[-1]
    new.in[as.character(vertices)] = new.in[as.character(vertices)] + in.bps.child[i]
    new.out[as.character(vertices)] = new.out[as.character(vertices)] + out.bps.child[i]
    nodes[,"up.syns.in"] = nodes[,"up.syns.in"]+ new.in
    nodes[,"up.syns.out"] = nodes[,"up.syns.out"] + new.out
  }
  # Calculate flow centrality
  in.total = nodes[1,"up.syns.in"] = length(point.no.in)
  out.total = nodes[1,"up.syns.out"] = length(point.no.out)
  if(mode[1] == "centrifugal"){nodes[,"flow.cent"] = (in.total - nodes[,"up.syns.in"])*nodes[,"up.syns.out"]}
  if(mode[1] == "centripetal"){nodes[,"flow.cent"] = (out.total - nodes[,"up.syns.out"])*nodes[,"up.syns.in"]}
  if(mode[1] == "sum"){nodes[,"flow.cent"] = ((in.total - nodes[,"up.syns.in"])*nodes[,"up.syns.out"]) + ((out.total - nodes[,"up.syns.out"])*nodes[,"up.syns.in"])}
  nodes = nodes[order(as.numeric(rownames(nodes))),]
  # Calculate flow centrality at branch points
  if(bending.flow){
    for(bp in bps){
      # We need to add the 'bending flow' to all the branchpoints if looking at centripetal flow centrality
      down = unlist(igraph::ego(n, 1, nodes = bp, mode = "in",mindist=0))[-1]
      bending.flow = centrifugal.bending.flow = c()
      for(u in down){
        this.seg.posts = nodes[u,]$down.syns.in
        other.segs.pre = nodes[down[!down==u],]$down.syns.out
        bending.flow = c(bending.flow,sum(this.seg.posts*other.segs.pre))
      }
      nodes[bp,"flow.cent"] = nodes[bp,"flow.cent"] + bending.flow
    }
  }
  ais = which(apply(nodes, 1, function(x) x["flow.cent"] == max(nodes[,"flow.cent"])))
  if (length(ais)>0){
    runstosoma = unlist(lapply(ais, function(x) length(unlist(igraph::shortest_paths(n, x, to = root)$vpath))))
    ais = ais[match(min(runstosoma),runstosoma)]
  }
  downstream = suppressWarnings(unique(unlist(igraph::shortest_paths(n, ais, to = leaves, mode = "in")$vpath)))
  upstream = rownames(nodes)[!rownames(nodes)%in%downstream]
  if(bending.flow){
    # If the split point is the primary branch point then both axon and dendrite are upstream!
    if(nodes[ais,]$up.syns.in==0&nodes[ais,]$up.syns.out==0){
      down = unlist(igraph::ego(n, 1, nodes = bp, mode = "in",mindist=0))[-1]
      ais = down[1]
      downstream = suppressWarnings(unique(unlist(igraph::shortest_paths(n, ais, to = leaves, mode = "in")$vpath)))
      upstream = rownames(nodes)[!rownames(nodes)%in%downstream]
    }
  }
  igraph::V(n)$name = igraph::V(n)
  # Work out the primary neurite and empty leaf Labels
  x.pruned = nat::prune_strahler(x=x,orderstoprune = 1:2)
  pnt = suppressWarnings(primary.neurite(x.pruned,keep.pnt = TRUE, resample = FALSE))
  p.n = rownames(nodes)[match(pnt$d$PointNo,nodes$PointNo)]
  ### OLD WAY OF DETECTING PNT ###
  # zeros = subset(rownames(nodes),nodes[,"flow.cent"]==0)
  # remove = rownames(nodes)[!rownames(nodes)%in%zeros]
  # nn = igraph::delete_vertices(n, v = as.character(remove))
  # clust = igraph::clusters(nn)
  # soma.group = clust$membership[names(clust$membership)==root]
  # p.n = names(clust$membership)[clust$membership==soma.group]
  nodes[p.n,"Label"] = 7
  if(!is.null(primary.dendrite)){
    highs = subset(rownames(nodes),nodes[,"flow.cent"]>=primary.dendrite*max(nodes[,"flow.cent"]))
    nodes[as.character(highs),"Label"] = 4
  }else{
    primary.dendrite = 0.9
    highs = subset(rownames(nodes),nodes[,"flow.cent"]>=primary.dendrite*max(nodes[,"flow.cent"]))
  }
  ### Find primary and secondary branch points ###
  # Find tract parent node for downstream
  downstream.unclassed = downstream[!downstream%in%c(p.n,highs)]
  remove = rownames(nodes)[!rownames(nodes)%in%downstream.unclassed]
  downstream.g = igraph::delete_vertices(n, v = as.character(remove))
  main1 = igraph::components(downstream.g) # Minimally the 'two' main branches, as a primary dendrite node will bisect the singular main branch
  main1 = names(main1$membership[main1$membership%in%1])
  nodes.downstream = nodes[as.character(main1),]
  tract.parent = unique(nodes.downstream$Parent[!nodes.downstream$Parent%in%nodes.downstream$PointNo])
  downstream.tract.parent = match(tract.parent,nodes$PointNo)
  if(sum(match(tract.parent,nodes$PointNo)%in%p.n)>0){
    bps.all = rownames(nodes)[match(as.numeric(nat::branchpoints(nodes)),nodes$PointNo)]
    bps.downstream = bps.all[bps.all%in%downstream.unclassed]
    runstoprimarybranchpoint = unlist(lapply(bps.downstream, function(x) length(unlist(suppressWarnings(igraph::shortest_paths(n, to = downstream.tract.parent, from = x)$vpath)))))
    downstream.tract.parent = bps.downstream[which.min(runstoprimarybranchpoint)]
  }
  downstream.tract.parent = nodes[downstream.tract.parent,]
  # Find tract parent node for upstream
  upstream.unclassed = upstream[!upstream%in%c(p.n,highs)]
  remove = rownames(nodes)[!rownames(nodes)%in%upstream.unclassed]
  upstream.g = igraph::delete_vertices(n, v = as.character(remove))
  main1 = igraph::components(upstream.g) # Minimally the 'two' main branches, as a primary dendrite node will bisect the singular main branch
  main1 = names(main1$membership[main1$membership%in%1])
  nodes.upstream = nodes[as.character(main1),]
  tract.parent = unique(nodes.upstream$Parent[!nodes.upstream$Parent%in%nodes.upstream$PointNo])
  upstream.tract.parent = match(tract.parent,nodes$PointNo)
  if(sum(match(tract.parent,nodes$PointNo)%in%p.n)>0){
    bps.all = rownames(nodes)[match(as.numeric(nat::branchpoints(nodes)),nodes$PointNo)]
    bps.upstream = bps.all[bps.all%in%upstream.unclassed]
    runstoprimarybranchpoint = unlist(lapply(bps.upstream, function(x) length(unlist(suppressWarnings(igraph::shortest_paths(n, to = upstream.tract.parent, from = x)$vpath)))))
    upstream.tract.parent = bps.upstream[which.min(runstoprimarybranchpoint)]
  }
  upstream.tract.parent = nodes[upstream.tract.parent,]
  # Find primary branchpoint
  neurite.nodes = nodes[!rownames(nodes)%in%p.n,]
  p.n.PointNo = nodes[p.n,"PointNo"]
  primary.branch.point = p.n[p.n.PointNo%in%neurite.nodes$Parent]
  ### Assign putative axonic and dendritic compartments ###
  if(grepl("synapses",split)){
    synapse.choice = gsub("synapses","",split)
    message(split)
    choice = sum(nodes[as.character(downstream.unclassed),synapse.choice]) < sum(nodes[as.character(upstream.unclassed),synapse.choice])
    if (choice){
      nodes[as.character(downstream.unclassed),"Label"] = 2
    }else if(!choice) {
      nodes[as.character(upstream.unclassed),"Label"] = 2
    } else{
      split=="distance"
      warning("synapse numbers are the same, splitting based on branch point distances to primary branchpoint")
    }
  }
  if(split=="distance"){ # Distance from primary branchpoint
    primary.branch.point.xyz = as.matrix(nat::xyzmatrix(nodes[primary.branch.point,]))
    secondary.branch.points.xyz = nat::xyzmatrix(rbind(upstream.tract.parent,downstream.tract.parent))
    dist.upstream.to.primary.branchpoint = length(unlist(igraph::shortest_paths(n, to = primary.branch.point, from = rownames(upstream.tract.parent))$vpath))# Or euclidean distance: nabor::knn(query=primary.branch.point.xyz,data= nat::xyzmatrix(upstream.tract.parent),k=1)$nn.dists
    dist.downstream.to.primary.branchpoint = length(unlist(igraph::shortest_paths(n, to = primary.branch.point, from = rownames(downstream.tract.parent))$vpath))# Or euclidean distance: nabor::knn(query=primary.branch.point.xyz,data= nat::xyzmatrix(downstream.tract.parent),k=1)$nn.dists
    if(dist.upstream.to.primary.branchpoint<dist.downstream.to.primary.branchpoint){
      nodes[as.character(downstream.unclassed),"Label"] = 2
    }else if(dist.upstream.to.primary.branchpoint>dist.downstream.to.primary.branchpoint){
      nodes[as.character(upstream.unclassed),"Label"] = 2
    }else{
      warning("branch point distances are the same, splitting based on postsynapses")
      choice = sum(nodes[as.character(downstream.unclassed),"post"]) < sum(nodes[as.character(upstream.unclassed),"post"])
      if (choice){
        nodes[as.character(downstream.unclassed),"Label"] = 2
      }else{
        nodes[as.character(upstream.unclassed),"Label"] = 2
      }
    }
  }
  ### Calculate segregation score ###
  dendrites = subset(nodes, nodes$Label == 3)
  dendrites.post = sum(subset(dendrites$post,dendrites$post>0))
  dendrites.pre = sum(subset(dendrites$pre,dendrites$pre>0))
  dendrites.both = dendrites.post + dendrites.pre
  dendrites.pi = dendrites.post/dendrites.both
  dendrites.si = -(dendrites.pi*log(dendrites.pi)+(1-dendrites.pi)*log(1-dendrites.pi))
  if(is.nan(dendrites.si)){dendrites.si = 0}
  axon = subset(nodes, nodes$Label == 2)
  axon.post = sum(subset(axon$post,axon$post>0))
  axon.pre = sum(subset(axon$pre,axon$pre>0))
  axon.both = axon.post + axon.pre
  axon.pi = axon.post/axon.both
  axon.si = -(axon.pi*log(axon.pi)+(1-axon.pi)*log(1-axon.pi))
  if(is.nan(axon.si)){axon.si = 0}
  entropy.score = (1/(dendrites.both+axon.both))*((axon.si*axon.both)+(dendrites.si*dendrites.both))
  both.comps = (dendrites.post+axon.post)/(dendrites.both+axon.both)
  # both.comps = sum(nodes$post)/(sum(nodes$pre+sum(nodes$post)))
  control.score = -(both.comps*log(both.comps)+(1-both.comps)*log(1-both.comps))
  segregation.index = 1 - (entropy.score/control.score)
  if(is.na(segregation.index)) { segregation.index = 0 }
  # Add new data to object
  x$d = nodes
  x$AD.segregation.index = segregation.index
  x$type = ifelse(segregation.index < 0.05, "interneuron", "PN")
  x$primary.branch.point = as.numeric(primary.branch.point)
  x$secondary.branch.points = as.numeric(c(downstream.tract.parent$PointNo,upstream.tract.parent$PointNo))
  x$max.flow.centrality = as.numeric(ais)
  x
}

#' @export
#' @rdname flow.centrality
flow.centrality.neuronlist <- function(x, mode = c("sum","centrifugal","centripetal"), polypre = T, primary.dendrite = 0.9, bending.flow = FALSE,split = c("postsynapses","presynapses","distance"),...){
  neurons = nat::nlapply(x, flow.centrality, mode = mode, polypre = polypre, primary.dendrite = primary.dendrite, OmitFailures = T, split = split, ...)
  neurons
}


#' Plot neurons split up by flow centrality
#'
#' @param someneuronlist a neuronlist or neuron object that has been modified by flow.centrality
#' @param col colours of sections. Defaults to orange or axons, green for primary dendrite, blue for dendrites and pink for nodes with no flow.
#' @param WithConnectors whether ot plot the anatomical location of pre (red) and post (cyan) synapses.
#' @param soma whether to plot a soma, and what the radius should be
#' @param WithNodes whether to plot branch points
#' @param lwd Line width (default 1)
#' @param radius For connectors and axon-dendrite split node (default 1)
#' @param highflow whether to plot the nodes of highest (with in one standard deviation less than maximum) flow centrality (pink points)
#' @param Verbose logical indicating that info about each selected neuron should be printed (default TRUE)
#' @param Wait logical indicating that there should be a pause between each displayed neuron
#' @param sleep time to pause between each displayed neuron when Wait=TRUE
#' @param extrafun an optional function called when each neuron is plotted, with two arguments: the current neuron name and the current selected neurons
#' @param selected_file an optional path to a yaml file that already contains a selection
#' @param selected_col the color in which selected neurons (such as those specified in selected_file) should be plotted
#' @param yaml a logical indicating that selections should be saved to disk in (human-readable) yaml rather than (machine-readable) rda format
#' @param ... additional arguments passed to methods.
#'
#' @return Plots coloured neuron(s)
#' @export
#' @seealso \code{\link{flow.centrality}} \code{\link{get.synapses}}
#' @importFrom stats sd
seesplit3d = function(someneuronlist, col = c("blue", "orange", "purple","green", "grey", "pink"), splitnode = FALSE,WithConnectors = TRUE, WithNodes = F, soma = 100, highflow = F, lwd = 1, radius = 1, ...){
  someneuronlist = nat::as.neuronlist(someneuronlist)
  for (n in 1:length(someneuronlist)){
    neuron = someneuronlist[[n]]
    if(is.null(neuron$d$flow.cent)){
      warning("No flow centrality calculated, dropping neuron")
      break
    }
    dendrites.v = subset(rownames(neuron$d), neuron$d$Label == 3)
    axon.v = subset(rownames(neuron$d), neuron$d$Label == 2)
    #nulls.v = subset(rownames(neuron$d), neuron$d$Label == 0)
    p.d.v = subset(rownames(neuron$d), neuron$d$Label == 4)
    p.n.v = subset(rownames(neuron$d), neuron$d$Label == 7)
    dendrites = nat::prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, p.d.v, p.n.v)))
    axon = nat::prune_vertices(neuron, verticestoprune = as.integer(c(dendrites.v, p.d.v, p.n.v)))
    #nulls = nat::prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, dendrites.v, p.d.v, p.n.v)))
    p.d = prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, dendrites.v, p.n.v)))
    p.n = prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, dendrites.v, p.d.v)))
    rgl::plot3d(dendrites, col = col[1], WithNodes = WithNodes, lwd = lwd,...)
    rgl::plot3d(axon, col = col[2], WithNodes = WithNodes, soma = FALSE, lwd = lwd,...)
    rgl::plot3d(p.n, col = col[3], WithNodes = WithNodes, soma = soma, lwd = lwd,...)
    rgl::plot3d(p.d, col = col[4], WithNodes = WithNodes, soma = FALSE, lwd = lwd,...)
    #rgl::plot3d(nulls, col = col[5], WithNodes = WithNodes, soma = FALSE, lwd = lwd)
    #rgl::plot3d(neuron, col = col[3], WithNodes = WithNodes, soma = soma)
    if (WithConnectors){
      rgl::spheres3d(subset(xyzmatrix(neuron$d),neuron$d$post>0), col = 'cyan', radius = radius,...)
      rgl::spheres3d(subset(xyzmatrix(neuron$d),neuron$d$pre>0), col = 'red', radius = radius,...)
    }
    if (highflow == T){
      highest = max(neuron$d[,"flow.cent"])
      s.d = sd(neuron$d[,"flow.cent"], na.rm = T)
      high = subset(neuron$d, neuron$d[,"flow.cent"] > (highest - s.d))
      rgl::points3d(nat::xyzmatrix(high), col = col[6],...)
    }
    if(splitnode==T){
      ais = which(apply(neuron$d, 1, function(x) x["flow.cent"] == max(neuron$d[,"flow.cent"])))
      rgl::spheres3d(nat::xyzmatrix(neuron$d[ais,]),radius=radius,col="magenta",...)
    }
  }
}


#' @export
#' @rdname seesplit3d
splitscan <- function (someneuronlist, col = c("blue", "orange", "purple","green", "grey", "pink"), WithConnectors = T, WithNodes = F, soma = 100, highflow = F, Verbose = T, Wait = T,
          sleep = 0.1, extrafun = NULL, selected_file = NULL, selected_col = "black",
          yaml = TRUE, ...)
{
  if(!requireNamespace('yaml', quietly = TRUE))
  stop("Suggested package yaml is required to use this function!")
  if (is.neuronlist(someneuronlist)) {
    db = someneuronlist
    neurons = as.data.frame(db)$name
  }
  frames <- length(neurons)
  selected <- character()
  i <- 1
  if (!is.null(selected_file) && file.exists(selected_file)) {
    selected <- yaml::yaml.load_file(selected_file)
    if (!all(names(selected) %in% neurons))
      stop("Mismatch between selection file and neurons.")
  }
  savetodisk <- function(selected, selected_file) {
    if (is.null(selected_file))
      selected_file <- file.choose(new = TRUE)
    if (yaml) {
      if (!grepl("\\.yaml$", selected_file))
        selected_file <- paste(selected_file, sep = "",
                               ".yaml")
      message("Saving selection to disk as ", selected_file,
              ".")
      writeLines(yaml::as.yaml(selected), con = selected_file)
    }
    else {
      if (!grepl("\\.rda$", selected_file))
        selected_file <- paste(selected_file, sep = "",
                               ".rda")
      save(selected, file = selected_file)
      message("Saving selection to disk as ", selected_file)
    }
    selected_file
  }
  chc <- NULL
  while (TRUE) {
    if (i > length(neurons) || i < 1)
      break
    n <- neurons[i]
    cat("Current neuron:", n, "(", i, "/", length(neurons),
        ")\n")
    pl <- seesplit3d(someneuronlist[i], col = col, WithConnectors = WithConnectors, WithNodes = WithNodes, soma = soma, highflow = highflow)
    print(someneuronlist[[i]]$segregation.index)
    more_rgl_ids <- list()
    if (!is.null(extrafun))
      more_rgl_ids <- extrafun(n, selected = selected)
    if (Wait) {
      chc <- readline("Return to continue, b to go back, s to select, d [save to disk], t to stop, c to cancel (without returning a selection): ")
      if (chc == "c" || chc == "t") {
        sapply(pl, rgl::rgl.pop, type = "shape")
        sapply(more_rgl_ids, rgl::rgl.pop, type = "shape")
        break
      }
      if (chc == "s") {
        if (n %in% selected) {
          message("Deselected: ", n)
          selected <- setdiff(selected, n)
        }
        else selected <- union(selected, n)
      }
      if (chc == "b")
        i <- i - 1
      else if (chc == "d")
        savetodisk(selected, selected_file)
      else i <- i + 1
    }
    else {
      Sys.sleep(sleep)
      i <- i + 1
    }
    sapply(pl, rgl::rgl.pop, type = "shape")
    sapply(more_rgl_ids, rgl::rgl.pop, type = "shape")
    rgl::clear3d()
  }
  if (is.null(chc) || chc == "c")
    return(NULL)
  if (!is.null(selected_file))
    savetodisk(selected, selected_file)
  selected
}


# nopen3d()
# for(i in 1:length(l)){
#   message(i)
#   seesplit3d(l[[i]],lwd=2,soma=500)
#   y = ll[[i]]
#   xyzmatrix(y) = xyzmatrix(y)+1000
#   seesplit3d(y,lwd=2,soma=500, col = c("cyan", "yellow", "magenta","darkgreen", "darkgrey", "deeppink"))
#   p = readline("Next? ")
#   clear3d()
# }

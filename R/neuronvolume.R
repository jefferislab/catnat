#' Get neuron volumes and stitch them into a cohesive neuron, or neuronal compartment
#'
#' @description Stitch together volumes for a single neuron, based on the Google FAFB segmentation by Peter Li. If the neuron given has its microtubules marked (catnat::mark.microtubules), Strahler orders (catnat::assign_strahler) or
#' axon/dendrite (catnat::flow.centrality) then this information can be used to estimate a volume that over these neuronal sub-compartments.
#' @param neuron a neuron object that matches the volumes given
#' @param volumes a list of mesh3d objects, retrieved using fafbseg::read_brainmaps_meshes
#' @param voxelSize a single number, which is used to downsample the starting meshes, resample neuron, and choose an alpha value for alphashape3d::ashape3d()
#' @param downsample.factor after the Google volumes have been resampled uniformly at voxelSize, points are randomly removed from the neuron cloud
#' (no.points=no.points/downsample.factor) in order to make alphashape generation by alphashape3d::ashape manageable
#' @param map if TRUE, instead of using the volumes argument, map_fafbsegs_to_neuron and fafbseg::read_brainmaps_meshes are used to retrieve Google 3D segmentations as mesh3d objects
#' @param soma if TRUE, the soma (root point for neuron) will be identified
#' @param node.match how many nodes of each neuron in someneuronlist, need to be within a auto segmented volume, for it to be said to match.
#' These nodes all need to be consecutive, in the sense that they must be in the same segment or a branch from that segment. I.e. If a neuron matches with a volume
#' 5 times at diverse points across it arbour, this is thought to be a non-match with a large, proximal auto-traced segment.
#' need be in the volumetric Google FAFB segmentation for a Neuroglancer fragment, for that fragment to be returned
#' @param resample.neuron if TRUE, neuron is resampled (nat:resample), stepsize = voxelSize
#' @param resample.volume if TRUE, the final mesh3D object is uniformly resampled using Rvcg::vcgUniformRemesh, voxelSize = voxelSize
#' @param smooth if TRUE, the final mesh will be smoothed, using
#' @param smooth.type character: select smoothing algorithm. Available are "taubin", "laplace", "HClaplace","fujiLaplace", "angWeight" (and any sensible abbreviations). See Rvcg::vcgSmooth
#' @param lambda numeric: parameter for Taubin smooth. See Rvcg::vcgSmooth
#' @param mu numeric:parameter for Taubin smooth. See Rvcg::vcgSmooth
#' @param delta numeric: parameter for Scale dependent laplacian smoothing and maximum allowed angle (in radians) for deviation between normals Laplacian (surface preserving). See Rvcg::vcgSmooth
#' @param conn CATMAID connection
#' @param ... methods passed to Rvcg functions: vcgUniformRemesh, vcgSmooth, vcgClost
#' @return a 'neuronvolume' object, with the same structure as a neuron object (see the nat package), but which contains a mesh3D object for the neuron and other volume related data as a list at neuron$volume
#' @export
#' @rdname fafb_segs_stitch_volumes
fafb_segs_stitch_volumes <- function(neuron, volumes = NULL, map = TRUE, voxelSize = 50, downsample.factor = 12,
                                     soma = TRUE, node.match = 4, smooth = FALSE, resample.neuron = TRUE, resample.volume = FALSE,
                                     smooth.type=c("taubin", "laplace", "HClaplace", "fujiLaplace","angWeight", "surfPreserveLaplace"),
                                     lambda = 0.5, mu = -0.53, delta = 0.1, conn = NULL, ...){
  if(!requireNamespace('pbapply', quietly = TRUE))
    stop("Please install suggested fafbseg pbapply")
  smooth.type = match.arg(smooth.type)
  downsample_vol <- function(vol,voxelSize = 50, ...){
    v = Rvcg::vcgUniformRemesh(vol, voxelSize = voxelSize, offset = 0, discretize = FALSE,
                               multiSample = TRUE, absDist = FALSE, mergeClost = FALSE,
                               silent = TRUE)
    v = tryCatch(Rvcg::vcgIsolated(v, facenum = NULL, diameter = NULL, split = FALSE,
                                   keep = 0, silent = TRUE),error = function(e) v)
    v
  }
  if(!nat::is.neuron(neuron)){
    if(nat::is.neuronlist(neuron)){
      message("Neuronlist object provided, only using first neuron in neuronlist")
      neuron = neuron[[1]]
    } else if(length(neuron)==1){
      message("Reading neuron from CATMAID")
      neuron = catmaid::read.neuron.catmaid(neuron, conn = conn)
    }else{
      stop("Please provide a valid neuron object or ID that can be read using catmaid::read.neuron.catmaid")
    }
  }
  if(map){ # Map neuron points onto FAFBseg volumes
    mapping = map_fafbsegs_to_neuron(neuron, node.match = node.match)
    unmapped = subset(mapping,ngl_id=="0")
    unmapped = sum(unmapped$node_hits)/nrow(nat::xyzmatrix(neuron))
    nglids = unique(mapping$ngl_id)
    nglids = nglids[nglids!=0]
    volumes = list()
    nams = c()
    pb <- utils::txtProgressBar(min = 0, max = length(nglids), style = 3)
    while(length(nglids)>0){
      message("Getting volumes for ",length(nglids), " Neuroglancer segments")
      for(ng in 1:length(nglids)){
        s = fafbseg::find_merged_segments(nglids[ng])
        vol = tryCatch(fafbseg::read_brainmaps_meshes(s), error = function(e) "read error")
        if(is.character(vol)){
          message(vol)
          next
        }
        else if(is.null(vol)|length(vol)<1){
          message("no volume to read")
          vol = NULL
          nams = c(nams,nglids[ng])
          volumes[[ng]] = vol
        }else{
          message(paste0("volume read: ",ng))
          nams = c(nams,nglids[ng])
          volumes[[ng]] = vol
        }
      }
      nglids = setdiff(nglids,nams)
      message("re-trying: ", length(nglids))
      utils::setTxtProgressBar(pb, ng)
    }
    close(pb)
  }
  if(is.null(volumes)){
    stop("Error: please provide a list of mesh3D objects associated with your neuron, or use map = TRUE")
  }
  if(resample.neuron){
    message("Resampling neuron")
    neuron = nat::resample(neuron,stepsize = voxelSize)
  }
  neuron.points = nat::xyzmatrix(neuron)
  neuron$d$in.segment = neuron$d$cut = FALSE
  neuron$d$volume = NA
  neuron$volume = list() # For collecting volume data
  message("Downsampling meshes")
  downvolumes = pbapply::pblapply(volumes,function(s) tryCatch(downsample_vol(s),error = function(e) s))
  ### FIND SOMA ###
  if(soma){
    message("Identifying soma")
    find.soma = c()
    pb <- utils::txtProgressBar(min = 0, max = length(downvolumes), style = 3)
    for(i in 1:length(downvolumes)){
      dv = downvolumes[[i]]
      p = nabor::knn(query=nat::xyzmatrix(neuron$d[nat::rootpoints(neuron),]),data=nat::xyzmatrix(dv),k=1,radius=1000)
      find.soma = c(find.soma,p$nn.dists)
      utils::setTxtProgressBar(pb, i)
    }
    find.soma[is.infinite(find.soma)] = 0
    if(sum(find.soma)==0|sum(find.soma)>1){
      children = neuron$d$PointNo[neuron$d$Parent==nat::rootpoints(neuron)]
      som.points = nat::xyzmatrix(neuron$d[c(nat::rootpoints(neuron),children),])
      soma.dv = tryCatch(fafbseg::read_brainmaps_meshes(som.points), error = function(e) "soma read error")
    }else{
      soma.dv = downvolumes[[which.max(find.soma)]]
    }
    if(is.character(soma.dv)|is.null(soma.dv)){
      soma = FALSE
    }else{
      soma.points = nat::xyzmatrix(soma.dv)
      p = nabor::knn(data=neuron.points[-1:-3,],query=soma.points,k=1,radius=1000)
      soma.ashape = alphashape3d::ashape3d(unique(soma.points[p$nn.dists>1000,]),pert = TRUE,alpha = voxelSize*50)
      soma.mesh = tryCatch(ashape2mesh3d(soma.ashape, remove.interior.points = TRUE),error = function(e) ashape2mesh3d(soma.ashape, remove.interior.points = FALSE))
      soma.dv = tryCatch(Rvcg::vcgUniformRemesh(x=soma.mesh,voxelSize = voxelSize*4),error = function(e) soma.mesh)
      downvolumes[[length(downvolumes)+1]] = soma.dv
      soma.ashape.volume = tryCatch(alphashape3d::volume_ashape3d(as3d=soma.ashape, byComponents = FALSE, indexAlpha = 1),error = function(e) NA)
      neuron$volume$volume.estimation$soma.ashape.volume = soma.ashape.volume
      neuron$volume$mesh3d$soma = soma.dv
    }
    close(pb)
  }
  ### RADIUS ESTIMATION ###
  message("Estimating node radii")
  neuron$d$radius = NA
  pb <- utils::txtProgressBar(min = 0, max = length(downvolumes),
                              style = 3)
  for (i in 1:length(downvolumes)) {
    dv = downvolumes[[i]]
    p = tryCatch(nat::pointsinside(neuron.points, surf = dv),
                 error = function(e) FALSE)
    neuron$d$in.segment = (neuron$d$in.segment + p) > 0
    neuron$d$volume[p] = i
    dv = nat::xyzmatrix(dv)
    if (sum(p) > 0) {
      query = neuron.points[p, ]
      if (sum(p) == 1) {
        query = matrix(neuron.points[p, ], ncol = 3)
      }
      near = nabor::knn(query = query, data = dv, k = ifelse(300 >
                                                               nrow(dv), nrow(dv), 300), radius = 5000)
      r = apply(near$nn.dists, 1, mean)
      r[is.infinite(r)] = 5000
      neuron$d$radius[p] = r
    }
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  ### INTERPLATE OVER UNMAPPED SEGMENTS ###
  message("Finding transitions between mapped and unmapped segments")
  pb <- utils::txtProgressBar(min = 0, max = nrow(neuron$d)-1, style = 3)
  neuron$d$in.segment[1] = TRUE
  for(s in 1:nrow(neuron$d)){
    neuron$d$cut[s] = !identical(neuron$d$in.segment[neuron$d$Parent==neuron$d$PointNo[s]][1],neuron$d$in.segment[s])
    utils::setTxtProgressBar(pb, s)
  }
  close(pb)
  message("Interpolating radii of unmapped points")
  neuron$d$cut.upstream = (neuron$d$cut+neuron$d$in.segment)==2 # Each unmapped length, get the most upstream node that is in a segment.
  neuron$d$cut.downstream = neuron$d$cut&!neuron$d$in.segment # And the most downstream that is not
  pb <-  utils::txtProgressBar(min = 0, max = nrow(neuron$d), style = 3)
  for(c in 1:nrow(neuron$d)){
    if (neuron$d$cut.downstream[c]){
      r = NA
      while(is.na(r)){
        children = subset(neuron$d,Parent==c)
        if(nrow(children)==0){
          r = 40 # Data z-resolution
        }else{
          r = max(children$radius,na.rm=TRUE)
          if(is.na(r)){
            children = subset(neuron$d,Parent%in%children$PointNo)
          }
        }
        neuron$d$radius[c] = r
      }
    }
    utils::setTxtProgressBar(pb, c)
  }
  close(pb)
  ng = nat::as.ngraph(neuron$d)
  u = neuron$d$PointNo[neuron$d$cut.upstream]
  d = neuron$d$PointNo[neuron$d$cut.downstream]
  eps = nat::endpoints(neuron)
  u = u[!u%in%eps]
  d = c(d,u[u%in%eps])
  distances =  igraph::distances(graph=ng, v = u, to=d, mode = "out",weights = NULL)
  distances[is.infinite(distances)]= NA
  pb <-  utils::txtProgressBar(min = 0, max = nrow(distances), style = 3)
  done.d = c()
  for(i in 1:nrow(distances)){
    from = as.numeric(rownames(distances)[i])
    to = as.numeric(names(which.min(distances[i,])))
    done.d = c(done.d,to)
    if(length(from)>0){
      path = as.vector(igraph::get.shortest.paths(graph=ng, from = from, to=to, mode = "out",weights = NULL)$vpath[[1]])
      radii = neuron$d$radius[path]
      if(sum(is.na(radii))>0){
        radii.first = min(which(is.na(radii)))-1
        r = radii[radii.first]
        path = path[radii.first:length(path)]
        radii = radii[radii.first:length(radii)]
      }else{
        r = neuron$d$radius[from]
      }
      radii[is.na(radii)] = 0
      u.propagate = seq(from=r,to=0, length.out = length(path))
      d.propagate = seq(to=neuron$d$radius[to],from=0, length.out = length(path))
      propagate = u.propagate+d.propagate
      propagate = sapply(seq_along(propagate),function(x) max(c(propagate[x],radii[x]),na.rm=TRUE))
      neuron$d$radius[path] = propagate
    }
    utils::setTxtProgressBar(pb, i)
  }
  if(sum(!d%in%done.d)>0){
    for(i in d[!d%in%done.d]){
      to = i
      from = as.numeric(names(which.min(distances[,as.character(i)])))
      if(length(from)>0){
        path = as.vector(igraph::get.shortest.paths(graph=ng, from = from, to=to, mode = "out",weights = NULL)$vpath[[1]])
        radii = neuron$d$radius[path]
        if(sum(is.na(radii))>0){
          radii.first = min(which(is.na(radii)))-1
          r = radii[radii.first]
          path = path[radii.first:length(path)]
          radii = radii[radii.first:length(radii)]
        }else{
          r = neuron$d$radius[from]
        }
        radii[is.na(radii)] = 0
        u.propagate = seq(from=r,to=0, length.out = length(path))
        d.propagate = seq(to=neuron$d$radius[to],from=0, length.out = length(path))
        propagate = u.propagate+d.propagate
        propagate = sapply(seq_along(propagate),function(x) max(c(propagate[x],radii[x]),na.rm=TRUE))
        neuron$d$radius[path] = propagate
      }
    }
  }
  close(pb)
  message("Creating estimated 3D points for unmapped segments")
  interpolated.points = list()
  pointnumbers = subset(neuron$d,in.segment==FALSE)$PointNo
  pb <-  utils::txtProgressBar(min = 0, max = length(pointnumbers), style = 3)
  for(i in 1:length(pointnumbers)){
    utils::setTxtProgressBar(pb, i)
    p = pointnumbers[i]
    r = subset(neuron$d,PointNo==p)$radius
    xyz <- nat::xyzmatrix(subset(neuron$d,PointNo==p))
    if(is.na(r)){
      r = 40
    } else if (is.infinite(r)){
      r = 5000
    }
    ### Create fake points ###
    randomp <- matrix(rnorm(30), ncol=3) # generate points on a sphere
    randomp <- sqrt(3)*randomp/drop(sqrt((randomp^2) %*% rep(1, 3)))
    randomp <- (randomp / sqrt(sum(randomp^2)))
    randomp <- r*randomp
    ps <- t(apply(randomp,1,function(x) x+xyz))
    interpolated.points[[i]] = ps
  }
  interpolated.points = do.call(rbind,interpolated.points)
  close(pb)
  ### CREATE VOLUME ###
  message("Creating cohesive mesh")
  pb <-  utils::txtProgressBar(min = 0, max = 5, style = 3)
  points = do.call(rbind,lapply(downvolumes,function(s) nat::xyzmatrix(s)))
  points = unique(points)
  if(nrow(points)>15000){
    downsample.factor = 1.5*downsample.factor
  }
  edge.points = nabor::knn(query=points,data=points,k=10,radius=100)
  e = points[apply(edge.points$nn.dists,1,function(e) max(e)>voxelSize*2),]
  #outside.points = nabor::knn(query=points,data=neuron.points,k=10,radius=3000)
  #points = points[apply(outside.points$nn.dists,1,function(e) min(e)>3000|is.infinite(min(e))),] # Remove points unlikely to be associated with neuron
  p = points[sample(1:nrow(points),nrow(points)/downsample.factor),] # Downsample
  p = unique(rbind(p,e,
                   as.matrix(neuron.points),
                   as.matrix(interpolated.points)))
  utils::setTxtProgressBar(pb, 1)
  p.na = apply(p,1,function(x) sum(is.na(x)>0))
  p = p[!p.na,]
  a = alphashape3d::ashape3d(p,pert = TRUE,alpha = voxelSize*4)
  utils::setTxtProgressBar(pb, 2)
  triangles = a$triang[apply(a$triang, 1, function(x) {( any(as.numeric(x[9]) > 1))} ),][,1:3]
  vertices = unique(as.vector(unique(triangles)))
  vert = a$x[vertices,]
  utils::setTxtProgressBar(pb, 3)
  neuron.ashape = alphashape3d::ashape3d(vert,pert = TRUE,alpha = voxelSize*4)
  utils::setTxtProgressBar(pb, 4)
  mesh = ashape2mesh3d(neuron.ashape, remove.interior.points = FALSE)
  utils::setTxtProgressBar(pb, 5)
  close(pb)
  ### DOWNSAMPLE MESH ###
  if(resample.volume){
    message("Uniformly downsampling final mesh")
    mesh = Rvcg::vcgUniformRemesh(mesh, voxelSize = voxelSize, offset = 0, discretize = FALSE,
                                  multiSample = TRUE, absDist = FALSE, mergeClost = FALSE,
                                  silent = TRUE)
  }
  ### SMOOTH ###
  if(smooth){
    message("Smoothing mesh")
    mesh = Rvcg::vcgSmooth(mesh, type = smooth.type, iteration = 10, lambda = lambda,
                           mu = mu, delta = delta)
  }
  ### CALCULATE VOLUME ###
  message("Calculating alphashape volume")
  whole.ashape.volume = tryCatch(alphashape3d::volume_ashape3d(as3d=a, byComponents = FALSE, indexAlpha = 1),error = function(e) NA)
  ### ATTACH TO NEURON OBJECT ###
  neuron$volume$mesh3d$whole = mesh
  neuron$volume$volume.estimation$whole.ashape.volume = whole.ashape.volume
  neuron$volume$mapped.proportion = 1 - unmapped
  ### LABEL VERTICES ###
  message("Labelling mesh vertices")
  mesh.vertices = as.data.frame(nat::xyzmatrix(mesh))
  colnames(mesh.vertices) = c("X","Y","Z")
  match.points = nabor::knn(query=mesh.vertices,data=neuron.points,k=1)
  matched.points = neuron$d[match.points$nn.idx,]
  matched.points = matched.points[,colnames(matched.points)%in%c("PointNo","Label","microtubule","strahler_order","in.segment")]
  mesh.vertices = cbind(as.data.frame(mesh.vertices),matched.points)
  message("Creating sub-volumes...")
  if(2%in%mesh.vertices$Label){
    axon.ashape = alphashape3d::ashape3d(as.matrix(subset(mesh.vertices,Label==2)[,1:3]),pert = TRUE,alpha = voxelSize*4)
    axon.ashape.volume = tryCatch(alphashape3d::volume_ashape3d(as3d=axon.ashape, byComponents = FALSE, indexAlpha = 1),error = function(e) NA)
    neuron$volume$volume.estimation$axon.ashape.volume = axon.ashape.volume
    neuron$volume$mesh3d$axon = ashape2mesh3d(axon.ashape, remove.interior.points = FALSE)
    message("axon volume created")
  }
  if(3%in%mesh.vertices$Label){
    dendrite.ashape = alphashape3d::ashape3d(as.matrix(subset(mesh.vertices,Label==3)[,1:3]),pert = TRUE,alpha = voxelSize*4)
    dendrite.ashape.volume = tryCatch(alphashape3d::volume_ashape3d(as3d=dendrite.ashape, byComponents = FALSE, indexAlpha = 1),error = function(e) NA)
    neuron$volume$volume.estimation$dendrite.ashape.volume = dendrite.ashape.volume
    neuron$volume$mesh3d$dendrite = ashape2mesh3d(dendrite.ashape, remove.interior.points = FALSE)
    message("dendrite volume created")
  }
  if(7%in%mesh.vertices$Label){
    primary.neurite.ashape = alphashape3d::ashape3d(as.matrix(subset(mesh.vertices,Label==7)[,1:3]),pert = TRUE,alpha = voxelSize*4)
    primary.neurite.ashape.volume = tryCatch(alphashape3d::volume_ashape3d(as3d=primary.neurite.ashape, byComponents = FALSE, indexAlpha = 1),error = function(e) NA)
    neuron$volume$volume.estimation$primary.neurite.ashape.volume = primary.neurite.ashape.volume
    neuron$volume$mesh3d$primary.neurite = ashape2mesh3d(primary.neurite.ashape, remove.interior.points = FALSE)
    message("primary neurite volume created")
  }
  if(4%in%mesh.vertices$Label){
    primary.dendrite.ashape = alphashape3d::ashape3d(as.matrix(subset(mesh.vertices,Label==4)[,1:3]),pert = TRUE,alpha = voxelSize*4)
    primary.dendrite.ashape.volume = tryCatch(alphashape3d::volume_ashape3d(as3d=primary.dendrite.ashape, byComponents = FALSE, indexAlpha = 1),error = function(e) NA)
    neuron$volume$volume.estimation$primary.dendrite.ashape.volume = primary.dendrite.ashape.volume
    neuron$volume$mesh3d$primary.dendrite = ashape2mesh3d(primary.dendrite.ashape, remove.interior.points = FALSE)
    message("primary dendrite volume created")
  }
  if(TRUE%in%mesh.vertices$microtubule){
    mt.ashape = alphashape3d::ashape3d(as.matrix(subset(mesh.vertices,microtubule==TRUE)[,1:3]),pert = TRUE,alpha = voxelSize*4)
    mt.ashape.volume = tryCatch(alphashape3d::volume_ashape3d(as3d=mt.ashape, byComponents = FALSE, indexAlpha = 1),error = function(e) NA)
    neuron$volume$volume.estimation$microtubule.ashape.volume = mt.ashape.volume
    neuron$volume$mesh3d$microtubule = ashape2mesh3d(mt.ashape, remove.interior.points = FALSE)
    message("microtubule volume created")
  }
  if(FALSE%in%mesh.vertices$microtubule){
    twig.ashape = alphashape3d::ashape3d(as.matrix(subset(mesh.vertices,microtubule==FALSE)[,1:3]),pert = TRUE,alpha = voxelSize*4)
    twig.ashape.volume = tryCatch(alphashape3d::volume_ashape3d(as3d=twig.ashape, byComponents = FALSE, indexAlpha = 1),error = function(e) NA)
    neuron$volume$volume.estimation$twig.ashape.volume = twig.ashape.volume
    neuron$volume$mesh3d$twig = ashape2mesh3d(twig.ashape, remove.interior.points = FALSE)
    message("twig volume created")
  }
  neuron$volume$vertices = mesh.vertices
  class(neuron) = c(class(neuron),"neuronvolume")
  return(neuron)
}

#' Visualise neuronal meshes
#'
#' @description Visualise neuron meshes or point clouds, retrieved and downsampled from the Google brainmaps segmentation by Peter Li.
#' @param neuronvolume an object of class neuronvolume, as returned by catnat::fafb_segs_stitch_volumes
#' @param volumes a list of mesh3d objects, retrieved using fafbseg::read_brainmaps_meshes
#' @param type whether to plot meshes of a point cloud
#' @param split whether to plot the whole neuron or some subset of it
#' @param alpha mesh transparency
#' @param synapse.radius radius information passed to rgl:spheres3d
#' @param cols colours for different parts of the neuron
#' @param WithConnectors whether or not to plot synaptic connectors using spheres3D()
#' @param ... methods passed to Rvcg functions: vcgUniformRemesh, vcgSmooth, vcgClost
#' @export
#' @rdname fafb_segs_stitch_volumes
neuronvolume3d <- function(neuronvolume,
                           type = c("volume","points"),
                           split = c("whole","split", "soma", "primary neurite","axon", "dendrite","microtubule", "primary dendrite"),
                           alpha = 0.3, WithConnectors = TRUE, synapse.radius=500,
                           cols = c(neuron = "grey",
                                    dendrite = "cyan",
                                    soma = "magenta",
                                    primary.neurite = "purple",
                                    primary.dendrite = "chartreuse",
                                    axon = "orange",
                                    microtubule = "green",
                                    twig = "brown")){
  type = match.arg(type)
  split = match.arg(split)
  names(cols) = c("neuron","dendrite","soma","primary.neurite","primary.dendrite","axon","microtubule","twig")[1:length(cols)]
  ### plot mesh3d objects ###
  if(type=="volume"){
    if(split=="whole"){
      rgl::plot3d(neuronvolume$volume$mesh3d$whole, col = cols["neuron"], alpha = alpha, add = TRUE)
      if(WithConnectors){
        rgl::spheres3d(get.synapses(neuronvolume,"POST"),col="cyan",radius = synapse.radius/3)
        rgl::spheres3d(get.synapses(neuronvolume,"PRE", polypre = FALSE),col="red",radius = synapse.radius)
      }
    }
    if (grepl("axon|split",split)){
      rgl::plot3d(neuronvolume$volume$mesh3d$axon, col = cols["axon"], alpha = alpha, add = TRUE)
      if(WithConnectors){
        rgl::spheres3d(get.synapses(axonic_cable(neuronvolume),"POST"),col="cyan",radius = synapse.radius/3)
        rgl::spheres3d(get.synapses(axonic_cable(neuronvolume),"PRE", polypre = FALSE),col="red",radius = synapse.radius)
      }
    }
    if (grepl("dendrite|split",split)){
      rgl::plot3d(neuronvolume$volume$mesh3d$dendrite, col = cols["dendrite"], alpha = alpha, add = TRUE)
      if(WithConnectors){
        rgl::spheres3d(get.synapses(dendritic_cable(neuronvolume),"POST"),col="cyan",radius = synapse.radius/3)
        rgl::spheres3d(get.synapses(dendritic_cable(neuronvolume),"PRE", polypre = FALSE),col="red",radius = synapse.radius)
      }
    }
    if (grepl("primary.dendrite|split",split)){
      rgl::plot3d(neuronvolume$volume$mesh3d$primary.dendrite, col = cols["primary.dendrite"], alpha = alpha, add = TRUE)
    }
    if(grepl("split|soma",split)){
      rgl::plot3d(neuronvolume$volume$mesh3d$soma, col = cols["soma"], alpha = alpha, add = TRUE)
    }
    if(grepl("split|primary.neurite",split)){
      rgl::plot3d(neuronvolume$volume$mesh3d$primary.neurite, col = cols["primary.neurite"], alpha = alpha, add = TRUE)
    }
    if(split=="microtubule"){
      rgl::plot3d(neuronvolume$volume$mesh3d$microtubule, col = cols["microtubule"], alpha = alpha, add = TRUE)
      rgl::plot3d(neuronvolume$volume$mesh3d$twig, col = cols["twig"], alpha = alpha, add = TRUE)
      if(WithConnectors){
        rgl::spheres3d(get.synapses(neuronvolume,"POST"),col="cyan",radius = synapse.radius/3)
        rgl::spheres3d(get.synapses(neuronvolume,"PRE", polypre = FALSE),col="red",radius = synapse.radius)
      }
    }
  }
  ### plot 3D points ###
  if(type=="points"){
    if(split=="whole"){
      rgl::points3d(nat::xyzmatrix(neuronvolume$volume$vertices), col = cols["neuron"], alpha = alpha, add = TRUE)
      if(WithConnectors){
        rgl::spheres3d(get.synapses(neuronvolume,"POST"),col="cyan",radius = synapse.radius/3)
        rgl::spheres3d(get.synapses(neuronvolume,"PRE", polypre = FALSE),col="red",radius = synapse.radius)
      }
    } else if (grepl("axon|split",split)){
      rgl::points3d(nat::xyzmatrix(subset(neuronvolume$volume$vertices, Label==2)), col = cols["axon"], alpha = alpha, add = TRUE)
      if(WithConnectors){
        rgl::spheres3d(get.synapses(axonic_cable(neuronvolume),"POST"),col="cyan",radius = synapse.radius/3)
        rgl::spheres3d(get.synapses(axonic_cable(neuronvolume),"PRE", polypre = FALSE),col="red",radius = synapse.radius)
      }
    } else if (grepl("dendrite|split",split)){
      rgl::points3d(nat::xyzmatrix(subset(neuronvolume$volume$vertices, Label==3)), col = cols["dendrite"], alpha = alpha, add = TRUE)
      if(WithConnectors){
        rgl::spheres3d(get.synapses(dendritic_cable(neuronvolume),"POST"),col="cyan",radius = synapse.radius/3)
        rgl::spheres3d(get.synapses(dendritic_cable(neuronvolume),"PRE", polypre = FALSE),col="red",radius = synapse.radius)
      }
    }
    if(grepl("split|soma",split)){
      rgl::points3d(nat::xyzmatrix(subset(neuronvolume$volume$vertices, Label==1)), col = cols["soma"], alpha = alpha, add = TRUE)
    }
    if(grepl("split|primary.neurite",split)){
      rgl::points3d(nat::xyzmatrix(subset(neuronvolume$volume$vertices, Label==7)), col = cols["primary.neurite"], alpha = alpha, add = TRUE)
    }
    if(grepl("split|primary.dendrite",split)){
      rgl::points3d(nat::xyzmatrix(subset(neuronvolume$volume$vertices, Label==7)), col = cols["primary.dendrite"], alpha = alpha, add = TRUE)
    }
    if(split=="microtubule"){
      rgl::points3d(nat::xyzmatrix(subset(neuronvolume$volume$vertices, microtubule==TRUE)), col = cols["microtubule"], alpha = alpha, add = TRUE)
      rgl::points3d(nat::xyzmatrix(subset(neuronvolume$volume$vertices, microtubule==FALSE)), col = cols["twig"], alpha = alpha, add = TRUE)
      if(WithConnectors){
        rgl::spheres3d(get.synapses(neuronvolume,"POST"),col="cyan",radius = synapse.radius/3)
        rgl::spheres3d(get.synapses(neuronvolume,"PRE", polypre = FALSE),col="red",radius = synapse.radius)
      }
    }
  }
}


#' Update radius information for skeletons in a CATMAID instance
#'
#' @description  Update radius information for skeletons in a CATMAID instance using Peter Li's segmentation of FAFB-v14. Brainmaps API access required
#' @param x the treenode ids to edit
#' @param radii a vector the same length as tnids, giving the new radius for each treenode id in that vector
#' @param max.dist the radius is calculated as the mean distance of the nearest 10 mesh vertices for a 3D volume to a each treenode the mesh encompasses. Max.dist sets the maximum distance fro which to search for the ten closest nodes. If exceeded, the radius is set to max.dist.
#' @param method whether to use raycast (casts 10 rays perpendicular to the path of the neuron for each point, slower but more accurate) or the nearest point on the bounding mesh, to estimate node radius
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_update_radius
fafbseg_update_node_radii <- function(x, max.dist = 2000, method = c("nearest.mesh.point","ray.cast"),
                                      pid = 1, conn = NULL, ...){
  if(!requireNamespace('fafbseg', quietly = TRUE))
    stop("Please install suggested fafbseg package")
  method = match.arg(method)
  if(!is.neuronlist(x)){
    message("Reading neurons from ", catmaid_get_server(conn))
    neurons = catmaid::read.neurons.catmaid(x, OmitFailures = TRUE, pid=pid,conn=conn, ...)
  }else{
    neurons = x
  }
  for(i in 1:length(neurons)){
    message("Finding the 3D segments that correspond to neuron ", i, " of ", length(neurons))
    print(neurons[i,])
    neuron = neurons[i][[1]]
    mapping = map_fafbsegs_to_neuron(neuron, node.match = node.match)
    nglids = unique(mapping$ngl_id)
    nglids = nglids[nglids!=0]
    volumes = fafbseg_get_volumes(nglids)
    message("Calculating radius information")
    radii = neuronvolume_get_radius(neuron, volumes, max.dist = max.dist)
    message("Updating CATMAID")
    catmaid_update_radius(tnids = radii$PointNo,radii = radii$radius,pid=pid,conn=conn, ...)
    message("Radii updated for ", names(neurons[i,]))
  }
}


# Hidden
fafbseg_get_volumes <- function(nglids){
  volumes = list()
  nams = c()
  pb <- utils::txtProgressBar(min = 0, max = length(nglids), style = 3)
  while(length(nglids)>0){
    message("Getting volumes for ",length(nglids), " Neuroglancer segments")
    for(ng in 1:length(nglids)){
      s = fafbseg::find_merged_segments(nglids[ng])
      vol = tryCatch(fafbseg::read_brainmaps_meshes(s), error = function(e) "read error")
      if(is.character(vol)){
        message(vol)
      }else if(is.null(vol)|length(vol)<1){
        message("no volume to read")
        vol = NULL
        nams = c(nams,nglids[ng])
        volumes[[ng]] = vol
      }else{
        message(paste0("volume read: ",ng))
        nams = c(nams,nglids[ng])
        volumes[[ng]] = vol
      }
    }
    nglids = setdiff(nglids,nams)
    message("re-trying: ", length(nglids))
    utils::setTxtProgressBar(pb, ng)
  }
  close(pb)
  volumes
}

# Hidden
neuronvolume_get_radius <- function(neuron, volumes, max.dist = 2000, method = c("nearest.mesh.point","ray.cast")){
  method = match.arg(method)
  message("Estimating node radii")
  neuron.points = nat::xyzmatrix(neuron)
  neuron$d$radius =  NA
  pb <- utils::txtProgressBar(min = 0, max = length(volumes), style = 3)
  for(i in 1:length(volumes)){
    v = volumes[[i]]
    p = tryCatch(nat::pointsinside(neuron.points,surf=v),error = function(e) FALSE)
    dv= nat::xyzmatrix(v)
    if(sum(p)>0){
      if(method=="ray.cast"){
        query = neuron.points[p,]
        if(sum(p)==1){query=matrix(neuron.points[p,],ncol=3)}
        near = dv[nabor::knn(query=query,data=dv,
                             k=1,radius=max.dist)$nn.idx,]
        query.parent = subset(neuron$d,PointNo%in%neuron$d$Parent[p])
        query.parent = nat::xyzmatrix(query.parent[match(neuron$d$Parent[p],query.parent$PointNo),])
        query.vectors = query - query.parent
        near.vectors = near - query.parent
        tangent = do.call(rbind,lapply(1:nrow(query.vectors), function(q)
          Morpho::crossProduct(c(query.vectors[q,]),c(near.vectors[q,]))))
        intersect <- lapply(1:nrow(query.vectors), function(t)
          Morpho::meshPlaneIntersect(v,query[t,],near[t,],tangent[t,]))
        r = c()
        for(j in 1:nrow(query)){
          int = intersect[[j]]
          int = int[sample(1:nrow(int),10),]
          if(length(int)){
            rays = Rvcg::setRays(coords = do.call("rbind", replicate(nrow(int), query[j,], simplify = FALSE)),
                                 dirs = int-query[j,])
            raycast = Rvcg::vcgRaySearch(x=rays, mesh = v, mintol = 0, maxtol = 1e+15, mindist = FALSE,threads = 1)
            r = c(r,median(raycast$distance))
          }else{
            r  = c(r, max.dist)
          }
        }
      }else{
        query = neuron.points[p,]
        clost = Rvcg::vcgClost(x = query,mesh = v)
        r  = c(r, clost$quality)
      }
      neuron$d$radius[p] = r
    }
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  neuron$d
}

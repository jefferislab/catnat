# Increase sigma, increase randomness ofgrowth direction
dendritic.morphogenisis <- function(start = c(173.038, 76.40729, 70.57396), root = c(212.2850, 54.08630, 50.54488), sigma = 0.25, inertia.with.microtubules = 1, inertia.without.microtubules = 0.5, root.tropism = 1, self.repulsion = 5, growth.rate = 1, boundary.repulsion = 10, root.tropism.decay = 1e-2, self.repulsion.decay = 1e-2, growth.rate.decay = 1e-2, boundary.repulsion.decay = 1, microtubule.termination = 0.001, microtubule.branch.probability = 0.1, no.microtubule.branch.probability = 0.2, microtubule.inheritance = 0.3, intersection.proximity = 0.1, mean.bifurcations.without.microtubules = 10, sd.bifurcations.without.microtubules = 5, mean.geodesic.length.without.microtubules = 15, sd.geodesic.length.without.microtubules = 10, cable.length.mean = 1500, cable.length.sd = 100, cable.length.microtubules.mean = 1500, cable.length.microtubules.sd = 100, arbourisation.zone = subset(FCWBNP.surf,"LH_R"), boundary.shape = FCWBNP.surf, resample.boundary.shape = FALSE, visualise = FALSE) {
  # Are we going to see where we are going?
  if(visualise){
    #if(!is.null(boundary.shape)){plot3d(boundary.shape,alpha=0.3,col="pink",add=TRUE)}
    if(!is.null(arbourisation.zone)){plot3d(arbourisation.zone,alpha=0.1,col="red",add=TRUE)}
  }
  # Define some functions that we will be using
  marsaglia.polar.method <- function(sigma){
    w=0
    while (w<=0|w>1) {
      v1 = runif(1, min = -1, max = 1)
      v2 = runif(1, min = -1, max = 1)
      w = v1^2*v2^2
    }
    sigma*v1*sqrt((-2*log(v1^2*v2^2))/(v1^2*v2^2))
  }
  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  norm_vec <- function(x) sqrt(sum(x^2))
  unit.vector <- function(x) {x / sqrt(sum(x^2))}
  # Function to calculate the self repulsive force in the dendrite
  self.repulsion.force <- function(growth.front,df, self.repulsion, self.repulsion.decay){
    points = nat::xyzmatrix(df)
    euc = nabor::knn(data=points,query=matrix(growth.front,ncol=3),k=nrow(points)) # Repulsion away from the closest 10 points?
    euc.distances = euc$nn.dists
    points = as.matrix(points[euc$nn.idx,],ncol=3)
    direction.vectors = t(apply(points,1,function(x) unit.vector(growth.front-x)))
    direction.vectors[is.nan(direction.vectors)]=0
    # Need to also create a model for this...
    self.repulsion.magnitude = cbind(direction.vectors,t(self.repulsion*exp(-euc.distances*self.repulsion.decay))) # Need a better model????
    self.repulsion.forces = t(apply(self.repulsion.magnitude,1,function(x) x[1:3]*x[4]))
    c(colSums(self.repulsion.forces))
  }
  boundary.repulsive.force <- function(growth.front,boundary.shape, boundary.repulsion.decay, boundary.repulsion){
    # Shall we just propel away from the points in the mesh? Points in mesh need to be evenly spaced
    if(is.null(boundary.shape)){return(0)
    }else{
      # We need to resample the mesh to get an even spacing of points
      #points= nat::xyzmatrix(boundary.shape) # These are the evenly spaced mesh points
      # Let's add in the closest point as well, or maybe just the closest point?
      #points = rbind(points,t(Morpho::closemeshKD(matrix(growth.front,ncol=3), mesh = as.mesh3d(boundary.shape))$vb)[,-4])
      points = matrix(Morpho::closemeshKD(matrix(growth.front,ncol=3), mesh = boundary.shape)$vb[-4],ncol=3)
      # Calculate euclidean distance
      euc = nabor::knn(data=points,query=matrix(growth.front,ncol=3),k=nrow(points)) # Repulsion away from the cloest 10 points
      euc.distances = euc$nn.dists
      points = points[euc$nn.idx,]
      # Calculate diection vectors
      direction.vectors = growth.front-points
      direction.vectors[is.nan(direction.vectors)]=0
      # Need to also create a model for this...
      boundary.repulsion.magnitude = as.numeric(boundary.repulsion*exp(-euc.distances*boundary.repulsion.decay)) # Need a better model????
      c(boundary.repulsion.magnitude*unit.vector(direction.vectors))
    }
  }
  # Function to determining ending growth
  mark.end <- function(df, arbourisation.zone=NULL, intersection.proximity = 0.1, mean.bifurcations.without.microtubules = 3, sd.bifurcations.without.microtubules = 1, mean.geodesic.length.without.microtubules = 15, sd.geodesic.length.without.microtubules = 5){
    leaves = nat::endpoints(df)[nat::endpoints(df)%in%subset(df,end==FALSE)$PointNo]
    # End because it is getting too close to another dendritic point - extra: UNLESS the two points share the same parent
    proximity.end = c(nabor::knn(data=nat::xyzmatrix(df),query=nat::xyzmatrix(subset(df,PointNo%in%leaves)),k=2)$nn.dists[,2]<intersection.proximity)
    df[df$PointNo%in%leaves[proximity.end],"end"] = TRUE
    # End because you are exceeding the boundary.
    if(!is.null(arbourisation.zone)){
      boundary.end = !nat::pointsinside(x=nat::xyzmatrix(subset(df,PointNo%in%leaves)),surf=arbourisation.zone)
      df[df$PointNo%in%leaves[boundary.end],"end"] = TRUE
    }
    leaves.mt.loss = leaves[leaves%in%subset(df,microtubules==FALSE)$PointNo]
    if(length(leaves.mt.loss)!=0){
      #ng = nat::as.ngraph(df)
      #paths.to.mt.loss = lapply(leaves.mt.loss, function(x) unlist(igraph::shortest_paths(ng,from=1,to=x)$vpath))
      #path.from.mt = lapply(paths.to.mt.loss, function(x) x[which(!x%in%subset(df,microtubules==TRUE)$PointNo)])
      #branching.end = sapply(path.from.mt, function(x) sum(nat::branchpoints(df)%in%x) > rnorm(1, mean = mean.bifurcations.without.microtubules, sd = sd.bifurcations.without.microtubules))
      #df[df$PointNo%in%leaves.mt.loss[branching.end],"end"] = TRUE
      #length.end = sapply(path.from.mt, function(x) length(x) > rnorm(1, mean = mean.geodesic.length.without.microtubules, sd = sd.geodesic.length.without.microtubules))
      #df[df$PointNo%in%leaves.mt.loss[length.end],"end"] = TRUE
      no.branch.points.mt.loss = subset(df,PointNo%in%leaves.mt.loss)$branch.points.since.microtubules
      branching.end = sapply(no.branch.points.mt.loss,function(x) x > rnorm(1, mean = mean.bifurcations.without.microtubules, sd = sd.bifurcations.without.microtubules))
      df[df$PointNo%in%leaves.mt.loss[branching.end],"end"] = TRUE
      geodesic.distance.mt.loss = subset(df,PointNo%in%leaves.mt.loss)$distance.from.microtubules
      length.end = sapply(geodesic.distance.mt.loss, function(x) x > rnorm(1, mean = mean.geodesic.length.without.microtubules, sd = sd.geodesic.length.without.microtubules))
      df[df$PointNo%in%leaves.mt.loss[length.end],"end"] = TRUE
    }
    df
  }
  if(resample.boundary.shape){
    boundary.shape = nat::as.hxsurf(Rvcg::vcgUniformRemesh(as.mesh3d(boundary.shape), mergeClost= TRUE, multiSample=TRUE, silent=TRUE, offset=2, absDist=TRUE)) # Resampled boundary shape evenly,make it slightly larger (2um) so that processes can actually reach the boundary
  }
  if(!is.null(boundary.shape)){boundary.shape = as.mesh3d(boundary.shape)}
  # This generates an inner and an outer object, we need to get rid of the inner one
  cable.length = 0
  max.cable.length = rnorm(n = 1, mean = cable.length.mean, sd = cable.length.sd)
  max.cable.microtubules.length = rnorm(n = 1, mean = cable.length.microtubules.mean, sd = cable.length.microtubules.sd)
  max.cable.length = ifelse(max.cable.length<max.cable.microtubules.length,max.cable.length,max.cable.microtubules.length)
  # Get our starting neuron
  df = data.frame(PointNo = 1, label = 0, X = root[1], Y= root[2], Z = root[3], W = -2, Parent = -1,end=TRUE, microtubules = TRUE, distance.from.start = -1, distance.from.microtubules = 0, branch.points = 0, branch.points.since.microtubules = 0, stepsize = 0)
  df = rbind(df,data.frame(PointNo = 2, label = 0, X = start[1], Y = start[2], Z = start[3], W = -2, Parent = 1,end=FALSE, microtubules = TRUE, distance.from.start = 0, distance.from.microtubules = 0, branch.points = 0, branch.points.since.microtubules = 0, stepsize = 0))
  leaves = nat::endpoints(df)[nat::endpoints(df)%in%subset(df,end==FALSE)$PointNo]
  while(length(leaves)>0&cable.length<max.cable.length){
    for(l in leaves){
      # Let;s get our growth position
      growth.front = c(xyzmatrix(df[l,]))
      # A non-unifrom growth rate account better to tortuosity. We want the rate to decrease as we gain distance/branchpoints from the soma.         Let's determine our growth rate, which will decrease exponentially as a function of branch points generated
      #ng = as.ngraph(df)
      #path.from.start = unlist(igraph::shortest_paths(ng,from=1,to=l)$vpath)
      #no.branch.points = sum(nat::branchpoints(df)%in%path.from.start)
      #no.branch.points = ifelse(no.branch.points==0,1,no.branch.points)
      #geodesic.distance = length(path.from.start)
      no.branch.points = df[l,]$branch.points
      geodesic.distance = df[l,]$distance.from.start
      rate.of.growth = growth.rate*exp(-geodesic.distance*growth.rate.decay)
      rate.of.growth = growth.rate*exp(-geodesic.distance*growth.rate.decay)
      # If microtubular, is this going to be a terminal branch?
      if(no.branch.points>0){terminal = microtubule.termination*exp(-1/no.branch.points) > runif(1, min = 0, max = 1)}else{terminal=FALSE}
      #Continue with current segment
      if(df[l,]$microtubules){
        root.tropic.force = root.tropism*exp(-(geodesic.distance+1)*root.tropism.decay)*unit.vector(c(growth.front-root))
        inertial.force = inertia.without.microtubules*unit.vector(c(t(growth.front-subset(df,PointNo==df[l,]$Parent)[,c("X","Y","Z")])))
      }else{
        root.tropic.force=0
        inertial.force = inertia.without.microtubules*unit.vector(c(t(growth.front-subset(df,PointNo==df[l,]$Parent)[,c("X","Y","Z")])))
      }
      repulsive.force = self.repulsion.force(self.repulsion=self.repulsion,df=df,growth.front=growth.front,self.repulsion.decay=self.repulsion.decay)
      boundary.force = boundary.repulsive.force(growth.front=growth.front,boundary.shape=boundary.shape, boundary.repulsion.decay=boundary.repulsion.decay, boundary.repulsion=boundary.repulsion)
      growth.force = unit.vector(c(marsaglia.polar.method(sigma=sigma),marsaglia.polar.method(sigma=sigma),marsaglia.polar.method(sigma=sigma)))
      growth.front = growth.front + rate.of.growth*unit.vector(growth.force + inertial.force + root.tropic.force + boundary.force)
      # Are there microtubules in this continuing brtanch?
      mt = ifelse(terminal,FALSE,df[l,]$microtubules)
      #Is this segment going to branch before continuing?
      if(!df[l,]$microtubules){branch.probability=no.microtubule.branch.probability}else{branch.probability=microtubule.branch.probability}
      branch = runif(1, min = 0, max = 1)<branch.probability|terminal # Need to set a model for branching chance, given a.microtubule presence and b.distance from start
      # Add this new data to the data.frame
      df = rbind(df,data.frame(PointNo = max(df$PointNo)+1, label = 0, X = growth.front[1], Y= growth.front[2], Z = growth.front[3], W = -2,
                               Parent = l, end = FALSE, microtubules = mt, distance.from.start = df[l,]$distance.from.start + rate.of.growth, distance.from.microtubules = ifelse(!mt,df[l,]$distance.from.microtubules + rate.of.growth,0), branch.points = df[l,]$branch.points + branch, branch.points.since.microtubules = ifelse(terminal|!mt,df[l,]$branch.points.since.microtubules + branch,0), stepsize = rate.of.growth))
      # Have we exceeded the mamimum cable length?
      if(sum(df$stepsize)>max.cable.length){
        if(visualise){message("Maximum cable length reached")}
        break
      }
      # let's branch?
      if(!df[l,]$microtubules){branch.probability=no.microtubule.branch.probability}else{branch.probability=microtubule.branch.probability}
      if(branch){ # Need to set a model for branching chance, given a.microtubule presence and b.distance from start
        #Will our new branch contain a microtubule?
        if(df[l,]$microtubules){
          # Need to fit a model to the microtubule branching (given distance from start) chance in the EM data to better fromulate this
          microtubule.branch.chance = microtubule.inheritance*exp(-no.branch.points/10)
          mt = runif(1, min = 0, max = 1)<microtubule.branch.chance
        }else{mt=FALSE}
        growth.front = c(xyzmatrix(df[l,]))
        rate.of.growth = growth.rate*exp(-(geodesic.distance+1)*growth.rate.decay)
        # Inertia different with an without microtubule?
        repulsive.force = self.repulsion.force(self.repulsion=self.repulsion,df=df,growth.front=growth.front,self.repulsion.decay=self.repulsion.decay)
        boundary.force = boundary.repulsive.force(growth.front=growth.front,boundary.shape=boundary.shape, boundary.repulsion.decay=boundary.repulsion.decay, boundary.repulsion=boundary.repulsion)
        growth.force = unit.vector(c(marsaglia.polar.method(sigma=sigma),marsaglia.polar.method(sigma=sigma),marsaglia.polar.method(sigma=sigma)))
        growth.front = growth.front + rate.of.growth*unit.vector(growth.force + inertial.force + root.tropic.force + boundary.force)
        df = rbind(df,data.frame(PointNo = max(df$PointNo)+1, label = 0, X = growth.front[1], Y = growth.front[2], Z =
                                   growth.front[3], W = -2,Parent = l, end = FALSE, microtubules = mt, distance.from.start = df[l,]$distance.from.start + rate.of.growth, distance.from.microtubules = ifelse(!mt,df[l,]$distance.from.microtubules + rate.of.growth,0), branch.points = df[l,]$branch.points + branch, branch.points.since.microtubules = ifelse(terminal|!mt,df[l,]$branch.points.since.microtubules + branch,0),stepsize = rate.of.growth))
      }
    }
    rownames(df) = 1:nrow(df)
    leaves = nat::endpoints(df)[nat::endpoints(df)%in%subset(df,end==FALSE)$PointNo]
    df = mark.end(df, arbourisation.zone=arbourisation.zone, intersection.proximity = intersection.proximity, mean.bifurcations.without.microtubules = mean.bifurcations.without.microtubules, sd.bifurcations.without.microtubules = sd.bifurcations.without.microtubules, mean.geodesic.length.without.microtubules = mean.geodesic.length.without.microtubules, sd.geodesic.length.without.microtubules = sd.geodesic.length.without.microtubules) # One round of looking at growing ends complete
    # Have we exceeded the mamimum cable length?
    if(sum(df$stepsize)>max.cable.length){
      if(visualise){message("Maximum cable length reached")}
      break
    }
    if(visualise){
      neuron = nat::as.neuron(df)
      #leaves.parents = subset(df,PointNo%in%leaves)$Parent
      #nn = nat::prune_vertices(neuron, c(leaves.parents,leaves,1,2,3), invert = TRUE)
      #plot3d(nn,col=rand_color(),WithNodes = FALSE);Sys.sleep(sum(visualise)-1)
      df.mt.loss = subset(df,microtubules==FALSE)
      leaves.mt.loss = leaves[leaves%in%df.mt.loss$PointNo]
      leaves.mt = leaves[!leaves%in%leaves.mt.loss]
      leaves.mt.loss.parents = subset(df,PointNo%in%leaves.mt.loss)$Parent
      leaves.mt.parents = subset(df,PointNo%in%leaves.mt)$Parent
      if(length(leaves.mt)>0){
        mt.cols <- colorRampPalette(c("darkred", "lightpink"))
        mt.col = mt.cols(200)[ifelse(max(leaves.mt)>200,200,max(leaves.mt))]
        nn.mt = nat::prune_vertices(neuron, c(leaves.mt.parents,leaves.mt,1,2,3), invert = TRUE)
        plot3d(nn.mt,col=mt.col,WithNodes = FALSE);Sys.sleep(sum(visualise)-1)
        rgl::points3d(nat::xyzmatrix(subset(df,end==TRUE&PointNo%in%leaves.mt)),col="red")
      }
      if(length(leaves.mt.loss)>0){
        mt.loss.cols <- colorRampPalette(c("chartreuse4", "palegreen1"))
        col.no = max(subset(df,PointNo%in%leaves.mt.loss)$distance.from.microtubules)
        mt.loss.col = mt.loss.cols(mean.geodesic.length.without.microtubules)[ifelse(max(leaves.mt.loss)>mean.geodesic.length.without.microtubules,mean.geodesic.length.without.microtubules,col.no)]
        nn.mt.loss = nat::prune_vertices(neuron, c(leaves.mt.loss.parents,leaves.mt.loss,1,2,3), invert = TRUE)
        plot3d(nn.mt.loss,col=mt.loss.col,WithNodes = FALSE);Sys.sleep(sum(visualise)-1)
        rgl::points3d(nat::xyzmatrix(subset(df,end==TRUE&PointNo%in%leaves.mt.loss)),col="cyan")
      }
      cable.length = sum(df$stepsize)
      # We want to just grow the microtubular portions first, and then the non-microtubular poritons (?)
      if(max.cable.microtubules.length>sum(subset(df,microtubules= TRUE)$stepsize)){
        leaves = nat::endpoints(subset(df,microtubules = TRUE))[nat::endpoints(df)%in%subset(df,end==FALSE)$PointNo]
      }else{
        leaves = nat::endpoints(df)[nat::endpoints(df)%in%subset(df,end==FALSE)$PointNo]
      }
    }
  }
  nat::as.neuron(df)
}


dendritic.morphogenisis.neuronlist<-function(n=10, names = paste0("simulatd-dendrite-",1:n), start = c(173.038, 76.40729, 70.57396), root = c(212.2850, 54.08630, 50.54488), sigma = 0.25, inertia.with.microtubules = 1, inertia.without.microtubules = 0.5, root.tropism = 1, self.repulsion = 5, growth.rate = 1, boundary.repulsion = 1, root.tropism.decay = 1e-2, self.repulsion.decay = 1e-2, growth.rate.decay = 1e-2, boundary.repulsion.decay = 100, max.cable = 100, branch.probability = 0.1, microtubule.branch.probability = 0.1, intersection.proximity = 0.1, average.bifurcations.without.microtubules = 5, sd.bifurcations.without.microtubules = 4, average.geodesic.length.without.microtubules = 20, sd.geodesic.length.without.microtubules = 5, max.cable.length = 1500, arbourisation.zone = subset(FCWBNP.surf,"LH_R"), boundary.shape = FCWBNP.surf, resample.boundary.shape = FALSE, visualise = FALSE){
  nl = nlapply(names, function(x) dendritic.morphogenisis(start = start, root = root, sigma = sigma, inertia.with.microtubules = inertia.with.microtubules, inertia.without.microtubules = inertia.without.microtubules, root.tropism = root.tropism, self.repulsion = self.repulsion, growth.rate = growth.rate, boundary.repulsion = boundary.repulsion, root.tropism.decay = root.tropism.decay, self.repulsion.decay = self.repulsion.decay, growth.rate.decay = growth.rate.decay, boundary.repulsion.decay = boundary.repulsion.decay, max.cable = max.cable, branch.probability = branch.probability, microtubule.branch.probability = microtubule.branch.probability, intersection.proximity = intersection.proximity, average.bifurcations.without.microtubules = average.bifurcations.without.microtubules, sd.bifurcations.without.microtubules = sd.bifurcations.without.microtubules, average.geodesic.length.without.microtubules = average.geodesic.length.without.microtubules, sd.geodesic.length.without.microtubules = sd.geodesic.length.without.microtubules, max.cable.length = max.cable.length, arbourisation.zone = arbourisation.zone, boundary.shape = boundary.shape, resample.boundary.shape = resample.boundary.shape, visualise = visualise) )
  attr(nl,"df") = match.call()[-1:4]
  nl
}







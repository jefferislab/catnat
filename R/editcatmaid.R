#' Upload neuron(s) to CATMAID
#'
#' @description  Uploads neurons to CATMAID, names them and annotates them.
#' Please use with caution, as you could be heavily adding to a live tracing environment.
#' @param swc a neuron or neuronlist object, or (a) local file path(s) to your saved .swc file(s).
#' @param name whatever you want to name your uploaded neurons. If a single character, then it will be added to all uploaded neurons. Else, can be a character vector the same length as swc.
#' @param annotations a character vector of annotations, to be added to all of the uploaded neurons
#' @param return.new.skids if TRUE, the new skids created in the CATMAID instance specificed by conn, are returned
#' @param include.tags whether of not, if swc is a CATMAID neuron/neuronlist, to transfer its tags to its newly uploaded cognate
#' @param include.connectors whether of not, if swc is a CATMAID neuron/neuronlist, to transfer its connectors to its newly uploaded cognate
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param max.upload the maximum number of files that the function will allow you to upload at once
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_upload_neurons
catmaid_upload_neurons <- function (swc, name ="neuron SWC upload", annotations = "SWC upload",
                                    include.tags = TRUE,include.connectors = TRUE,
                                    search.range.nm = 100, return.new.skids = FALSE,
                                    pid = 1, conn = NULL, max.upload = 10, ...) {
  if(is.neuronlist(swc)|is.neuron(swc)){
    temp.files = paste0(names(swc),".swc")
    nat::write.neurons(swc, format= "swc",dir = tempdir(), files = temp.files, ...)
    swc = paste0(tempdir(),temp.files)
  }
  if (length(swc)>max.upload){
    stop(paste0('You are uploading a large number of SWC files.
                Are you sure you want to do this?
                If so change the value of the max.upload argument'))
  }
  if(length(name)==1&length(swc)>1){
    name = paste0(name,"_",seq_along(swc))
  } else if (length(name)>length(swc)|length(swc)>length(name)){
    stop(paste0('The names argument mut be a vector or length 1, to be applied to all neurons
                or a vector the same length as swc'))
  }
  skids = c(),
  for(file in 1:length(swc)){
    post_data = list()
    post_data["name[1]"] = as.list(name[file])
    post_data["file[1]"] = list(upload_file(swc[file]))
    path = sprintf("/%d/skeletons/import", pid)
    res = catmaid::catmaid_fetch(path, body = post_data, include_headers = F,
                                 simplifyVector = T, conn = conn, ...)
    invisible(catmaid:::catmaid_error_check(res))
    if(include.tags&is.neuron(swc[[file]])){
      fafbseg_transfer_tags(new.neuron = res$skeleton_id, old.neuron = swc[[file]], search.range.nm = search.range.nm, pid = pid, conn = conn, ...)
    }
    if(include.connectors&is.neuron(swc[[file]])){
      fafbseg_transfer_connectors(new.neuron = res$skeleton_id, old.neuron = swc[[file]], links = TRUE, search.range.nm = search.range.nm, pid = pid, conn = conn, ...)
    }
    skids = c(skids,res$skeleton_id)
  }
  if(!is.null(annotations)){
    sapply(seq_along(skids), function(x) catmaid_rename_neuron(skids = skids[x], names = name[x], pid = pid, conn = conn))
    catmaid::catmaid_set_annotations_for_skeletons(skids = as.numeric(skids),
                                                   annotations = annotations,
                                                   pid = pid, conn = conn, ...)
    message(paste0("Annotations ", annotations," set to new skeletons: ", res$skeleton_id))
  }
  if(return.new.skids){
    return(skids)
  }
}

# fafbseg.skid = 210077691
# skid = 2428780

#' Transfer information between v14 and v14-seg neurons
#'
#' @description  Functions to enable the mapping of connectors and tags between cognate FAFB v14 and v14-seg neurons, to help enable neuronal reconstruction using the partial auto-segmentation in v14-seg space
#' @param new.neuron a neuron object, or the skeleton ID of a neuron in v14 space
#' @param fafbseg.skid a neuron object, or the skeleton ID of a neuron in v14-seg space
#' @param search.range.nm the maimum distance between points in the v14 and v14-seg skeleton, that is acceptable for the transfer ot tag/connector information. Not set to 0 be default to accept small variations in node position that might occur between the two tracing environments, due to tracer edit.s
#' @param links whether or not to link transferred connectors  to post-pre synaptic sites on the neuron identified with the new.neuron argument
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param max.upload the maximum number of files that the function will allow you to upload at once
#' @param ... methods passed to catmaid::catmaid_set_labels
#' @export
#' @rdname fafbseg_transfer
fafbseg_transfer_tags<- function(new.neuron, old.neuron, search.range.nm = 100, pid = 1, conn = NULL, ...){
  if(!is.neuron(new.neuron)){
    new.neuron = catmaid::read.neuron.catmaid(new.neuron)
  }
  if(!is.neuron(old.neuron)){
    old.neuron = read.neurons.fafbseg(old.neuron)[[1]]
  }
  labels = unlist(sapply(1:length(old.neuron$tags),function(t) rep(names(old.neuron$tags[t]),length(old.neuron$tags[[t]]))))
  tagged.points = unlist(old.neuron$tags)
  tagged.points = tagged.points[tagged.points%in%old.neuron$d$PointNo]
  tagged.points = subset(old.neuron$d,PointNo%in%tagged.points)
  tagged.points$tag = labels
  near = nabor::knn(query = tagged.points[,c('X','Y', 'Z')], data =nat::xyzmatrix(new.neuron), k=1, radius=search.range.nm)$nn.idx
  near = near[near!=0]
  labels = tagged.points$tag[near!=0]
  nodes = new.neuron$d$PointNo[near]
  if(length(nodes)!=length(labels)){
    stop("error with assigning tag information")
  }
  if(length(nodes)>0){
    for(i in 1:length(nodes)){
      catmaid::catmaid_set_labels(node = nodes[i] , labels = labels[i], type= "treenode", pid = pid, conn = conn, ...)
    }
  }else{
    NULL
  }
}
#' @export
#' @rdname fafbseg_transfer
fafbseg_transfer_connectors<- function(new.neuron, old.neuron, search.range.nm = 100, links = TRUE, pid = 1, conn = NULL, ...){
  if(!is.neuron(new.neuron)){
    new.neuron = catmaid::read.neuron.catmaid(new.neuron)
  }
  if(!is.neuron(old.neuron)){
    old.neuron = read.neurons.fafbseg(old.neuron)[[1]]
  }
  c.df = old.neuron$connector
  near = nabor::knn(query = nat::xyzmatrix(c.df), data =nat::xyzmatrix(new.neuron), k=1, radius=search.range.nm)$nn.idx
  new.c.df = c.df[near!=0,]
  new.c.df$near_id = new.neuron$d[near,"PointNo"]
  if(nrow(new.c.df)>0){
    for(i in 1:nrow(new.c.df)){
      post_data = list()
      post_data["x"] = new.c.df[i,"x"]
      post_data["y"] = new.c.df[i,"y"]
      post_data["z"] = new.c.df[i,"z"]
      path = sprintf("/%d/connector/create", pid)
      res = catmaid::catmaid_fetch(path, body = post_data, include_headers = F,
                                   simplifyVector = T, conn = conn, ...)
      invisible(catmaid:::catmaid_error_check(res))
      if(links){
        from_id = ifelse(new.c.df[i,"prepost"]==0,new.c.df[i,"near_id"],res$connector_id)
        to_id = ifelse(new.c.df[i,"prepost"]!=0,new.c.df[i,"near_id"],res$connector_id)
        link_type = ifelse(new.c.df[i,"prepost"]==0,"presynaptic_to","postsynaptic_to")
        state = sprintf('{"%s":"%s"}, {"%s":"%s"}',from_id, res$connector_edition_time, to_id, res$connector_edition_time)
        post_data = list()
        post_data["from_id"] = from_id
        post_data["to_id"] = to_id
        post_data["link_type"] = link_type
        post_data["state"] = state
        path = sprintf("/%d/link/create", pid)
        res2 = catmaid::catmaid_fetch(path, body = post_data, include_headers = F,
                                     simplifyVector = T, conn = conn, ...)
        invisible(catmaid:::catmaid_error_check(res2))
      }
    }
  }
}

#' Find the location of specified tags for a CATMAID neuron
#'
#' @description  Find the location of tags in a CATMAID neuron, either as URLs to the location of a TODO tag in CATMAID or as a data.frame reporting the location and skeleton treenode locations of specified tags.
#' @param x a neuron or neuronlist object
#' @param tag a single character specifying which tag to look for. Defaults to TODO.
#' @param url if TRUE (defualt) a list of URLs pertaining to specified tag locations are returned. If FALSE, a data.frame subsetted from x$d is returned, reporting treenode ID and X,Y,Z positions for specified tags
#' @param server the CATMAID server for which URLs should be created if url == TRUE
#' @export
#' @rdname catmaid_get_tag
catmaid_get_tag<-function(x, tag = "TODO", url = TRUE, server = "https://neuropil.janelia.org/tracing/fafb/v14/") UseMethod("catmaid_get_tag")
catmaid_get_tag.neuron <- function(x, tag = "TODO", url = TRUE, server = "https://neuropil.janelia.org/tracing/fafb/v14/"){
  TODO = unique(unlist(x$tags[[tag]]))
  if(is.null(TODO)){
    NULL
  }else{
    df = subset(x$d,PointNo%in%TODO)
    if(url){
      catmaid_url = paste0(server, "?pid=1")
      catmaid_url = paste0(catmaid_url, "&zp=", df[["Z"]])
      catmaid_url = paste0(catmaid_url, "&yp=", df[["Y"]])
      catmaid_url = paste0(catmaid_url, "&xp=", df[["X"]])
      catmaid_url = paste0(catmaid_url, "&tool=tracingtool")
      catmaid_url = paste0(catmaid_url, "&sid0=5&s0=0")
      invisible(catmaid_url)
    }
    else{
      df
    }
  }
}
catmaid_get_tag.neuronlist <- function(x, tag = "TODO", url = TRUE, server = "https://neuropil.janelia.org/tracing/fafb/v14/"){
  if(url){
    unlist(lapply(x,catmaid_get_tag.neuron, url=url, tag= tag, server = server))
  }else{
    do.call(rbind,lapply(x,catmaid_get_tag.neuron, url=url, tag = tag, server = server))
  }
}

#' Find the TODO tagged merge sites between a given neuron and the CATMAID database at large
#'
#' @description  Find the location of tags in a CATMAID neuron, either as URLs to the location of a TODO tag in CATMAID or as a data.frame reporting the location and skeleton treenode locations of specified tags.
#' @param TODO a data.frame, as returned by catmaid_get_tag with tag = "TODO" (or a different tag indicative of a merge point) and url = FALSE
#' @param skid the skeleton ID of a neuron for which merges are to be found
#' @param fafbseg whether or not to use fafbseg::read_brainmaps_meshes on the TODO tag locations and restrict possible merges to neurons within the search range and within the cognate auto-segmented volume
#' @param search.range.nm the maimum distance from which to search from the TODO point tag to find potential mergers
#' @param min_nodes the minimum number of nodes a potential merger skeleton needs to have
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @export
#' @rdname catmaid_find_likely_merge
catmaid_find_likely_merge <- function(TODO, skid, fafbseg = FALSE, min_nodes = 2, search.range.nm = 100, pid=1, conn = conn, ...){
  possible.merges = data.frame()
  for(i in 1:nrow(TODO)){
    todo = TODO[i,]
    seg = tryCatch(fafbseg::brainmaps_xyz2id(todo[,c('X','Y', 'Z')]),error=function(e)NULL)
    if(is.null(seg)){
      warning("cannot find auto-segmented volume containing merge tag")
    }else{
      bbx = rbind(todo[,c('X','Y', 'Z')]-search.range.nm,todo[,c('X','Y', 'Z')]+search.range.nm)
      skids.bbx = catmaid_skeletons_in_bbx(bbx=bbx,pid=pid,min_nodes=min_nodes)
      neuron.bbx = catmaid::read.neurons.catmaid(skids.bbx, OmitFailures = TRUE)
      neuron.bbx = neuron.bbx[setdiff(names(neuron.bbx),skid)]
      if(fafbseg){
        require(fafbseg)
        s = fafbseg::find_merged_segments(seg)
        vol = tryCatch(fafbseg::read_brainmaps_meshes(s), error = function(e) "read error")
        in.vol = sapply(neuron.bbx,function(n) sum(nat::pointsinside(nat::xyzmatrix(n),surf=vol))>0)
        neuron.bbx = neuron.bbx[in.vol]
        message(length(neuron.bbx), " fragments found in fafb auto-traced segment corresponding with merge tag")
      }else{
        message(length(neuron.bbx), " fragments found within ", search.range.nm," nm of merge tag")
      }
      neuron.bbx.d = do.call(rbind,lapply(1:length(neuron.bbx),function(n) cbind(neuron.bbx[[n]]$d,skid=names(neuron.bbx[n]))))
      near = nabor::knn(query = todo[,c('X','Y', 'Z')], data =nat::xyzmatrix(neuron.bbx.d),k=10, radius=search.range.nm)$nn.idx
      near = near[near!=0]
      near = neuron.bbx.d[near,]
      near$TODO = todo$PointNo
      near$merger = skid
      possible.merges = rbind(possible.merges,near)
    }
  }
  possible.merges
}

#' Interactively choose to join neurons in CATMAID, visualising with rgl in R
#'
#' @description Interactively choose to join neurons in CATMAID via rgl in R. Use with caution in live tracing environments.
#' @param possible.merges a data frame of possible mergers between tree nodes for neurons in a CATMAID database, as returned via catmaid_find_likely_merge
#' @param brain the brain to plot while visualising potential mergers using rgl. Defaults to NULL, no brain plotted.
#' @param search.range.nm the maimum distance from which to search from the TODO point tag to find potential mergers
#' @param min_nodes the minimum number of nodes a potential merger skeleton needs to have
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @export
#' @rdname catmaid_interactive_join
fafbseg_interactive_join <- function(possible.merges, brain = NULL, pid = 1, conn = conn, ...){
  neurons = read.neurons.fafbseg(skids=unique(possible.merges$merger))
  for(todo in unique(possible.merges$TODO)){
    TODO.possible = subset(possible.merges,TODO==todo)
    continue = FALSE
    i = 1
    nat::no
    while(!progress){
      rgl::clear3d()
      if(!is.null(brain)){
        rgl::plot3d(brain,alpha=0.1,col="grey")
      }
      skid = possible.merges[i,"skid"]
      neuron = catmaid::read.neurons.catmaid(skid, pid = pid, conn = conn, ...)
      merger = catmaid::catmaid_get_treenodes_detail(tnids = possible.merges[i,"TODO"], pid = pid, conn = conn, ...)
      merger.neuron = neurons[as.character(merger$skid)]
      rgl::plot3d(neuron,col="black",lwd=2,soma=T)
      rgl::plot3d(merger.neuron,col="red",lwd=2,soma=T)
      rgl::spheres3d(nat::xyzmatrix(merger), col="orange", radius = 10)
      rgl::spheres3d(nat::xyzmatrix(possible.merges[i,]), col="grey", radius = 10)
      progress = readline("Make merge? y = yes, n = no, c = cycle : ")
      if(progress=="y"){
        nat::npop3d()
        rgl::plot3d(merger.neuron,col="green",lwd=2,soma=T)
        sure = readline("Sure? y = yes, n = no  :")
        if(sure=="y"){
          catmaid_join_skeletons(from_treenode_id = possible.merges[i,"merger"], to_treenode_id = possible.merges[i,"TODO"], pid = pid, conn = conn, ...)
          continue=TRUE
        }else{
          i = ifelse(i==nrow(TODO.possible),i,i+1)
        }
      }else if(progress=="n"){
        TODO.possible = TODO.possible[-i,]
        i = ifelse(i==nrow(TODO.possible),i,i+1)
      }else{
        i = ifelse(i==nrow(TODO.possible),i,i+1)
      }
    }
  }
}


# to_treenode_id = 1254672855
# from_treenode_id = 676710716

#' Programmatically join CATMAID skeletons
#'
#' @description  Programmatically join CATMAID skeletons. Use with caution in live tracing environments.
#' @param from_treenode_id the treenode ID of the downstream neuron
#' @param from_treenode_id the treenode ID of the upstream neuron
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_join_skeletons
catmaid_join_skeletons <- function(from_treenode_id, to_treenode_id, pid = 1, conn = NULL, ...){
  post_data = list()
  post_data["from_id"] = to_treenode_id
  post_data["to_id"] = from_treenode_id
  path = sprintf("/%d/skeleton/join", pid)
  res = catmaid::catmaid_fetch(path, body = post_data, include_headers = F,
                               simplifyVector = T, conn = conn, ...)
  invisible(catmaid:::catmaid_error_check(res))
  message(res$message)
}

#' Search for CATMAID skeletons within a volume
#'
#' @description  Programmatically search for skeleton IDs pertaining to neurons within a search volume defined by a bounding box.
#' @param bbx the bounding box (a matrix of 2 rows and 3 columns) describing a search volume
#' @param min_nodes the minimum number of nodes a neuron in the search area must have (includes ndoes outside search area)
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_skeletons_in_bbx
catmaid_skeletons_in_bbx <- function(bbx, min_nodes = 2, pid = 1, conn = NULL, ...){
  post_data = list()
  post_data["minx"] = bbx[1,1]
  post_data["miny"] = bbx[1,2]
  post_data["minz"] = bbx[1,3]
  post_data["maxx"] = bbx[2,1]
  post_data["maxy"] = bbx[2,2]
  post_data["maxz"] = bbx[2,3]
  post_data["min_nodes"] = min_nodes
  path = sprintf("/%d/skeletons/in-bounding-box", pid)
  res = catmaid::catmaid_fetch(path, body = post_data, include_headers = F,
                               simplifyVector = T, conn = conn, ...)
  res
}

#' Try to detect whether a neuron your are about to import into CATMAID might cause a duplication in that CATMAID instance
#'
#' @description  Programmatically search for skeleton IDs pertaining to neurons within a search volume defined by a bounding box.
#' @param neuron the neuron object you are thinking about uploading
#' @param tolerance how many potentially duplictated nodes you will tolerate
#' @param search.range.nm determines the size of the bounding box around each node in neuron to search for a duplicated. Defaults to 1 nm
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_skeletons_in_bbx
catmaid_duplicated <- function(neuron, tolerance = 0.5, search.range.nm = 100, pid = 1, conn = NULL){
  duplicated = c()
  pb <- utils::txtProgressBar(min = 0, max = nrow(neuron$d), style = 3)
  for(i in 1:nrow(neuron$d)){
    bbx = rbind(nat::xymatrix(neuron$d[i,])-search.range.nm, nat::xymatrix(neuron$d[i,])+search.range.nm)
    skids.bbx = catmaid_skeletons_in_bbx(bbx=bbx,pid=pid, min_nodes = 2, pid = pid, conn = conn)
    duplicated = c(duplicated,length(skids.bbx)>0)
    utils::setTxtProgressBar(pb, n)
  }
  close(pb)
  if(sum(duplicated)/nrow(neuron$d)>tolerance){
    TRUE
  }else{
    FALSE
  }
}

#' Function used to delete single neurons from CATMAID based on a given skeleton ID
#'
#' @description  Delete a single neuron from a given CATMAID instance.
#' Deletes a neuron if and only if two things are the case: 1. The user
#' owns all treenodes of the skeleton modeling the neuron in question and
#' 2. The neuron is not annotated by other users.
#' Use with extreme caution as you may be significantly affecting others' work.
#' You will first be shown the neuron you want to delete and who has worked on it, and then asked whether or not you want to continue.
#' @param skid the skeleton ID for a single neuron you want to delete
#' @param plot whether or not the plot the neuron you are considering deleting, before deciding to end it
#' @param brain the brain you want to plot alongside the neuron you are considering deleting.
#' @param max.nodes the maximum number of nodes a neuron can have, and still be successfully deleted. Helps prevent the accidental deletion of large neurons.
#' @param server the CATMAID server from which to delete the neuron identified by skid
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_delete_neuron
catmaid_delete_neuron <- function(skid,
                                  plot = TRUE, brain = NULL,
                                  max.nodes = 100,
                                  server = "https://neuropil.janelia.org/tracing/fafb/v14/",
                                  pid = 1, conn = NULL, ...){
  if(length(skid)>1){
    stop("length(skid)>1 - You may only attempt to delete one neuron at a time.")
  }
  message("You are about to delete a NEURON from a CATMAID instance.
          This may affect people's hard work.
          You will first be shown the neuron you want to delete and who has worked on it,
          and then asked whether or not you want to continue.")
  progress = readline(paste0("Do you want to continue with: ", skid," in CATMAID instance ", server, " y=yes, n=no : "))
  if(progress=="y"){
    nid = catmaid_get_neuronid(skids = skid, pid = pid, conn = conn, ...)
    neuron = catmaid::read.neurons.catmaid(skid, pid = pid, conn = conn, ...)
    xyz = nat::xyzmatrix(neuron)
    meta = summary(neuron)
    if(meta$nodes>max.nodes){
      stop("The neuron flagged for deletion has ", meta$nodes,
           " which is greater than the given argument max.nodes, ", max.nodes, ". Aborting.")
    }
    catmaid_url = paste0(server, "?pid=1")
    catmaid_url = paste0(catmaid_url, "&zp=", xyz[1,"Z"])
    catmaid_url = paste0(catmaid_url, "&yp=", xyz[1,"Y"])
    catmaid_url = paste0(catmaid_url, "&xp=", xyz[1,"X"])
    catmaid_url = paste0(catmaid_url, "&tool=tracingtool")
    catmaid_url = paste0(catmaid_url, "&active_skeleton_id=", skid)
    catmaid_url = paste0(catmaid_url, "&sid0=5&s0=0")
    ### Plot ###
    if(plot){
      rgl::clear3d()
      if(!is.null(brain)){
        rgl::plot3d(brain,alpha=0.1,col="grey")
      }
      rgl::plot3d(neuron,lwd=3,soma=T,WithConnectors = TRUE, col ="black")
    }
    ### Decide ###
    print(cbind(neuron[,], meta))
    print(paste0("See neuron in CATMAID at: ", catmaid_url))
    decide = readline(paste0("Are you SURE you want to delete: ", skid," in CATMAID instance ", server, "? y=yes, n=no : "))
    if(decide=="y"){
      path = sprintf("/%d/neuron/%d/delete", pid, nid)
      res = catmaid_fetch(path, body = NULL, include_headers = F,
                          conn = conn, ...)
      if(!is.null(res$success)){
        message(res$success)
      }else{
        message(res$error)
      }
      message("This activity has been logged.")
    }
  }
}

#' Get the CATMAID neuron ID that corresponds to the skeleton ID
#'
#' @description Retrieve the neuron IDs for given skeleton IDs. This is typically the skeleton ID + 1, and is often, but not always accurately, kept by CATMAID tracers in the name of a neuron.
#' @param skids a vector of skeleton IDs or argument applicable to catmaid::catmaid_get_neuronid
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_get_neuronid
catmaid_get_neuronid <- function(skids, pid = 1, conn = NULL, ...){
  skids = catmaid_skids(skids, conn = conn, pid = pid, ...)
  if (any(duplicated(skids))) {
    uskids = unique(skids)
    unids = catmaid_get_neuronid(uskids, pid = pid, conn = conn,...)
    res = unids[match(skids, uskids)]
    return(res)
  }
  skids[is.na(skids)] = -1L
  res = lapply(skids,function(skid)
    catmaid::catmaid_fetch(sprintf("/%d/skeleton/%s/neuronname", pid, skid), body = NULL, include_headers = F,
                           conn = conn, ...)$neuronid)
  res = sapply(res,function(r) ifelse(is.null(r),NA,r))
  names(res) = skids
  res
}

#' Upload neuron(s) to CATMAID
#'
#' @description  Uploads neurons to CATMAID, names them and annotates them.
#' Please use with caution, as you could be heavily adding to a live tracing environment.
#' @param neurons neuronlist object to upload to the CATMAID instance specified when you use catmaid::catmaid_login(), if conn is NULL, else specified by conn$server
#' @param tolerance how many potentially duplictated nodes you will tolerate
#' @param name whatever you want to name your uploaded neurons. If a single character, then it will be added to all uploaded neurons. Else, can be a character vector the same length as swc.
#' @param annotations a character vector of annotations, to be added to all of the uploaded neurons
#' @param include.tags whether of not to transfer each neuron's tags to its newly uploaded cognate
#' @param include.connectors whether of not to transfer each neuron's connectors to its newly uploaded cognate
#' @param join whether or not to attempt to join each uploaded neuron to neurons within the new CATMAID instance, based on the placement of join.tags specified using the next argument
#' @param join.tag a single character specifying a tag that has been used to signify a potential merge point
#' @param min_nodes the minimum number of nodes a potential merger skeleton needs to have
#' @param search.range.nm the maimum distance from which to search from the TODO point tag to find potential mergers
#' @param fafbseg whether or not to use fafbseg::read_brainmaps_meshes on the TODO tag locations and restrict possible merges to neurons within the search range and within the cognate auto-segmented volume
#' @param brain the brain to plot while visualising potential mergers using rgl. Defaults to NULL, no brain plotted.
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param max.upload the maximum number of files that the function will allow you to upload at once
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_controlled_upload
catmaid_controlled_upload <- function(neurons, tolerance = 0.5, name = "v14-seg neuron upload", annotations = "v14-seg upload",
                                      include.tags = TRUE, include.connectors = TRUE,
                                      search.range.nm = 100, join = FALSE, join.tag = "TODO",
                                      fafbseg = FALSE, min_nodes = 2, search.range.nm = 1000,
                                      brain = NULL,
                                      pid = 1, conn = NULL, ...){
  if(!is.neuronlist){
    stop("A neuronlist object must be given to the neurons argument")
  }
  if(length(name)!=1&length(name)!=length(neurons)){
    stop("The name argument must be a chracter vector the same length as neurons,
         or if all uploaded neurons are to be named the same, a character vector of length one")
  }
  nat::nopen3d()
  for(i in 1:length(neurons)){
    neuron = neurons[[i]]
    old.skid = neurons[i,"skid"]
    if(length(name)==1){
      nam = names
    }else{
      nam = names[i]
    }
    message("Assessing whether upload of neuron ",i, " may cause a duplication:" )
    dupe = catmaid_duplicated(neuron, tolerance = tolerance, pid = pid, conn = conn, ...)
    if(dupe){
      warning("Neuron ", i, " with skid ", old.skid, " appears to already exist in the CATMAID instance to which you are seeking to upload.
              Upload for neuron ", i, " aborted.")
      next
    }
    message("Uploading neuron ", i)
    new.skid = catmaid_upload_neurons(swc=neuron,name=nam,annotations=annotations,
                                    include.tags=include.tags,include.connectors=include.connectors,
                                    search.range.nm=10, return.new.skids = TRUE,
                                    conn = conn, pid = pid, max.upload = 1, ...)
    new.neuron = read.neurons.catmaid(new.skid, conn=conn, pid=pid, ...)
    message("Upload successful, neuron ", new.skid, " created, named: ", new.neuron[,"name"])
    if(join){
      message("Finding tagged potential join points ...")
      TODO = catmaid_get_tag(x = new.neuron, tag = join.tag, url = FASE)
      if(nrow(TODO)>1){
        possible.merges = catmaid_find_likely_merge <- function(TODO, skid, pid=1, fafbseg = fafbseg,
                                              min_nodes = min_nodes, search.range.nm = search.range.nm)
        if(nrow(possible.merges)>1){
          catmaid_interactive_join <- function(possible.merges=possible.merges, brain = brain, pid = pid, conn = conn, ...)
        }else{
          message("No potential join sites found for new neuron ", new.skid)
          next
        }
      }else{
        message("No potential join sites found for new neuron ", new.skid)
        next
      }
    }
  }
}


#' Lock or unlock a CATMAID neuron reconstruction
#'
#' @description  Lock or unlock a CATMAID neuron reconstruction by adding or removing a 'locked' annotation to a set of skeleton IDs (skids). A locked neuron cannot be edited until it is unlocked.
#' @param skids the skeleton IDs neurons you wish to lock / unlock
#' @param server the CATMAID server from which to delete the neuron identified by skid
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_lock_neurons
catmaid_lock_neurons <- function(skids, pid = 1, conn = NULL, ...){
  catmaid::catmaid_set_annotations_for_skeletons(skids, annotations = "locked", pid = 1,
                                        conn = NULL, ...)
}
#' @export
#' @rdname catmaid_lock_neurons
catmaid_unlock_neurons <- function(skids, pid = 1, conn = NULL, ...){
  catmaid::catmaid_remove_annotations_for_skeletons(skids, annotations = "locked", pid = 1,
                                                 conn = NULL, ...)
}

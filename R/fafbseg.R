# Function to use the FAFB segmentations in R

#' Calculate a connectivity similarity score between two connectivity profiles
#'
#' @description  Set the local path to the location where you have .zip files of FAFB segmented skeletons
#' @param path path to FAFB segmentation skeleton .zip files
#' @export
#' @rdname set_segmentation_location
set_segmentation_location <- function(path){
  options(fafbseg.skelziproot = path)
}

#' Read neurons from the FAFB segmentation instance (temporary)
#'
#' @description  Reads neurons (or node count data) from the temporary FAFB segmentation instance.
#' containing auto-segmented fragments from Peter Li. Neurons can also be read using a neuroglancer server fafbseg::read.neurons.brainmaps
#' or fafbseg::read_segments2 if from local .zip files
#' @param x typically a skeleton ID in from the CATMAID FAFB v14-seg instance neuroglancer ID.
#' @param google if TRUE, we treat x as the name of a Google FAFB segment. If name = FALSE, then x can be skeleton ids (numeric or characters) or the name / annotation for a neuron,
#' if given as a character starting with 'name:' or 'annotation:' respectively. When Google FAFB fragments are merged, one name is maintained and the other may be found as an annotation. In order to check for this, annotations can also be read.
#' @param read.from from where to read FAFB segmented skeletons to generate node count data
#' @param conn a CATMAID connection object. If NULL catmaid::catmaid::catmaid_login() is called.
#' @param ... methods passed to catmaid::read.neurons.catmaid
#' @export
#' @rdname read.neurons.fafbseg
read.neurons.fafbseg <- function(x, google = FALSE, conn = NULL, ...){
  if(google){
    x = paste0("name:",x)
    y = paste0("annotation:",x)
    skids = unique(unlist(sapply(c(x,y),catmaid::catmaid_skids,several.ok=FALSE,conn=conn, ...)))
  }else{
    skids = x
  }
  if(is.null(conn)){
    conn = catmaid::catmaid_login()
  }
  if(conn$server != "https://neuropil.janelia.org/tracing/fafb/v14/"){
    message("You need to log into CATMAID: https://neuropil.janelia.org/tracing/fafb/v14/")
    message("See ?catmaid_login")
  }
  conn$server = "https://neuropil.janelia.org/tracing/fafb/v14-seg/"
  n = catmaid::read.neurons.catmaid(skids, conn=conn, OmitFailures = TRUE,...)
  n[,"nodes"] = nat:::summary.neuronlist(n)$nodes
  n[,"skeleton.type"] = "FAFB-seg"
  n
}

#' Use any CATMAID function in the v14-seg environment, without logging into separate CATMAID instance
#'
#' @description Use any CATMAID function in the v14-seg environment, without logging into separate CATMAID instance
#' @param FUN catmaid function from rcatmaid or catnat
#' @param ... methods passed to FUN
#' @export
#' @rdname fafb_seg
fafb_seg <- function(FUN, ...){
  conn = catmaid::catmaid_login()
  if(conn$server != "https://neuropil.janelia.org/tracing/fafb/v14/"){
    message("You need to log into CATMAID: https://neuropil.janelia.org/tracing/fafb/v14/")
    message("See ?catmaid_login")
  }
  conn$server = "https://neuropil.janelia.org/tracing/fafb/v14-seg/"
  FUN(conn=conn, ...)
}

#' Log into the v14-seg CATMAID instance for FAFB v14 Adult flybrain segmented skeletonisations
#'
#' @description Log into the v14-seg CATMAID instance for FAFB v14 Adult flybrain segmented skeletonisations from Peter Li at Google
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname fafb_seg_conn
fafb_seg_conn <- function(pid = 1, conn = NULL, ...){
  if(is.null(conn)){
    conn = catmaid::catmaid_login()

  }
  if(conn$server != "https://neuropil.janelia.org/tracing/fafb/v14/"){
    warning("You need to log into CATMAID: https://neuropil.janelia.org/tracing/fafb/v14/")
    warning("See ?catmaid_login")
  }
  conn$server = "https://neuropil.janelia.org/tracing/fafb/v14-seg/"
  conn
}

#' @export
#' @rdname read.neurons.fafbseg
fafbseg_get_node_count <-function(x, read.from = c("CATMAID","Neuroglancer","local"), ...){
  read.from=match.arg(read.from)
  if(read.from=="CATMAID"){
    y = paste0("name:",x)
    conn = catmaid::catmaid_login()
    if(conn$server != "https://neuropil.janelia.org/tracing/fafb/v14/"){
      message("You need to log into CATMAID: https://neuropil.janelia.org/tracing/fafb/v14/")
      message("See ?catmaid_login")

    }
    conn$server = "https://neuropil.janelia.org/tracing/fafb/v14-seg/"
    skids = sapply(y,catmaid_skids,several.ok=FALSE,conn=conn, ...)
    skids[sapply(skids,length)==0] = 0
    skids = unlist(skids)
    n = catmaid::catmaid_get_node_count(skids, conn=conn,OmitFailures = TRUE,...)
    nc = data.frame(ngl_id = x, skid = skids, nodes = n)
  }else if (read.from=="Neuroglancer"){
    rs = fafbseg::read.neurons.brainmaps(x, OmitFailures = TRUE, ...)
    nc = nat:::summary.neuronlist(rs)
  }else{
    rs= fafbseg::read_segments2(x,  OmitFailures = TRUE, ...)
    nc = nat:::summary.neuronlist(rs)
    rownames(nc) = rs[,"segment"]
  }
  nc
}

#' Find out what is already known about a neuron's connectivity profile and connected FAFB segments
#'
#' @description  Set the local path to the location where you have .zip files of FAFB segmented skeletons
#' @param skids neuron skeleton ids
#' @param direction whether to fetch putative incoming or outgoing partners
#' @param volume volume to which to restrict connector locations
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid_get_connector_table
#' @export
#' @rdname fafb_neuron_details
fafb_neuron_details <- function(skids, direction = c("incoming","outgoing"), connector_ids = connector_ids, volume = NULL, pid = 1, conn = NULL, ...) {
  if(!requireNamespace('fafbseg', quietly = TRUE))
    stop("Please install suggested fafbseg package")
  direction=match.arg(direction)
  if(direction=="incoming"){
    connected=catmaid::catmaid_get_connectors_between(post_skids = skids, pid = 1, conn = NULL,...)
    if(length(skids)>15){
      connected=do.call(rbind,lapply(skids,function(skid) catmaid::catmaid_get_connectors_between(post_skids = skid, pid = 1, conn = NULL,...)))
    }
    connected = connected[, c("connector_id", "pre_node_x",
                              "pre_node_y", "pre_node_z")]
  }else{
    connected=catmaid::catmaid_get_connectors_between(pre_skids = skids, pid = 1, conn = NULL,...)
    if(length(skids)>15){
      connected=do.call(rbind,lapply(skids,function(skid) catmaid::catmaid_get_connectors_between(pre_skids = skid, pid = 1, conn = NULL,...)))
    }
    connected = connected[, c("connector_id", "post_node_x",
                              "post_node_y", "post_node_z")]
  }
  colnames(connected) = c("connector_id","X","Y","Z")
  if(!is.null(volume)){
    i = nat::pointsinside(nat::xyzmatrix(connected),volume,rval = "logical")
    connected = connected[i,]
    if(is.null(connector_ids)){
      connector_ids = connected[,"connector_id"]
    }
  }
  known = catmaid::catmaid_get_connector_table(skids, direction = direction,get_partner_names = TRUE, get_partner_nodes = TRUE, pid = pid, conn = conn, ...)
  known = known[match(connected$connector_id,known$connector_id),]
  df = cbind(known,connected[,-1])
  if(!is.null(connector_ids)){
    df = subset(df,connector_id%in%connector_ids)
  }
  df[is.na(df$X),c("X","Y","Z")] = df[is.na(df$X),c("x","y","z")] # In the cases where there are no pre nodes
  df$ngl_id = fafbseg::brainmaps_xyz2id(df[,c('X','Y', 'Z')])
  df
}

#' Fetch the up or outgoing auto-traced FAFB fragments from a CATMAID FAFB neuron
#'
#' @description  Fetch skeleton ids or read skeletons from the FAFb Google Segmentation by Peter Li.
#' This is done by mapping the location of connector nodes in FAFB to the volumetric Google segmentations,
#' and then their cognate skeletons. Relies on package fafbseg and brainmaps authentication, user list is curated.
#' @param skids neuron skeleton ids
#' @param ids ids for FAFB segmentations to be read. If skids are given, fafb_frags_ids is called and ids is overlooked.
#' @param direction whether to fetch putative incoming or outgoing partners
#' @param volume volume to which to restrict connector locations
#' @param connector_ids restrict your search to only certain connectors. Use if, for example, you want to spatially restrict your search.
#' Default set to NULL, searches all incoming or outgoing connectors, as specified by direction.
#' @param max.nodes the maximum number of nodes that an extant FAFB partner can have, before we consider using the segmentation for our tracing list.
#' @param add.links add links to CATMAID in our tracing list
#' @param treat.skids.separately create hitlist, with hits pooled for all skids given (FALSE, default) or per skid given (TRUE)
#' @param read.from from where to read FAFB segmented skeletons
#' @param unique if TRUE, fafb_seg_tracing_list gives each segment only once in the tracing list, the rest of the information is just one example of a putative connection out of the total number, given in the 'hits' column
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch, catmaid::catmaid_get_connector_table for fafb_frags_ids, and read read methods for fafb_frags_skeletons
#' @details fafb_frags_ids returns Neuroglancer IDs for FAFB segments up or downstream of the specified FAFB CATMAID skeleton IDs.
#' fafb_frags_skeletons reads neurons from Neuroglancer IDs or calls fafb_frags_ids, using either saved skeletons, Neuroglancer brainmaps access or CATMAID access.
#' fafb_seg_hitlist generates a ranked hitlist of fragments from the given skids, either up or downstream.
#' fafb_seg_tracing_list goes a bit further and supplies emt information and links to \code{FAFBv14} and the \code{FAFBv14} segmentation instance.
#' @export
#' @rdname fafb_frags
fafb_frags_ids <- function(skids, direction = c("incoming","outgoing"), connector_ids = NULL, volume = NULL, pid = 1, conn = NULL, ...) {
  if(!requireNamespace('fafbseg', quietly = TRUE))
    stop("Please install suggested fafbseg package")
  direction=match.arg(direction)
  if(direction=="incoming"){
    connected=catmaid::catmaid_get_connectors_between(post_skids = skids, pid=pid, conn = conn, ...)
    if(length(skids)>15){
      connected=do.call(rbind,lapply(skids,function(skid) catmaid::catmaid_get_connectors_between(post_skids = skid, pid=pid, conn = conn, ...)))
    }
    connected = connected[,c("connector_id", "pre_node_x", "pre_node_y", "pre_node_z")]
  }else{
    connected=catmaid::catmaid_get_connectors_between(pre_skids = skids, ...)
    if(length(skids)>15){
      connected=do.call(rbind,lapply(skids,function(skid) catmaid::catmaid_get_connectors_between(pre_skids = skid, pid=pid, conn = conn, ...)))
    }
    connected = connected[,c("connector_id", "post_node_x", "post_node_y", "post_node_z")]
  }
  if(!is.null(connector_ids)){
    connected = subset(connected,connector_id%in%connector_ids)
  }
  colnames(connected) = c("connector_id","X","Y","Z")
  if(!is.null(volume)){
    i = nat::pointsinside(nat::xyzmatrix(connected),volume,rval = "logical")
    connected = connected[i,]
  }
  if(!is.null(connector_ids)){
    connected = subset(connected,connector_id%in%connector_ids)
  }
  connected_ids= fafbseg::brainmaps_xyz2id(connected[,c('X','Y', 'Z')])
  connected_ids
}
#' @export
#' @rdname fafb_frags
fafb_frags_skeletons <- function(ids, skids = NULL, direction = c("incoming","outgoing"),
                                 connector_ids = NULL, read.from = c("CATMAID","Neuroglancer","local"),
                                 pid = 1, conn = NULL, ...) {
  direction=match.arg(direction)
  read.from=match.arg(read.from)
  if(skids){
    ids = fafb_frags_ids(skids = skids, direction = direction, connector_ids = connector_ids, pid=pid, conn = conn, ...)
  }
  uids=setdiff(ids,0)
  if(read.from=="CATMAID"){
    rs = read.neurons.fafbseg(uids, google = TRUE, OmitFailures = TRUE, ...)
  }else if (read.from=="Neuroglancer"){
    rs = fafbseg::read.neurons.brainmaps(x, OmitFailures = TRUE, ...)
  }else{
    rs= fafbseg::read_segments2(uids,  OmitFailures = TRUE, ...)
  }
  tt=table(ids)
  conndf=data.frame(ngl_id=as.numeric(names(tt)), hits=unname(unclass(tt)))
  m=merge(rs[,], conndf, by='ngl_id', sort = FALSE)
  rs[,]=m
  rs
}
#' @export
#' @rdname fafb_frags
fafb_seg_hitlist <- function(skids, direction = c("incoming","outgoing"),
                             connector_ids = NULL, treat.skids.separately = FALSE,
                             pid=pid, conn = conn, ...){
  direction=match.arg(direction)
  if(!requireNamespace('reshape2', quietly = TRUE))
    stop("Please install suggested reshape2 package")
  func <- function(skids = skids, direction = direction, connector_ids = connector_ids, pid=pid, conn = conn, ...){
    ids = fafb_frags_ids(skids = skids, direction = direction, connector_ids = connector_ids, pid=pid, conn = conn, ...)
    df = reshape2::melt(table(ids))
    df = df[order(df$value, decreasing = TRUE),]
    colnames(df) = c("ngl_id","hits")
    df
  }
  if(treat.skids.separately){
    hitlist = lapply(skids, function(skid) func(skids = skid, direction = direction, connector_ids = connector_ids, pid=pid, conn = conn, ...))
    nams = rep(skids,sapply(hitlist,nrow))
    df = do.call(rbind, hitlist)
    df$skid = nams
  }else{
    df = func(skids = skids, direction = direction, connector_ids = connector_ids, pid=pid, conn = conn, ...)
  }
  df
}

#' @export
#' @rdname fafb_frags
fafb_seg_tracing_list <- function(skids, direction = c("incoming","outgoing"),
                               connector_ids = NULL, max.nodes = 10,
                               add.links = TRUE, unique = TRUE, ...){
  direction=match.arg(direction)
  df = fafb_neuron_details(skids = skids, direction = direction, connector_ids = connector_ids)
  if(!is.null(max.nodes)){
   df = subset(df,as.numeric(partner_nodes)<=max.nodes)
  }
  if (nrow(df)>0){
    if(add.links){
      df$FAFB.link = connector_URL(df, server = "https://neuropil.janelia.org/tracing/fafb/v14/")
      df$FAFBseg.link = connector_URL(df, server = "https://neuropil.janelia.org/tracing/fafb/v14-seg/")
    }
    df$ngl_id[df$ngl_id==0] = paste0(df$ngl_id[df$ngl_id==0],"_",1:sum(df$ngl_id==0))
    hits = table(df$ngl_id)
    hits["0"] = NA
    df$hits = hits[as.character(df$ngl_id)]
    df = df[order(df$hits, decreasing = TRUE),]
    if(unique){
      df = df[!duplicated(df$ngl_id),]
      df$entry = "example_connection_to_ngl_id"
    }
  }
  df
}

#' Find out what is already known about a neuron's connectivity profile and connected FAFB segments
#'
#' @description  Set the local path to the location where you have .zip files of FAFB segmented skeletons. This function also very roughly estimates
#' whether fragment is microtubule containing,  what Strahler order it is at, or if it is axonic or dendritic, if the neuron has this marked
#' @param someneuronlist a neuronlist or neuron object
#' @param node.match how many nodes of each neuron in someneuronlist, need to be within a auto segmented volume, for it to be said to match.
#' These nodes all need to be consecutive, in the sense that they must be in the same segment or a branch from that segment. I.e. If a neuron matches with a volume
#' 5 times at diverse points across it arbour, this is thought to be a non-match with a large, proximal auto-traced segment.
#' need be in the volumetric Google FAFB segmentation for a Neuroglancer fragment, for that fragment to be returned.
#' @param return.unmatched defaults to FALSE. If TRUE, then a data frame of unmatched Point Numbers for the neuron in question are returned, and their Strahler order.
#' @param ... methods passed to fafbseg::brainmaps_xyz2id
#' @export
#' @rdname map_fafbsegs_to_neuron
map_fafbsegs_to_neuron <- function(someneuronlist, node.match = 5, return.unmatched = FALSE, ...){
  if(is.neuron(someneuronlist)){
    someneuronlist = as.neuronlist(someneuronlist)
  }
  t = data.frame()
  if(nat::is.neuronlist(someneuronlist)){
    pb <- utils::txtProgressBar(min = 0, max = length(someneuronlist), style = 3)
    for(n in 1:length(someneuronlist)){
      message(names(someneuronlist)[n])
      neuron = someneuronlist[[n]]
      segs = fafbseg::brainmaps_xyz2id(nat::xyzmatrix(neuron), ...)
      m = reshape2::melt(table(segs))
      colnames(m) = c("ngl_id","node_hits")
      m$skid = names(someneuronlist)[n]
      # Get a rough idea if a fragment if microtubule containing, or what Strahler order it is at, if the neuron has this marked
      if(!is.null(neuron$d$microtubules)){
        mt = aggregate(neuron$d$microtubules,list(ngl_id = segs),function(x) (sum(x)/length(x))>=0.5)
        colnames(mt) = c("ngl_id","microtubules")
        m = merge(m,mt,all = TRUE)
      }else{
        m$microtubules = NA
      }
      if(!is.null(neuron$d$strahler_order)){
        so = aggregate(neuron$d$strahler_order,list(ngl_id = segs),function(x) names(sort(-table(x)))[1])
        colnames(so) = c("ngl_id","strahler_order")
        m = merge(m,so,all = TRUE)
      }else{
        m$strahler_order = NA
      }
      if(!is.null(neuron$d$Label)){
        lab = aggregate(neuron$d$Label,list(ngl_id = segs),function(x) names(sort(-table(x)))[1])
        colnames(lab) = c("ngl_id","Label")
        m = merge(m,lab,all = TRUE)
      }else{
        m$Label = NA
      }
      mm = subset(m,node_hits>=node.match)
      s = neuron$SegList
      keep = c()
      for(nid in m$ngl_id){
        if(nid!=0){
          pnos = which(segs==nid)
          in.segs = lapply(s,function(y) pnos%in%y)
          in.segs.sum = sapply(in.segs,sum)
          for(pos in which(in.segs.sum>0)){
            branches = sapply(s,function(x) sum(s[[pos]]%in%x))
            score = sum(in.segs.sum[which(branches>0)])
            if(score>=node.match){
              keep = c(keep,nid)
            }
          }
        }
      }
      mm = subset(mm,ngl_id%in%keep)
      if(return.unmatched){
        as = assign_strahler(neuron)
        unmatched.pnos = neuron$d$PointNo[which(!segs%in%mm$ngl_id)]
        mm = subset(as$d,PointNo%in%unmatched.pnos)
        mm$index = rownames(mm)
        mm$skid = neuron$skid
      }
      t = rbind(t,mm)
      utils::setTxtProgressBar(pb, n)
    }
    close(pb)
  }else {
    t = NULL
  }
  t
}


# write_marked_swc <- function(nl, files,
#                              dir = '/GD/LMBD/Papers/2017pns/fig/Alex/images/fafbseg/data/whole/',
#                              format="swc"){
#   for(n in 1:length(nl)){
#     file = paste0(dir,files[n])
#
#   }
#   write(nl[[n]], file = file,
#         ncolumns = if(is.character(x)) 1 else 5,
#         append = FALSE, sep = " ")
#
# }


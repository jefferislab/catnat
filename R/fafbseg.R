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

#' Read neurons from the FAFB segmentation instance (temportary)
#'
#' @description  Reads neurons (or node count data) from the temporary FAFB segmentation instance.
#' containing auto-segmented fragments from Peter Li. Neurons can also be read using a neuroglancer server fafbseg::read.neurons.brainmaps
#' or fafbseg::read_segments2 if from local .zip files
#' @param x typically a neuroglancer ID. By default, we treat it as the name of a Google FAFB segment. If name = FALSE, then x can be skeleton ids (numeric or characters) or the name / annotation for a neuron,
#' if given as a character starting with 'name:' or 'anotation:' respectively.
#' @param name when TRUE, "name:" is added in front of ever element in x
#' @param and.annotation when Google FAFB fragments are merged, one name is maintained and the other may be found as an annotation. In order to check for this, annotations can also be read
#' @param read.from from where to read FAFB segmented skeletons to generate node count data
#' @param ... methods passed to catmaid::read.neurons.catmaid
#' @export
#' @rdname read.neurons.fafbseg
read.neurons.fafbseg <- function(x, name = TRUE, and.annotation = TRUE, ...){
  require(catmaid)
  if(name){
    x = paste0("name:",x)
    y = paste0("annotation:",x)
    skids = unique(unlist(sapply(c(x,y),catmaid_skids,several.ok=FALSE,conn=conn, ...)))
  }else{
    skids = x
  }
  conn = catmaid::catmaid_login()
  if(conn$server != "https://neuropil.janelia.org/tracing/fafb/v14/"){
    message("You need to log into CATMAID: https://neuropil.janelia.org/tracing/fafb/v14/")
    message("See ?catmaid_login")

  }
  conn$server = "https://neuropil.janelia.org/tracing/fafb/v14-seg/"
  n = catmaid::read.neurons.catmaid(skids, conn=conn,OmitFailures = TRUE,...)
  n[,"nodes"] = nat:::summary.neuronlist(n)$nodes
  n[,"skeleton.type"] = "FAFB-seg"
  n
}

#' @export
#' @rdname read.neurons.fafbseg
fafbseg_get_node_count <-function(x, read.from = c("CATMAID","Neuroglancer","local"), ...){
  read.from=match.arg(read.from)
  if(read.from=="CATMAID"){
    require(catmaid)
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
#' @param ... methods passed to catmaid_get_connector_table
#' @export
#' @rdname set_segmentation_location
fafb_neuron_details <- function(skids, direction = c("incoming","outgoing"), connector_ids = connector_ids, ...) {
  require(fafbseg)
  direction=match.arg(direction)
  if(direction=="incoming"){
    connected=catmaid::catmaid_get_connectors_between(post_skids = skids, ...)
    connected = connected[,c("connector_id", "pre_node_x", "pre_node_y", "pre_node_z")]
  }else{
    connected=catmaid::catmaid_get_connectors_between(pre_skids = skids, ...)
    connected = connected[,c("connector_id", "post_node_x", "post_node_y", "post_node_z")]
  }
  colnames(connected) = c("connector_id","X","Y","Z")
  known = catmaid::catmaid_get_connector_table(skids, direction = direction,get_partner_names = TRUE, get_partner_nodes = TRUE)
  known = known[match(connected$connector_id,known$connector_id),]
  df = cbind(known,connected[,-1])
  if(!is.null(connector_ids)){
    df = subset(df,connector_id%in%connector_ids)
  }
  df[is.na(df$X),c("X","Y","Z")] = df[is.na(df$X),c("x","y","z")] # In the cases where there are no pre nodes
  df$segment = fafbseg::brainmaps_xyz2id(df[,c('X','Y', 'Z')])
  df
}

#' Fetch the up or outgoing auto-traced FAFB fragments from a CATMAID FAFB neuron
#'
#' @description  Fetch skeleton ids or read skeletons from the FAFV Google Segmentation by Peter Li.
#' This is done by mapping the location of connector nodes in FAFB to the volumetric Google segmentations,
#' and then their cognate skeletons. Relies on package fafbseg and brainmaps authentication, user list is curated.
#' @param skids neuron skeleton ids
#' @param ids ids for FAFB segmentations to be read. If skids are given, fafb_frags_ids is called and ids is overlooked.
#' @param direction whether to fetch putative incoming or outgoing partners
#' @param connector_ids restrict your search to only certain connectors. Use if, for example, you want to spatially restrict your search.
#' Default set to NULL, searches all incoming or outgoing connectors, as specified by direction.
#' @param max.nodes the maximum number of nodes that an extant FAFB partner can have, before we consider using the segmentation for our tracing list.
#' @param add.links add links to CATMAID in our tracing list
#' @param unique if TRUE, fafb_seg_tracing_list gives each segment only once in the tracing list, the rest of the information is just one example of a putative connection out of the total number, given in the 'hits' column
#' @param ... methods passed to catmaid::catmaid_get_connector_table for fafb_frags_ids, and read read methods for fafb_frags_skeletons
#' @details fafb_frags_ids returns Neuroglancer IDs for FAFB segments up or downstream of the specified FAFB CATMAID skeleton IDs.
#' fafb_frags_skeletons reads neurons from Neuroglancer IDs or calls fafb_frags_ids, using either saved skeletons, Neuroglancer brainmaps access or CATMAID access.
#' fafb_seg_hitlist generates a ranked hitlist of fragments from the given skids, either up or downstream.
#' fafb_seg_tracing_list goes a bit further and supplies emt information and links to FAFBv14 and the FAFBv14-segmentation instance.
#' @export
#' @rdname fafb_frags
fafb_frags_ids <- function(skids, direction = c("incoming","outgoing"), connector_ids = NULL, ...) {
  require(fafbseg)
  direction=match.arg(direction)
  if(direction=="incoming"){
    connected=catmaid::catmaid_get_connectors_between(post_skids = skids, ...)
    connected = connected[,c("connector_id", "pre_node_x", "pre_node_y", "pre_node_z")]
  }else{
    connected=catmaid::catmaid_get_connectors_between(pre_skids = skids, ...)
    connected = connected[,c("connector_id", "post_node_x", "post_node_y", "post_node_z")]
  }
  if(!is.null(connector_ids)){
    connected = subset(connected,connector_id%in%connector_ids)
  }
  colnames(connected) = c("connector_id","X","Y","Z")
  connected_ids= fafbseg::brainmaps_xyz2id(connected[,c('X','Y', 'Z')])
  connected_ids
}
#' @rdname fafb_frags
fafb_frags_skeletons <- function(ids, skids = NULL, direction = c("incoming","outgoing"),
                                 connector_ids = NULL, read.from = c("CATMAID","Neuroglancer","local"), ...) {
  direction=match.arg(direction)
  read.from=match.arg(read.from)
  if(skids){
    ids = fafb_frags_ids(skids = skids, direction = direction, connector_ids = connector_ids)
  }
  uids=setdiff(ids,0)
  if(read.from=="CATMAID"){
    rs = read.neurons.fafbseg(uids,name = TRUE, OmitFailures = TRUE, ...)
  }else if (read.from=="Neuroglancer"){
    rs = fafbseg::read.neurons.brainmaps(x, OmitFailures = TRUE, ...)
  }else{
    rs= fafbseg::read_segments2(uids,  OmitFailures = TRUE, ...)
  }
  tt=table(ids)
  conndf=data.frame(segment=as.numeric(names(tt)), hits=unname(unclass(tt)))
  m=merge(rs[,], conndf, by='segment', sort = FALSE)
  rs[,]=m
  rs
}
#' @rdname fafb_frags
fafb_seg_hitlist <- function(skids, direction = c("incoming","outgoing"),
                             connector_ids = NULL, ...){
  ids = fafb_frags_ids(skids = skids, direction = direction, connector_ids = connector_ids, ...)
  df = reshape2::melt(table(ids))
  df = df[order(df$value, decreasing = TRUE),]
  colnames(df) = c("ngl_ids","hits")
  df
}

#' @rdname fafb_frags
fafb_seg_tracing_list <- function(skids, direction = c("incoming","outgoing"),
                               connector_ids = NULL, max.nodes = 10,
                               add.links = TRUE, unique = TRUE, ...){
  df = fafb_neuron_details(skids, direction = direction, connector_ids = connector_ids)
  df = subset(df,as.numeric(partner_nodes)<=max.nodes)
  if (nrow(df)>0){
    if(add.links){
      df$FAFB.link = connector_URL(df, server = "https://neuropil.janelia.org/tracing/fafb/v14/")
      df$FAFBseg.link = connector_URL(df, server = "https://neuropil.janelia.org/tracing/fafb/v14-seg/")
    }
    df$segment[df$segment==0] = paste0(df$segment[df$segment==0],"_",1:sum(df$segment==0))
    hits = table(df$segment)
    hits["0"] = NA
    df$hits = hits[as.character(df$segment)]
    df = df[order(df$hits, decreasing = TRUE),]
    if(unique){
      df = df[!duplicated(df$segment),]
    }
  }
  df
}

upload_swc_to_catmaid <- function (swc, name ="ASB upload", pid = 1, conn = NULL, ...) {
  post_data = list()
  post_data[sprintf("file[%d]", seq_along(swc))] = as.list(swc)
  post_data[sprintf("name[%d]", seq_along(annotations))] = as.list(name)
  path = sprintf("/%d/skeletons/import", pid)
  res = catmaid_fetch(path, body = post_data, include_headers = F,
                      simplifyVector = T, conn = conn, ...)
  invisible(catmaid_error_check(res))
}
# swc = readLines("/GD/LMBD/Papers/2017pns/fig/Alex/Data/tracing/swc/370309397.swc")
# obj_list <- lapply(list(swc),paste,collapse="\r\n")
# obj_vec <- as.vector(obj_list)





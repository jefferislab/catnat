# Functions for sampling up-/downstream of neurons in CATMAID, making use of the Google FAFB segmentation to rank synapses

#' Ranked sampling
#'
#' @description  Generating a ranked sampling spreadsheet based on number of overlapping synapses with segments for CATMAID neurons or neuron lists. The function relies on the brainmaps API from Google. The result depends on marked up pre- or postsynapses.
#' @param x neuron or neuronlist. given by a CATMAID skid, skid list, or annotation. Synaptic partners of these neurons will be returned.
#' @param fileout character. defines the file to write out the generated spreadsheet in csv format
#' @param omitNAs logical. Should NAs be omitted? FALSE by default
#' @param ... catmaid::catmaid_get_connector_table for segments_ranked_pre, catmaid::catmaid_get_connectors_between for segments_ranked_post
#' @return A \code{list}, most importantly the columns: "segment"  - google segment ID, "n" - number of overlaps with query neuron connections, "FAFB.link" and "FAFBseg.link" URLs to position of the connection in the regular and the Google seg CATMAID environments, respectively.
#' @examples
#' \donttest{
#' DA2_up_segs = segments_ranked_pre("glomerulus DA2 right")
#' View(DA2_up_segs)
#' segments_ranked_post(38885, fileout = "single_DA2_downstream_segs.csv", omitNAs = TRUE)
#' }
#' @export
#' @rdname segments_ranked
segments_ranked_pre <- function(x,
                                fileout = NULL,
                                omitNAs = FALSE,
                                volume = "brainmaps://772153499790:fafb_v14:fafb-ffn1-20190805",
                                ...) {
  known_up=catmaid::catmaid_get_connector_table(x,
                                       direction = 'incoming',
                                       get_partner_names = T,
                                       get_partner_nodes = T, ...)
  known_up = known_up[!data.table::duplicated(known_up),]
  known_up$segment = NA
  for (j in 1:ceiling(nrow(known_up)/200)) {
    i = (j*200)-199
    if (i < nrow(known_up)-200) {
      vec = c(i:(i+199))
      known_up$segment[vec]=fafbseg::brainmaps_xyz2id(nat::xyzmatrix(known_up[vec, ]), volume = volume)
    }
    else {
      vec = c(i:nrow(known_up))
      known_up$segment[vec]=fafbseg::brainmaps_xyz2id(nat::xyzmatrix(known_up[vec, ]), volume = volume)
    }
  }

  known_up = dplyr::arrange(known_up, known_up$segment)
  known_up$n = as.vector(rep(table(known_up$segment), table(known_up$segment)))
  known_up[known_up==0] <- NA
  known_up = dplyr::arrange(known_up, dplyr::desc(known_up$n))

  if (isTRUE(omitNAs)){
    na.count = sum(is.na(known_up$segment))
    message(paste(na.count, "NAs omitted"))
    known_up = data.table::na.omit(known_up)
  }

  known_up$FAFB.link = connector_URL(known_up, server = "https://neuropil.janelia.org/tracing/fafb/v14/")
  known_up$FAFBseg.link = connector_URL(known_up, server = "https://neuropil.janelia.org/tracing/fafb/v14-seg-li-190805.0/")

  if(!is.null(fileout)){
    utils::write.csv(known_up, file = fileout)
  }
  else{
    return(known_up)
  }
}


#' @export
#' @rdname segments_ranked
segments_ranked_post <- function(x,
                                 fileout = NULL,
                                 omitNAs = FALSE,
                                 volume = "brainmaps://772153499790:fafb_v14:fafb-ffn1-20190805",
                                 ...) {
  known_down=catmaid::catmaid_get_connectors_between(pre_skids = x, get_names = T, ...)

  known_down$x = as.numeric(unlist(subset(known_down, select = post_node_x)))
  known_down$y = as.numeric(unlist(subset(known_down, select = post_node_y)))
  known_down$z = as.numeric(unlist(subset(known_down, select = post_node_z)))
  known_down <- subset(known_down, select = -c(6:14))

  known_down$segment = NA
  for (j in 1:ceiling(nrow(known_down)/200)) {
    i = (j*200)-199
    if (i < nrow(known_down)-200) {
      vec = c(i:(i+199))
      known_down$segment[vec]=fafbseg::brainmaps_xyz2id(nat::xyzmatrix(known_down[vec, ]), volume = volume)
    }
    else {
      vec = c(i:nrow(known_down))
      known_down$segment[vec]=fafbseg::brainmaps_xyz2id(nat::xyzmatrix(known_down[vec, ]), volume = volume)
    }
  }

  known_down = dplyr::arrange(known_down, known_down$segment)
  known_down$n = as.vector(rep(table(known_down$segment), table(known_down$segment)))
  known_down[known_down==0] <- NA
  known_down = dplyr::arrange(known_down, dplyr::desc(known_down$n))

  if (isTRUE(omitNAs)){
    na.count = sum(is.na(known_down$segment))
    message(paste(na.count, "NAs omitted"))
    known_down = data.table::na.omit(known_down)
  }

  known_down$FAFB.link = connector_URL(known_down, server = "https://neuropil.janelia.org/tracing/fafb/v14/")
  known_down$FAFBseg.link = connector_URL(known_down, server = "https://neuropil.janelia.org/tracing/fafb/v14-seg-li-190805.0/")

  if(!is.null(fileout)){
    utils::write.csv(known_down, file = fileout)
  }
  else{
    return(known_down)
  }
}


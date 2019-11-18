# Functions for sampling up-/downstream of neurons in CATMAID, making use of the Google FAFB segmentation to rank synapses


#' Sampling BLAH
#'
#' @description  BLAH
#' @param x a neuron or neuronlist object, or (a) local file path(s) to your saved .swc file(s).
#' @param fileout whatever you want to name your uploaded neurons. If a single character, then it will be added to all uploaded neurons. Else, can be a character vector the same length as swc.
#' @param omitNAs a character vector of annotations, to be added to all of the uploaded neurons
#' @param ... methods passed to catmaid::catmaid_fetch
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
  known_up = known_up[!duplicated(known_up),]
  known_up$segment = NA
  for (j in 1:ceiling(nrow(known_up)/200)) {
    i = (j*200)-199
    if (i < nrow(known_up)-200) {
      vec = c(i:(i+199))
      known_up$segment[vec]=fafbseg::brainmaps_xyz2id(xyzmatrix(known_up[vec, ]), chunksize=2000, volume = volume)
    }
    else {
      vec = c(i:nrow(known_up))
      known_up$segment[vec]=brainmaps_xyz2id(xyzmatrix(known_up[vec, ]), chunksize=2000, volume = volume)
    }
  }

  known_up = arrange(known_up, known_up$segment)
  known_up$n = as.vector(rep(table(known_up$segment), table(known_up$segment)))
  known_up[known_up==0] <- NA
  known_up = arrange(known_up, desc(known_up$n))

  if (isTRUE(omitNAs)){
    na.count = sum(is.na(known_up$segment))
    message(paste(na.count, "NAs omitted"))
    known_up = na.omit(known_up)
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

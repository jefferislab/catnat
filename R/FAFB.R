#' Read neurons from CATMAID with extra meta-data
#'
#' @description Use \code{read.neurons.catmaid} to acquire neurons from a CATMAID instance.
#' In addition, include as meta-data certain annotations that these neurons have. To select which
#' annotations appear as entries in the neurons' meta-data, select 'meta-annotations' to appear as the
#' data field. By default, meta-annotations for lineage identities are used, to get the lineage / hemilineage
#' to which a FAFB Drosophila neuron belongs.
#'
#' @inheritParams catmaid::read.neurons.catmaid
#' @param meta a vector of meta-annotations to query. The annotations labelled by these meta-annotations will appear as
#' entries in the \code{data.frame} attached to the returned \code{neuronlist}.
#' @param sub what to remove from pulled annotations. Can be set to NULL.
#' @export
#' @rdname read.neurons.catmaid.meta
read.neurons.catmaid.meta <- function(skids,
                                      meta = c("ItoLee_Lineage",
                                               "ItoLee_Hemilineage",
                                               "Hartenstein_Lineage",
                                               "Hartenstein_Hemilineage"),
                                      sub = " |:",
                                      OmitFailures = TRUE,
                                      ...){
  maddf <- data.frame()
  for(ma in meta){
    mad <- catmaid_query_meta_annotations(ma, ...)
    mad$meta <- ma
    mad$field <- gsub(ma, "", mad$name)
    mad$field <- gsub(sub, "", mad$field)
    maddf <- rbind(maddf, mad)
  }
  n <-  catmaid::read.neurons.catmaid(skids, OmitFailures = OmitFailures, ...)
  as <-  catmaid::catmaid_get_annotations_for_skeletons(names(n), ...)
  as <- as[as$id %in% maddf$id, ]
  m <-  merge(maddf, as[,c("skid","id")])
  n[, c(meta, "unique.assignment")] <- NA
  for(skid in n[,"skid"]){
    mm <- m[m$skid==skid,]
    unique.assignment <- TRUE
    mmm <- data.frame()
    for(ma in meta){
      field <- mm[mm$meta==ma,"field"]
      if(length(field)>1){
        field <- paste(field, collapse = "/")
        unique.assignment <-  FALSE
      }else if (length(field)==0){
        field <-  NA
      }
      mmm <- rbind(mmm, data.frame(meta = ma, field = field))
    }
    n[as.character(skid), "unique.assignment"] <- unique.assignment
    n[as.character(skid), meta] <- mmm[match(meta,mmm$meta),"field"]
  }
  n
}

# hidden
catmaid_query_meta_annotations <-function(meta_annotations,
                                          with_annotations = FALSE,
                                          pid=1, conn=NULL,...){
  if(!possibly.numeric(meta_annotations)){
    a <- catmaid::catmaid_get_annotationlist(pid=pid, conn=conn, ...)
    meta_annotations <- a$annotations[a$annotations$name%in%meta_annotations,"id"]
  }
  if(!length(meta_annotations)){
    stop("Please give at least one valid meta annotation or meta annotation ID for your chosen CATMAID instance.")
  }
  post_data <- list()
  post_data[sprintf("annotated_with[%d]", seq_along(meta_annotations))] <- as.list(meta_annotations)
  post_data["with_annotations"] <- with_annotations
  post_data["types"] <- 'annotation'
  path <- sprintf("/%d/annotations/query-targets", pid)
  res <- catmaid::catmaid_fetch(path, body = post_data, include_headers = F,
                       simplifyVector = T, conn = conn, ...)
  invisible(catmaid:::catmaid_error_check(res))
  res$entities
}

# hidden
possibly.numeric <- function(x) {
  stopifnot(is.atomic(x) || is.list(x))
  nNA <- sum(is.na(x))
  nNA.new <- suppressWarnings(sum(is.na(as.numeric(x))))
  nNA.new == nNA
}

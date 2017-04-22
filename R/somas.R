#' Functions for retrieving soma data
#'
#' @description Functions for retrieving data amout somata from skeletons, that has been assigned in CATMAID using the tag system
#' @param ... additional arguments passed to methods.
#'
#' @details CATMAID access required. Data collected and described in cited publication.
#'
#' @references Ohyama T, Schneider-Mizell CM, Fetter RD, Aleman JV, Franconville R, Rivera-Alba M, Mensh BD, Branson KM, Simpson JH, Truman JW, et al. (2015) A multilevel multimodal circuit enhances action selection in Drosophila. Nature.
#' @return Either the soma's 3D coorindates (soma) its id number (somaid) or its posiiton in the neuron's skeleton (somapos)
#' @export
soma<-function(x, ...) UseMethod("soma")

#' @export
#' @rdname soma
soma.neuronlist<-function(x, ...) {
  rdf=plyr::ldply(x, soma, ...)
  rownames(rdf)=rdf[[1]]
  rdf[-1]
}

#' @export
#' @rdname soma
soma.neuron<-function(x, ...) {
  r=if(length(somaid<-x$tags$soma[[1]])){
    x$d[match(somaid, x$d$PointNo),c("X","Y","Z")]
  } else {
    matrix(NA_real_, ncol = 3L, dimnames = list(NULL, c("X","Y","Z") ))
  }
  as.data.frame(r)
}

#' @export
#' @rdname soma
somaid<-function(x, ...) UseMethod("somaid")

#' @export
#' @rdname soma
somaid.neuronlist<-function(x, ...) {
  unlist(nlapply(x, somaid.neuron))
}

#' @export
#' @rdname soma
somaid.neuron <- function(n) unlist(n$tags$soma)

#' @export
#' @rdname soma
somapos<-function(x, ...) UseMethod("somaid")

#' @export
#' @rdname soma
somapos.neuronlist<-function(x, ...) {
  unlist(nlapply(x, somapos))
}

#' @export
#' @rdname soma
somapos.neuron <- function(n) match(somaid(n), n$d$PointNo)



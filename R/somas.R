# Functions concerning soma positions

soma<-function(x, ...) UseMethod("soma")

soma.neuronlist<-function(x, ...) {
  rdf=plyr::ldply(x, soma, ...)
  rownames(rdf)=rdf[[1]]
  rdf[-1]
}

soma.neuron<-function(x, ...) {
  r=if(length(somaid<-x$tags$soma[[1]])){
    x$d[match(somaid, x$d$PointNo),c("X","Y","Z")]
  } else {
    matrix(NA_real_, ncol = 3L, dimnames = list(NULL, c("X","Y","Z") ))
  }
  as.data.frame(r)
}

find.soma <- function (sel3dfun = select3d(), indices = names(db), db = getOption("nat.default.neuronlist"),
                       threshold = 0, invert = FALSE, rval = c("names", "data.frame",
                                                               "neuronlist"))
{
  if (is.null(db))
    stop("Please pass a neuronlist in argument db or set options",
         "(nat.default.neuronlist='myfavneuronlist'). See ?nat for details.")
  if (is.character(db))
    db = get(db)
  selfun = function(x) {
    pointsinside = sel3dfun(soma(x))
    sum(pointsinside, na.rm = T) > threshold
  }
  if (invert)
    selfun = Negate(selfun)
  rval = match.arg(rval)
  subset(db, subset = indices, filterfun = selfun, rval = rval)
}

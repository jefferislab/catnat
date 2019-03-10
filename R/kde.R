#' Convert a 3D Kernel Density Estimate (package ks) into rgl mesh3d object
#'
#' @param x The kde object
#' @param cont The contour level to use
#' @param method Which conversion method to use (see details)
#' @param ashape.alpha A length parameter (alpha) passed to
#'   \code{alphashape3d::\link[alphashape3d]{ashape3d}} (see details).
#' @param ... Additional rgl rendering parameters eventually passed to
#'   \code{rgl::\link[rgl]{tmesh3d}}
#' @param approx.cont Whether to compute an approximate contour level. (Default
#'   \code{FALSE}).
#'
#' @details There are two conversion methods. The default uses the output of the
#'   misc3d::contour3d function. This looks nice but seems to me to have a lot
#'   of redundant triangles. The alternative uses the \code{alphashape3d}
#'   library which produces a simpler output. \bold{But} you must provide a
#'   length parameter \code{ashape.alpha} which is passed to \code{alphashape3d}
#' @return An object of class mesh3d
#' @export
#' @seealso \code{\link[ks]{plot.kde}}, \code{\link[rgl]{as.mesh3d}},
#'   \code{\link[nat]{as.mesh3d.ashape3d}}
#' @examples
#' \dontrun{
#' # Assuming you have lh.voxels object from Alex
#' lh.voxels.m3d=lapply(lh.voxels, function(x) as.mesh3d(attr(x, 'fhat')))
#' shade3d(lh.voxels.m3d[[1]], col='red', alpha=.3)
#' # plot several at once
#' mapply(shade3d, lh.voxels.m3d[c(1,4,10)], col=rainbow(3), alpha=.3)
#'
#' # alternatively making a smoother/simpler alphashape
#' # note that the parameter
#' am1=as.mesh3d(attr(lh.voxels[[1]], 'fhat'), method='alpha', alpha=10)
#' lh.voxels.m3d2=lapply(, function(x) )
#' }
as.mesh3d.kde <- function(x, cont=50, approx.cont=TRUE, method=c("tri","ashape3d"),
                          ashape.alpha=NULL, ...) {
  if(!requireNamespace("ks", quietly = TRUE))
    stop("You must install suggested package ks in order to use as.mesh3d.kde!")
  fhat <- x
  d <- ncol(fhat$x)
  if(!length(d) || d!=3) stop("I need 3d data!")
  if(length(cont)!=1) stop("I need exactly one contour level!")
  if (!is.null(fhat$cont)) {
    cont.ind <- rep(FALSE, length(fhat$cont))
    cont.levels=as.numeric(sub("%", "", names(fhat$cont), fixed = T))
    cont.ind[which(cont == 100 - cont.levels)] <- TRUE
    if (all(!cont.ind))
      hts <- ks::contourLevels(fhat, prob = (100 - cont)/100, approx = approx.cont)
    else hts <- fhat$cont[cont.ind]
  }
  else hts <- ks::contourLevels(fhat, prob = (100 - cont)/100,
      approx = approx.cont)

  tris=misc3d::contour3d(fhat$estimate, level = hts,
                    x = fhat$eval.points[[1]],
                    y = fhat$eval.points[[2]],
                    z = fhat$eval.points[[3]],
                    draw=FALSE)
  method=match.arg(method)
  if(method=="tri") {
    nverts=nrow(tris[[1]])
    verts=do.call(rbind, tris[1:3])
    inds=rbind(1:nverts,(1:nverts)+nverts, (1:nverts)+nverts*2L)
    rgl::tmesh3d(t(verts),indices = inds, homogeneous = F, ...)
  } else {
    as=alphashape3d::ashape3d(xyzmatrix(tris[[1]]), alpha=ashape.alpha)
    rgl::as.mesh3d(as, ...)
  }
}

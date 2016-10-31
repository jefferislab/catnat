#' Jump from R into the L1 larval EM data in CATMAID
#'
#' @description Takes you to a region in the Em data chosen by an rgl selection window by the user
#'
#' @param x a 3D shape, neuronlist or neuron object that has been plotted in 3D whose coorinates can be accessed with nat::xyzmatrix()
#' @param s selection method. Defaults to rgl selection window.
#' @param zoom zoom level to be applied in CATMAID
#' @param open whether to open CATMAID window in browser
#' @param scalefac the scale factor that has been applied to the neurons/shapes plotted in R
#' @param ... additional arguments passed to methods.
#'
#' @details CATMAID access required. Data collected and described in cited publication.
#'
#' @references Ohyama T, Schneider-Mizell CM, Fetter RD, Aleman JV, Franconville R, Rivera-Alba M, Mensh BD, Branson KM, Simpson JH, Truman JW, et al. (2015) A multilevel multimodal circuit enhances action selection in Drosophila. Nature.
#' @return Appropriate L1 CATMAID url
#' @export
#' @rdname open_l1
open_l1=function(x, s = rgl::select3d(), zoom=1, open = interactive(), scalefac=1, ...){
  if (is.vector(x, mode = "numeric") && length(x) == 3) {
    xyz = matrix(x, ncol = 3)
  }
  else {
    xyz = nat::xyzmatrix(x)
    if (nrow(xyz) > 1) {
      xyz = colMeans(xyz[s(xyz), , drop = F])
      xyz = matrix(xyz, ncol = 3)
    }
  }
  xyzi = as.integer(xyz*scalefac)
  url=sprintf("https://neurocean.janelia.org/catmaidL1/?pid=1&zp=%d&yp=%d&xp=%d&tool=tracingtool&sid0=1&s0=%f",
              xyzi[3], xyzi[2], xyzi[1], zoom)
  if (open) {
    browseURL(url)
    invisible(url)
  }
  else {
    url
  }
}

#' Shift neurons in space
#'
#' @description Function for shifting neurons randomly in 3D space. Does not tilt the neuron.
#'
#' @param neurons a neuron/neuronlist object
#' @param shift the size of the shift
#' @param n the number of random shifts to perform
#' @export
#' @rdname shift.neurons
shift.neurons <- function(neurons, shift = 1, n = 100){
  neurons.shifted = neuronlist()
  for(i in 1:n){
    theta = runif(1, min = 0, max = 2*pi)
    z = runif(1, min = -1, max = 1)
    #mag = runif(1, min = 0, max = spatial.scale)
    shift.vector = shift*c(sqrt(1-z^2)*cos(theta),sqrt(1-z^2)*sin(theta),z)
    shift.m = matrix(shift.vector,ncol=3,nrow = nrow(xyzmatrix(neurons)),byrow=TRUE)
    shifted = neurons
    xyzmatrix(shifted) = xyzmatrix(shifted)+shift.m
    names(shifted) = paste0(names(shifted),"_",i)
    neurons.shifted = c(neurons.shifted,shifted)
  }
  neurons.shifted
}

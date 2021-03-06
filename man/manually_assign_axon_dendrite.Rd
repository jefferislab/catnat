% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prune.R
\name{manually_assign_axon_dendrite}
\alias{manually_assign_axon_dendrite}
\alias{manually_assign_axon_dendrite.neuron}
\alias{plot3d.split}
\alias{manually_assign_axon_dendrite.neuronlist}
\title{Manually assign the dendrite and axon to a neuron}
\usage{
manually_assign_axon_dendrite(x, ...)

\method{manually_assign_axon_dendrite}{neuron}(x, ...)

\method{plot3d}{split}(x, soma = TRUE, ...)

\method{manually_assign_axon_dendrite}{neuronlist}(x, ...)
}
\arguments{
\item{x}{a neuron/neuronlist object}

\item{...}{additional arguments passed to methods}

\item{soma}{whether or not to plot a soma, and at what radius}
}
\value{
The neuron/neuronlist object with axon/dendrite info assigned in SWC format to neuron$d
}
\description{
Manually assign the dendrite and axon to neurons / a neuron
}

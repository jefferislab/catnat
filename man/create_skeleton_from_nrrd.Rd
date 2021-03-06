% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_skeleton.R
\name{create_skeleton_from_nrrd}
\alias{create_skeleton_from_nrrd}
\title{Create a skeleton from a skeletonised neuron saved as a .Nrrd file}
\usage{
create_skeleton_from_nrrd(
  files,
  connection.distance = 25,
  k = "all",
  distance.steps = 0.5,
  ...
)
}
\arguments{
\item{files}{paths to saved .Nrrd files that have been skeletonised in Fiji}

\item{connection.distance}{maximal connection distance between points}

\item{k}{value used in nearest neighbours search to identify close points and draw lines between them}

\item{distance.steps}{each round of the algorithm connects the closest points to the start point / leaf node, that are within the search range, starting at distance.steps and increasing by this amount to the maximal connection.distance}

\item{...}{additional arguments passed to methods.}
}
\value{
A neuronlist object
}
\description{
Skeletonise a neuron in Fiji and then use this function to retrieve it as a skeleton
}

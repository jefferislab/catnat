% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alphashape.R
\name{transform3dalphashape}
\alias{transform3dalphashape}
\title{Transform the 3D vertices of an alphashape}
\usage{
transform3dalphashape(ashape, transformations, ...)
}
\arguments{
\item{ashape}{list of alphashape objects}

\item{transformations}{A transformation that can be accepted by nat::xform(). I.e. A registration defined by a matrix, a function, a cmtkreg object, or a character vector specifying a path to one or more registrations on disk}

\item{...}{additional arguments passed to methods}
}
\value{
An alpha shape object
}
\description{
Transform the 3D vertices of an alphashape
}
\references{
Lafarge T, Pateiro-López B, Possolo A, Dunkers J (2014) R Implementation of a Polyhedral Approximation to a 3D Set of Points Using the alpha-Shape. J Stat Softw 56
}
\seealso{
\code{\link{combine.alphashape}} \code{\link{WriteVTKalphashape}}
}

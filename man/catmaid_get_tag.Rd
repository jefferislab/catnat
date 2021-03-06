% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/editcatmaid.R
\name{catmaid_get_tag}
\alias{catmaid_get_tag}
\title{Find the location of specified tags for a CATMAID neuron}
\usage{
catmaid_get_tag(
  x,
  tag = "TODO",
  url = FALSE,
  only.leaves = TRUE,
  conn = NULL,
  pid = 1
)
}
\arguments{
\item{x}{a neuron or neuronlist object}

\item{tag}{a single character specifying which tag to look for. Defaults to TODO}

\item{url}{if TRUE (default) a list of URLs pertaining to specified tag locations are returned. If FALSE, a data.frame subsetted from x$d is returned, reporting treenode ID and X,Y,Z positions for specified tags}

\item{only.leaves}{whether or not to only return leaf nodes with the specified tag}

\item{conn}{CATMAID connection object, see ?catmaid::catmaid_login for details. For making the URL.}

\item{pid}{project id. Defaults to 1. For making the URL.}
}
\description{
Find the location of tags in a CATMAID neuron, either as URLs to the location of a TODO tag in CATMAID or as a data.frame reporting the location and skeleton treenode locations of specified tags.
}

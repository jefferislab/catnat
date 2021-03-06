% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/catmaid.R
\name{catmaid_set_meta_annotations}
\alias{catmaid_set_meta_annotations}
\title{Meta-annotate CATMAID annotations}
\usage{
catmaid_set_meta_annotations(
  meta_annotations,
  annotations,
  pid = 1,
  conn = NULL,
  ...
)
}
\arguments{
\item{meta_annotations}{meta-annotation to add}

\item{annotations}{annotations to meta-annotate}

\item{pid}{project id (default 1)}

\item{conn}{a catmaid_connection objection returned by catmaid_login. I
f NULL (the default) a new connection object will be generated using the values of the catmaid.* package options as described in the help for catmaid_login}

\item{...}{additional arguments passed to methods.}
}
\description{
Meta-annotate a group of CATMAID annotations
}

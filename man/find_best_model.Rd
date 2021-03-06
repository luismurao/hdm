% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_best_model.R
\name{find_best_model}
\alias{find_best_model}
\title{Function to find the best n-dimensional ellipsoid model using Partial Roc as a performance criteria.}
\usage{
find_best_model(this.species, cor_threshold = 0.9, ellipsoid_level = 0.975,
  nvars_to_fit = 3, E = 0.05, RandomPercent = 50, NoOfIteration = 1000)
}
\arguments{
\item{this.species, }{Species Temporal Environment "sp.temporal.env" object see \code{\link[hdm]{extract_by_year}}.}

\item{cor_threshold}{Threshold valuefrom which it is considered that the correlation is high see \code{\link[hdm]{correlation_finder}}.}

\item{ellipsoid_level}{The proportion of points to be included inside the ellipsoid see \code{\link[hdm]{ellipsoidfit}}.}

\item{nvars_to_fit}{Number of variables that will be used to model.}

\item{E}{Amount of error admissible for Partial Roc test (by default =.05). Value should range between 0 - 1. see \code{\link[hdm]{PartialROC}}}

\item{RandomPercent}{Occurrence points to be sampled in randomly for the boostrap of the Partial Roc test \code{\link[hdm]{PartialROC}}.}

\item{NoOfIteration}{Number of iteration for the bootstrapping of the Partial Roc test \code{\link[hdm]{PartialROC}}.}
}
\value{
A "sp.temp.best.model" object with metadata of the best model given the performance of the Partial Roc test.
}
\description{
Function to find the best n-dimensional ellipsoid model using Partial Roc as a performance criteria.
}

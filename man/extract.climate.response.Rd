% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QUALYPSOSS.r
\name{extract.climate.response}
\alias{extract.climate.response}
\title{Extract climate response for one time series z}
\usage{
extract.climate.response(ClimateProjections, predCont, predContUnique,
  nMCMC, lam, uniqueFit, parSmooth = 1, listCR = NULL)
}
\arguments{
\item{ClimateProjections}{matrix of climate projections}

\item{predCont}{matrix of continuous predictor corresponding to the climate projections}

\item{predContUnique}{vector of predictors for which we need fitted climate reponses}

\item{nMCMC}{number of MCMC samples}

\item{lam}{fixed smoothing parameter lambda}

\item{uniqueFit}{logical value indicating if only one fit is applied}

\item{parSmooth}{smoothing parameter \code{spar} in \code{\link[stats]{smooth.spline}}: varies in [0,1]}

\item{listCR}{list of objects for the extraction of the climate response}
}
\value{
list with the following fields:
\itemize{
  \item \strong{phi}: MCMC draws of climate response
  \item \strong{eta}: MCMC draws of deviation from the climate response
  \item \strong{sigma2}: MCMC draws of sigma2
  \item \strong{beta}: MCMC draws of beta
}
}
\description{
Extract climate response for one time series z
}
\author{
Guillaume Evin
}
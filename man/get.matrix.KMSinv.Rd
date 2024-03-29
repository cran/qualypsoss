% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QUALYPSOSS.r
\name{get.matrix.KMSinv}
\alias{get.matrix.KMSinv}
\title{get.matrix.KMSinv 
return the inverse of the square Kac-Murdock-Szego matrix for a rho correlation and n lines/colums}
\usage{
get.matrix.KMSinv(n, rho)
}
\arguments{
\item{n}{nummber of lines/columns of the square matrix}

\item{rho}{correlation parameter in (-1,1)}
}
\value{
n x n Kac-Murdock-Szego matrix
}
\description{
get.matrix.KMSinv 
return the inverse of the square Kac-Murdock-Szego matrix for a rho correlation and n lines/colums
}
\references{
Kac, M., W. L. Murdock, and G. Szego. 1953. 'On the Eigen-Values of Certain Hermitian Forms'
Journal of Rational Mechanics and Analysis 2: 767-800.
}
\author{
Guillaume Evin
}

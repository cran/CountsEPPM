\name{EPPMprob}
\alias{EPPMprob}
\title{
Calculation of vector of probabilities for a extended Poisson process model (EPPM).
}
\description{
Calculates a vector of probabilities given a vector of rates vlambda using the matrix exponential function from Bates and Maechler (2012).
}
\usage{
EPPMprob(vlambda)
}
\arguments{
  \item{vlambda}{
vlambda is a vector of rates of an extended Poisson process.
}
}
\value{
The value returned is a vector of probabilities.
}
\references{
Bates D, Maechler M (2016). Matrix: Sparse and Dense Matrix Classes and Methods. 
R package version 1.2-4, \url{https://CRAN.R-project.org/package=Matrix}. 
}
\author{
David M. Smith <smithdm1@us.ibm.com>
}
\keyword{ distribution }

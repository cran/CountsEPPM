\name{Faddyprob.limiting}
\alias{Faddyprob.limiting}
\title{
Calculation of vector of probabilities for the limiting form of the Faddy distribution.
}
\description{
Given a vector of parameters and a scalar of the maximum count the function calculates the vector of lambdas for the limiting form of a  Faddy distribution applicable to under-dispersed data and returns a vector of probabilities.
}
\usage{
Faddyprob.limiting(parameter, nmax)
}
\arguments{
  \item{parameter}{
A vector of the parameters of the limiting form of a Faddy distribution.
}
  \item{nmax}{
The value of the maximum count.
}
}
\value{
Vector of probabilities
}
\references{
Faddy M, Smith D (2011). Analysis of count data with covariate dependence in both mean and variance. Journal of Applied Statistics, 38, 2683-2694.
}
\author{
David M. Smith <david.m.smith@truvenhealth.com>
}
\examples{
all.counts=c(rep(0,5),352,479,530,291,101,17)

nmax1  <- length(all.counts)
nmax   <- nmax1 - 1

parameter <- c(1.8388023,0.6009881)
names(parameter) <- c('beta0 log(mean)','beta0 log(variance)') 
probability <- Faddyprob.limiting(parameter,nmax) 
print(probability)
}
\keyword{ distribution }

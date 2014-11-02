\name{Faddyprob.limiting}
\alias{Faddyprob.limiting}
\title{
Calculation of vector of probabilities for the limiting form of the Faddy distribution.
}
\description{
Given a vector of parameters and a scalar of the maximum count the function calculates the vector of lambdas for the limiting form of a  Faddy distribution applicable to under-dispersed data and returns a vector of probabilities. This limiting form is described in Faddy and Smith (2011) and it is appropriate for use on count data displaying under dispersion with respect to the Poisson. If the general model of Faddyprob.general is fitted to such under-dispersed data and a large value of b results, possibly with the hessian at the apparent maximum being poorly conditioned, it is possible that the limiting model having one less parameter than the general model will fit better.
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
names(parameter) <- c("beta0 log(mean)","beta0 log(variance)") 
probability <- Faddyprob.limiting(parameter,nmax) 
print(probability)
}
\keyword{ distribution }
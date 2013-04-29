\name{Faddyprob.general}
\alias{Faddyprob.general}
\title{
Calculation of vector of probabilities for a Faddy distribution.
}
\description{
Given a vector of parameters and a scalar of the maximum count the function calculates the vector of lambdas for a Faddy distribution and returns a vector of probabilities.
}
\usage{
Faddyprob.general(parameter, nmax)
}
\arguments{
  \item{parameter}{
A vector of the parameters of the Faddy distribution.
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

parameter <- c(1.8386079,0.6021198,6.0714071)
names(parameter) <- c('beta0 log(mean)','beta0 log(variance)','log(b)') 
probability <- Faddyprob.general(parameter,nmax) 
print(probability)
}
\keyword{ distribution }

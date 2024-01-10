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
Faddy M, Smith D. (2011). Analysis of count data with covariate dependence in both mean and variance. 
\emph{Journal of Applied Statistics}, \bold{38}, 2683-2694. \doi{10.1002/bimj.201100214}.
}
\author{
David M. Smith <dmccsmith@verizon.net>
}
\examples{
all.counts=c(rep(0,5),352,479,530,291,101,17)
nmax   <- length(all.counts) - 1
parameter <- c(exp(53.047752),exp(3.801599),-13.205655)
names(parameter) <- c("a","b","c") 
probability <- Faddyprob.general(parameter,nmax) 
print(probability)
}
\keyword{ distribution }

\name{Model.Faddy}
\alias{Model.Faddy}
\title{
Function for Faddy distribution with log link.
}
\description{
Returns probabilities for a Faddy distribution given inputs of model and parameters.
}
\usage{
Model.Faddy(parameter, model.name, link, covariates.matrix.mean,
 offset.mean, fixed.b, vnmax)
}
\arguments{
  \item{parameter}{ 
A vector of the parameters of the model which is set to initial estimates on function call.
}
  \item{model.name}{
The model being fitted is one of the five 'Poisson', 'negative binomial', 'negative binomial fixed b', 'Faddy distribution', 'Faddy distribution fixed b'. 
}
  \item{link}{ 
Takes one of one values i.e., 'log'. The default is 'log'. 
}
  \item{covariates.matrix.mean}{
A matrix of covariates for the mean where rows are the number of values in listcounts and columns the covariates. This matrix is extracted from the formulae in function CountsEPPM. However, in the accompanying example it is shown how it can be constructed 
independently of function CountsEPPM. 
}
  \item{offset.mean}{
An offset vector for the mean. The default is a vector of ones. This matrix is extracted from the formulae in function CountsEPPM.
}
  \item{fixed.b}{
Set to the value of the parameter b if a fixed.b model is being used.
}
  \item{vnmax}{
A vector of the maximum counts for each vector in list.counts i.e. the list of grouped counts.
}
}
\value{
The list output with elements
\item{model}{The model being fitted}
\item{estimate}{Estimates of parameters}
\item{probabilities}{Vector of probabilities}
}
\author{
David M. Smith <dmccsmith@verizon.net>
}
\references{
Faddy M, Smith D. (2011). Analysis of count data with covariate dependence in both mean and variance. 
\emph{Journal of Applied Statistics}, \bold{38}, 2683-2694. \doi{10.1002/bimj.201100214}.
}
\examples{
all.counts=c(rep(0,5), 352, 479, 530, 291, 101, 17)
nmax1  <- length(all.counts)
nmax   <- nmax1 - 1
cnum   <- 0:nmax
ncount <- sum(all.counts)  
all.mean <- t(cnum) %*% all.counts / ncount 
alldata <- data.frame(all.mean)
mf <- model.frame(formula = all.mean ~ 1, data = alldata)
covariates.matrix.mean <- model.matrix(attr(mf, "terms"), data = mf)
list.counts <- list(all.counts)
parameter <- c(53.047752, -13.205655, 3.801599)
names(parameter) <- c('log(a)', 'c', 'log(b)') 
model.name <- 'Faddy distribution'
link <- "log"
attr(link, which="mean") <- make.link(link)
offset.mean <- matrix(c(rep(0, nrow(covariates.matrix.mean))), ncol=1)
output <- Model.Faddy(parameter, model.name, link,
 covariates.matrix.mean, offset.mean, fixed.b = NA,
 vnmax = c(10))
print(output)
}
\keyword{ models }
\name{Model.FaddyJMV.limiting}
\alias{Model.FaddyJMV.limiting}
\title{
Function to fit the limiting form of a Faddy distribution for under-dispersed counts.
}
\description{
Outputs probabilities for the limiting form of a Faddy distribution modeled by means and variances i.e. with the design matrices for mean and variance input together with data, offsets and scale factor model identifier.
}
\usage{
Model.FaddyJMV.limiting(parameter, covariates.matrix.mean, 
covariates.matrix.variance, offset.mean, offset.variance, 
scale.factor.model,vnmax)
}
\arguments{
  \item{parameter}{ 
A vector of the parameters of the model which is set to initial estimates on function call.
}
  \item{covariates.matrix.mean}{
A matrix of covariates for the mean where rows are the number of values in listcounts and columns the covariates. This matrix is extracted from the formulae in function CountsEPPM. However, in the accompanying example it is shown how it can be constructed 
independently of function CountsEPPM. 
}
  \item{covariates.matrix.variance}{
A matrix of covariates for the variance where rows are the number of values in listcounts and columns the covariates. The default is a vector of ones. This matrix is extracted from the formulae in function CountsEPPM. However, in the accompanying example it is shown how it can be constructed independently of function CountsEPPM. 
}
  \item{offset.mean}{
An offset vector for the mean. The default is a vector of ones.
}
  \item{offset.variance}{
An offset vector for the variance. The default is a vector of ones.
}
  \item{scale.factor.model}{
Takes one of the two values 'yes' or 'no' when 'yes' means that a model where the variance is expressed as a scale factor of the mean is desired. It only applies to 'mean and variance models' and the default is 'no'.
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
\references{
Faddy M, Smith D. (2011). Analysis of count data with covariate dependence in both mean and variance. 
\emph{Journal of Applied Statistics}, \bold{38}, 2683-2694. \doi{10.1002/bimj.201100214}.
}
\author{
David M. Smith <david.m.smith@truvenhealth.com>
}
\examples{
all.counts=c(rep(0,5),352,479,530,291,101,17)
nmax1  <- length(all.counts)
nmax   <- nmax1 - 1
cnum   <- 0:nmax
ncount <- sum(all.counts)  
all.mean <- sum(cnum*all.counts)/ncount 
all.variance <- (sum(cnum*cnum*all.counts) - ncount*all.mean*all.mean) / (ncount - 1)
alldata <- data.frame(all.mean,all.variance)
mf  <- model.frame(formula=all.mean~1,data=alldata)
covariates.matrix.mean <- model.matrix(attr(mf,"terms"),data=mf)
mf  <- model.frame(formula=all.variance~1,data=alldata)
covariates.matrix.variance <- model.matrix(attr(mf,"terms"),data=mf)
list.counts <- list(all.counts)
parameter <- c(1.8388023,0.6009881)
names(parameter) <- c("beta0 log(mean)","beta0 log(variance)") 
vnmax <- c(10)
offset.mean     <- matrix(c(rep(0,nrow(covariates.matrix.mean))),ncol=1)
offset.variance <- matrix(c(rep(0,nrow(covariates.matrix.mean))),ncol=1)
output <- Model.FaddyJMV.limiting(parameter,covariates.matrix.mean,
                    covariates.matrix.variance,offset.mean,offset.variance, 
                    scale.factor.model="no",vnmax)
print(output)
}
\keyword{ models }
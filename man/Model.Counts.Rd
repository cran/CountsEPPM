\name{Model.Counts}
\alias{Model.Counts}
\title{
Function for obtaining output from distributional models.
}
\description{
Produces output of model, parameters and probabilities from the various models.
}
\usage{
Model.Counts(parameter, model.type, model.name, link,
 covariates.matrix.mean, covariates.matrix.scalef,
 offset.mean, offset.scalef, fixed.b, vnmax)
}
\arguments{
  \item{parameter}{ 
A vector of the parameters of the model which is set to initial estimates on function call.
}
  \item{model.type}{
Takes one of two values i.e. 'mean only' or 'mean and scale-factor'. The 'mean only' value fits a linear predictor function to the parameter 'a' in equation (3) of Faddy and Smith (2011). If the model type being fitted is Poisson modeling 'a' is the same as modeling the mean. The 'mean and scale-factor' value fits linear predictor functions to both the mean and the scale-factor.
}
  \item{model.name}{
If model.type is 'mean only' the model being fitted is one of the three 'Poisson', 'negative binomial', 'Faddy distribution'. If model.type is 'mean and scale-factor' the model being fitted is either 'general' i.e. as equations (4) and (6) of Faddy and Smith (2011), or 'limiting' i.e. as equations (9) and (10) of Faddy and Smith (2011).
}
  \item{link}{ 
Takes one of one values i.e., 'log'. The default is 'log'. 
}
  \item{covariates.matrix.mean}{
A matrix of covariates for the mean where rows are the number of values in listcounts and columns the covariates. This matrix is extracted from the formulae in function CountsEPPM. However, in the accompanying example it is shown how it can be constructed 
independently of function CountsEPPM. 
}
  \item{covariates.matrix.scalef}{
A matrix of covariates for the scale-factor where rows are the number of values in listcounts and columns the covariates. The default is a vector of ones. This matrix is extracted from the formulae in function CountsEPPM. However, in the accompanying example it is shown how it can be constructed independently of function CountsEPPM. 
}
  \item{offset.mean}{
An offset vector for the mean. The default is a vector of ones.
}
  \item{offset.scalef}{
An offset vector for the scale-factor. The default is a vector of ones.
}
  \item{fixed.b}{
Set to the value of the parameter b if a fixed.b model is being used.
}
  \item{vnmax}{
A vector of the maximum counts for each vector in list.counts i.e. the list of grouped counts.
}
}
\value{
Output which is the output from either Model.Faddy, Model.Faddy.general, or Model.Faddy.limiting
}
\author{
David M. Smith <dmccsmith@verizon.net>
}
\references{
Faddy M, Smith D. (2011). Analysis of count data with covariate dependence in both mean and variance. 
\emph{Journal of Applied Statistics}, \bold{38}, 2683-2694. \doi{10.1002/bimj.201100214}.
}
\examples{
all.counts=c(rep(0,5),352,479,530,291,101,17)
nmax1  <- length(all.counts)
nmax   <- nmax1 - 1
cnum   <- 0:nmax
ncount <- sum(all.counts)  
all.mean <- sum(cnum*all.counts)/ncount 
all.scalef <- ((sum(cnum*cnum*all.counts) - ncount*all.mean*all.mean) / (ncount - 1)) / all.mean
alldata <- data.frame(all.mean, all.scalef)
mf  <- model.frame(formula = all.mean~1 ,data=alldata)
covariates.matrix.mean <- model.matrix(attr(mf,"terms"), data=mf)
mf  <- model.frame(formula = all.scalef~1, data = alldata)
covariates.matrix.scalef <- model.matrix(attr(mf,"terms"), data = mf)
list.counts <- list(all.counts)
parameter <- c(1.8388023, 0.6009881)
names(parameter) <- c("beta0 log(mean)" ,"beta0 log(scale-factor)") 
offset.mean     <- matrix(c(rep(0, nrow(covariates.matrix.mean))), ncol=1)
offset.scalef <- matrix(c(rep(0, nrow(covariates.matrix.mean))), ncol=1)
link <- "log"
attr(link, which="mean") <- make.link(link)
output <- Model.Counts(parameter, model.type = "mean and scale-factor",
 model.name = "limiting", link, covariates.matrix.mean,
 covariates.matrix.scalef, offset.mean, offset.scalef,
 fixed.b = NA, vnmax = c(10))
print(output)
}
\keyword{ models }
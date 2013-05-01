\name{Model.Faddy}
\alias{Model.Faddy}
\title{
Function for Faddy distribution with log link.
}
\description{
Returns probabilities for a Faddy distribution given inputs of model and parameters.
}
\usage{
Model.Faddy(parameter, model.type, model, covariates.matrix.mean, 
offset.mean, vnmax)
}
\arguments{
  \item{parameter}{ 
A vector of the parameters of the model which is set to initial estimates on function call.
}
  \item{model.type}{
Takes one of two values i.e. 'mean only' or 'mean and variance'. The 'mean only' value fits a linear predictor function to the parameter 'a' in equation (3) of Faddy and Smith (2011). If the model type being fitted is Poisson modeling 'a' is the same as modeling the mean. The 'mean and variance' value fits linear predictor functions to both the mean and the variance.
}
  \item{model}{
If model.type is 'mean only' the model being fitted is one of the three 'Poisson', 'negative binomial', 'Faddy distribution'. If model.type is 'mean and variance' the model being fitted is either 'general' i.e. as equations (4) and (6) of Faddy and Smith (2011), or 'limiting' i.e. as equations (9) and (10) of Faddy and Smith (2011).
}
  \item{covariates.matrix.mean}{
A matrix of covariates for the mean where rows are the number of values in listcounts and columns the covariates. This matrix is extracted from the formulae in function CountsEPPM. However, in the accompanying example it is shown how it can be constructed 
independently of function CountsEPPM. 
}
  \item{offset.mean}{
An offset vector for the mean. The default is a vector of ones. This matrix is extracted from the formulae in function CountsEPPM.
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
David M. Smith <david.m.smith@truvenhealth.com>
}
\references{
Faddy M, Smith D (2011). Analysis of count data with covariate dependence in both mean and variance. Journal of Applied Statistics, 38, 2683-2694.
}
\examples{
all.counts=c(rep(0,5),352,479,530,291,101,17)

nmax1  <- length(all.counts)
nmax   <- nmax1 - 1
cnum   <- 0:nmax
ncount <- sum(all.counts)  
all.mean <- sum(cnum*all.counts)/ncount 
alldata <- data.frame(all.mean)

model.type <- 'mean only'
model <- 'Faddy distribution'

mf  <- model.frame(formula=all.mean~1,data=alldata)
covariates.matrix.mean <- model.matrix(attr(mf,'terms'),data=mf)
listcounts <- list(all.counts)

parameter <- c(53.047752,3.801599,-13.205655)
names(parameter) <- c('log(a)','log(b)','c') 
vnmax <- c(10)
offset.mean <- matrix(c(rep(0,nrow(covariates.matrix.mean))),ncol=1)
output <- Model.Faddy(parameter,model.type,model,
                      covariates.matrix.mean,offset.mean,vnmax)
print(output)
}
\keyword{ models }

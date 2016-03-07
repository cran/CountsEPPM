\name{Model.Counts}
\alias{Model.Counts}
\title{
Function for obtaining output from distributional models.
}
\description{
Produces output of model, parameters and probabilities from the various models.
}
\usage{
Model.Counts(parameter, model.type, model, covariates.matrix.mean, 
covariates.matrix.variance,offset.mean, offset.variance, 
scale.factor.model, fixed.b, vnmax)
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
David M. Smith <david.m.smith@truvenhealth.com>
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
all.variance <- (sum(cnum*cnum*all.counts) - ncount*all.mean*all.mean) / (ncount - 1)
alldata <- data.frame(all.mean,all.variance)
mf  <- model.frame(formula=all.mean~1,data=alldata)
covariates.matrix.mean <- model.matrix(attr(mf,"terms"),data=mf)
mf  <- model.frame(formula=all.variance~1,data=alldata)
covariates.matrix.variance <- model.matrix(attr(mf,"terms"),data=mf)
list.counts <- list(all.counts)
parameter <- c(1.8388023,0.6009881)
names(parameter) <- c("beta0 log(mean)","beta0 log(variance)") 
fixed.b <- NA
vnmax <- c(10)
offset.mean     <- matrix(c(rep(0,nrow(covariates.matrix.mean))),ncol=1)
offset.variance <- matrix(c(rep(0,nrow(covariates.matrix.mean))),ncol=1)
output <- Model.Counts(parameter,
            model.type="mean and variance",model="limiting",
            covariates.matrix.mean,covariates.matrix.variance,
            offset.mean,offset.variance,
            scale.factor.model="no",fixed.b,vnmax)
print(output)
}
\keyword{ models }
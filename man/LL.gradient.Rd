\name{LL.gradient}
\alias{LL.gradient}
\title{
Function used to calculate the first derivatives of the log likelihood with respect to the model parameters.
}
\description{
Function used to calculate the first derivatives of the log likelihood with respect to the model parameters. 
These are numerical derivatives calculated using the numerical derivative functions of Gilbert and Varadhan (2015). 
}
\usage{
LL.gradient(parameter, model.type, model.name, link, list.data,
   covariates.matrix.mean, covariates.matrix.scalef, 
   offset.mean, offset.scalef, ltvalue, utvalue, fixed.b,
   weights, grad.method)
}
\arguments{
  \item{parameter}{ 
A vector of the parameters of the model which is set to initial estimates on function call.
}
  \item{model.type}{
Takes one of two values i.e. 'mean only' or 'mean and scale-factor'. The 'mean only' value fits linear predictor functions to the mean 
as in Faddy and Smith (2012). The 'mean and scale-factor' value fits linear predictor functions to both the 'mean' and the scale-factor. The default is 'mean and scale-factor'.
}
  \item{model.name}{
If model.type is 'mean only' the model being fitted is one of the four 'binomial', 'generalized binomial', 'beta binomial' or 'correlated binomial'. If model.type is 'mean and scale-factor' the model being fitted is one of the three 'generalized binomial', 'beta binomial' or 'correlated binomial'. Information about these models is given in Faddy and Smith (2012). The default is 'generalized binomial'.
}
  \item{link}{ 
Takes one of one values i.e., 'log'. The default is 'log'. 
}
  \item{list.data}{
A list of vectors of the counts as grouped data i.e. number of observations for each possible count value.
}
  \item{covariates.matrix.mean}{
A matrix of covariates for the mean where rows are the number of values in list.data and columns the covariates. This matrix is extracted from the formulae in function CountsEPPM. However, in the accompanying example it is shown how it can be constructed independently of function CountsEPPM. 
}
  \item{covariates.matrix.scalef}{
A matrix of covariates for the variance where rows are the number of 
values in list.binary and columns the covariates. The default is a 
vector of ones. This matrix is extracted from the formulae in function 
CountsEPPM. However, in the accompanying example it is shown how it 
can be constructed independently of function CountsEPPM. 
}
  \item{offset.mean}{
An offset vector for the probability of success p. The default 
is a vector of ones.
}
  \item{offset.scalef}{
An offset vector for the scale-factor. The default is a 
vector of ones.
}
   \item{ltvalue}{
Lower truncation value.
}
  \item{utvalue}{
Upper truncation value.
}
  \item{fixed.b}{
Set to the value of the parameter b if a fixed.b model is being used.
}
 \item{weights}{
A vector or list of weights for the modeling of probability of success. The default is a vector of ones.
}
\item{grad.method}{
Set to the method to be used to calculate the gradients either "simple" or "Richardson".
}
}
\value{
A vector of numerical first derivatives.
}
\references{
Gilbert P, Varadhan R. (2016). numDeriv: Accurate Numerical Derivatives. R Package version 2016.8-1.1, 
\url{https://CRAN.R-project.org/package=numDeriv}.
}
\author{
David M. Smith <dmccsmith@verizon.net>
}
\examples{
\dontrun{
gradient <- grad( LL.Regression.Binary, x = parameter,
 model.type = model.type, model.name = model.name, link = link,
 ntrials = ntrials, nsuccess = nsuccess,
 covariates.matrix.mean = covariates.matrix.mean, 
 covariates.matrix.scalef = covariates.matrix.scalef,
 offset.mean = offset.mean, offset.scalef = offset.scalef,
 weights = weights, grad.method = "Richardson")
 return(gradient)
  }
}
\keyword{ model }
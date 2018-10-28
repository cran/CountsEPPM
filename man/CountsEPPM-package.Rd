\name{CountsEPPM-package}
\alias{CountsEPPM-package}
\docType{package}
\title{
Fitting of EPPM models to count and binary data.
}
\description{
Fits regression models to under- and over-dispersed count data using extended Poisson process models.
}
\details{
\tabular{ll}{
Package: \tab CountsEPPM\cr
Type: \tab Package\cr
Version: \tab 2.1\cr
Date: \tab 2016-03-04\cr
License: \tab GPL-2 \cr
}
Using Generalized Linear Model (GLM) terminology, the functions utilize linear predictors for mean and variance with log link functions
to fit the regression models. Smith and Faddy (2016) gives further details about the package as well as examples of its use.
}
\author{
David M. Smith <smithdm1@us.ibm.com>
}
\references{
Faddy M, Smith D. (2011). Analysis of count data with covariate dependence in both mean and variance. 
\emph{Journal of Applied Statistics}, \bold{38}, 2683-2694. \doi{10.1002/bimj.201100214}.

Smith D, Faddy M. (2016). Mean and Variance Modeling of Under- and Overdispersed Count Data. 
\emph{Journal of Statistical Software}, \bold{69}(6), 1-23. \doi{10.18637/jss.v069.i06}.

Zeileis A, Croissant Y. (2010). Extended Model Formulas in R: Multiple Parts and Multiple Responses.
\emph{Journal of Statistical Software}, \bold{34}(1), 1-13. \doi{10.18637/jss.v034.i01}.
}
\keyword{ package }
\examples{
data(herons.group)
initial <- c(1.9871533,1.9900881,3.6841305,0.4925816)
names(initial) <- c("Adult mean","Immature mean", "Variance","log(b)")
output.fn  <- CountsEPPM(number.attempts~0+group | 1, herons.group,initial=initial)
print(output.fn)
}

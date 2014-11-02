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
Version: \tab 2.0\cr
Date: \tab 2014-09-07\cr
License: \tab GPL-2 \cr
}
Using Generalized Linear Model (GLM) terminology, the functions utilize linear predictors for mean and variance with log link functions
to fit the regression models. 
}
\author{
David M. Smith <david.m.smith@truvenhealth.com>
}
\references{
Faddy M, Smith D (2011). Analysis of count data with covariate dependence in both mean and variance. 
Journal of Applied Statistics, 38, 2683-2694.

Zeileis A, Croissant Y (2010). Extended Model Formulas in R: Multiple Parts and Multiple Responses.
Journal of Statistical Software, 34(XYZ), 1-12. URL http://www.jstatsoft.org/v34/iXYZ/.
}
\keyword{ package }
\examples{
data(herons.group)
initial <- c(1.9871533,1.9900881,3.6841305,0.4925816)
names(initial) <- c("Adult mean","Immature mean",
                    "Variance","log(b)")
output.fn  <- CountsEPPM(number.attempts~0+group | 1,
                         herons.group,initial=initial)
CountsEPPM.summary(output.fn)
}

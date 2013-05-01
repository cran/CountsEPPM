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
Version: \tab 1.0\cr
Date: \tab 2013-04-24\cr
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
data(ceriodaphnia.group)
initial        <- c(3.14041013,0.17347830,-0.01961849,
                    4.43874024,-0.49266496,0.02774201,-0.14063983 )
names(initial) <- c('(Intercept)','vdose','vdose2',
                    '(Intercept)','vdose','vdose2','log(beta)')
output.fn <- 
  CountsEPPM(mean.obs | variance.obs~1+vdose+vdose2 | 1+vdose+vdose2,
                         initial=initial,ceriodaphnia.group)
CountsEPPM.summary(output.fn)
}

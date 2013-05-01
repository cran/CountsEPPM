\name{CountsEPPM}
\alias{CountsEPPM}
\title{
Fitting of EPPM models to count data.
}
\description{
Fits regression models to under- and over-dispersed count data using extended Poisson process models.
}
\usage{
CountsEPPM(formula, data, model.type = 'mean and variance', model = 'general', 
initial = c(NA), ltvalue = NA, utvalue = NA, optimization.method = 'optim', 
scale.factor.model = 'no')
}
\arguments{
  \item{formula}{
Formulae for the mean and variance. The package 'Formula' of Zeileis and Croissant (2010) which allows multiple parts and multiple responses is used. Only formulae with a single response variable and single linear predictor, or those with two response variables and two linear predictors are allowed. The response variables must be named 'mean.obs', 'variance.obs' or 'scalef.obs' according to the model type being fitted. The values for these response variables are not input as part of 'data' they are calculated within the function from the list of grouped count data input. If the 'model.type' is 'mean only' only a single formula consisting of a response variable of 'mean.obs' and its related linear predictor is required. If the 'model.type' is 'mean and variance' and 'scale.factor.model'='no' two formula are needed with response variables 'mean.obs' and 'variance.obs' together with their related linear predictors. If 'scale.factor.model'='yes' the second of these response variables should be 'scalef.obs'.     
}
  \item{data}{
Data is an object with four components these being 'covariates', 'listcounts', 'offset.mean', 'offset.variance'. The component 'covariates' is a data frame of vectors of covariates. The component 'listcounts' is a list of vectors of the grouped counts. The components 'offset.mean', 'offset.variance' are each a vector of offsets for the mean and variance respectively. 
}
  \item{model.type}{
Takes one of two values i.e. 'mean only' or 'mean and variance'. The 'mean only' value fits a linear predictor function to the parameter 'a' in equation (3) of Faddy and Smith (2011). If the model type being fitted is Poisson modeling 'a' is the same as modeling the mean. For the negative binomial the mean is 'b'(exp('a')-1), 'b' also being as in equation (3) of Faddy and Smith (2011). The 'mean and variance' value fits linear predictor functions to both the mean and the variance. 
}
  \item{model}{
If model.type is 'mean only' the model being fitted is one of the three 'Poisson', 'negative binomial', 'Faddy distribution'. If model.type is 'mean and variance' the model being fitted is either 'general' i.e. as equations (4) and (6) of Faddy and Smith (2011), or 'limiting' i.e. as equations (9) and (10) of Faddy and Smith (2011).
}
  \item{initial}{
This is a vector of initial values for the parameters. If this vector is not present (NA) then initial values are calculated within the function.
}
  \item{ltvalue}{
Lower truncation value.
}
  \item{utvalue}{
Upper truncation value.
}
  \item{optimization.method}{
Takes one of the two values 'optim' or 'nlm' these being R optimization functions.
}
  \item{scale.factor.model}{
Takes one of the two values 'yes' or 'no' where 'yes' means that a model where the variance is expressed as a scale factor of the mean is desired. It only applies to 'mean and variance models' and the default is 'no'.
}
}
\value{
\item{model.type}{The type of model being fitted}
\item{model}{The model being fitted}
\item{covariates.matrix.mean}{The design matrix for the means}
\item{covariates.matrix.variance}{The design matrix for the variances}
\item{offset.mean}{The offset vector for the means}
\item{offset.variance}{The offset vector for the variances}
\item{ltvalue}{The lower truncation value}
\item{utvalue}{The upper truncation value}
\item{scale.factor.model}{Whether a scale factor model is desired no:yes}
\item{estimates}{Estimates of model parameters}
\item{vnmax}{Vector of maximums of grouped count data vectors in list.counts}
\item{loglikelihood}{Loglikelihood}
}
\references{
Faddy M, Smith D (2011). Analysis of count data with covariate dependence in both mean and variance. 
Journal of Applied Statistics, 38, 2683-2694.

Zeileis A, Croissant Y (2010). Extended Model Formulas in R: Multiple Parts and Multiple Responses.
Journal of Statistical Software, 34(XYZ), 1-12. URL http://www.jstatsoft.org/v34/iXYZ/.
}
\author{
David M. Smith <david.m.smith@truvenhealth.com>
}
\examples{
Luningetal.all <- list(covariates=data.frame(c(1)),
                       list.counts=list(c(rep(0,5),352,479,530,291,101,17)))  
model.type <- 'mean only'
model <- 'Faddy distribution'
output.fn <- CountsEPPM(mean.obs~1,Luningetal.all,model.type,model,
                         ltvalue=4,utvalue=11,optimization.method='nlm')
CountsEPPM.summary(output.fn)
initial <- output.fn$estses$estimates
names(initial) <- c('log(a)','log(b)','c') 
output.fn <- CountsEPPM(mean.obs~1,Luningetal.all,model.type,model,initial,
                         ltvalue=4,utvalue=11,optimization.method='nlm')
CountsEPPM.summary(output.fn)
}
\keyword{ models }

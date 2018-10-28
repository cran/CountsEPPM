\name{CountsEPPM}
\alias{CountsEPPM}
\title{
Fitting of EPPM models to count data.
}
\description{
Fits regression models to under- and over-dispersed count data using extended Poisson process models.
}
\usage{
CountsEPPM(formula, data, subset=NULL, na.action=NULL, weights=NULL,
model.type = "mean and scale-factor", model.name = "general", 
link="log", initial = NULL, ltvalue = NA, utvalue = NA, 
method = "Nelder-Mead", control = NULL, fixed.b = NA)
}
\arguments{
  \item{formula}{
Formulae for the mean and variance. The package 'Formula' of Zeileis and Croissant (2010) which allows multiple parts and multiple responses is used. 'formula' should consist of a left hand side (lhs) of single response variable and a right hand side (rhs) of one or two sets of variables for the linear predictors for the mean and (if two sets) the variance. This is as used for the R function 'glm' and also, for example, as for the package 'betareg' (Cribari-Neto and Zeileis, 2010). The function identifies from the argument data whether a data frame (as for use of 'glm') or a list (as required in Version 1.0 of this function) has been input. The list should be exactly the same as for a data frame except that the response variable is a list of vectors of frequency distributions rather than a vector of single counts as for the data frame. As with version 1.0 of this function, the subordinate functions fit models where the response variables are 'mean.obs', 'variance.obs' or 'scalef.obs' according to the model type being fitted. The values for these response variables are not input as part of 'data', they are calculated within the function from a list of grouped count data input. If the 'model.type' is 'mean only' 'formula' consists of a lhs of the response variable and and a rhs of the terms of the linear predictor for the mean model. If the 'model.type' is 'mean and variance' and 'scale.factor.model'='no' there are two set of terms in the rhs of 'formula' i.e., 'mean.obs' and 'variance.obs' together with the two sets of terms for the linear predictors of mean and variance. If 'scale.factor.model'='yes' the second response variable used by the subordinate functions would be 'scalef.obs'.     
}
  \item{data}{
'data' should be either a data frame (as for use of 'glm') or a list (as required in Version 1.0 of this function). The list should be exactly the same as for a data frame except that the response variable is a list of vectors of frequency distributions rather than a vector of single counts as for the data frame. Within the function a working list 'listcounts' and data frames with components such as 'mean.obs', 'variance.obs', 'scalef.obs', 'covariates', 'offset.mean', 'offset.variance' are set up . The component 'covariates' is a data frame of vectors of covariates in the model. The component 'listcounts' is a list of vectors of the grouped counts, or the single counts in grouped form if 'data' is a data frame. 
}
  \item{subset}{
Subsetting commands.
}
  \item{na.action}{
Action taken for NAs in data.
}
  \item{weights}{
Vector of list of lists of weights.
}
  \item{model.type}{
Takes one of two values i.e. 'mean only' or 'mean and variance'. The 'mean only' value fits a linear predictor function to the parameter 'a' in equation (3) of Faddy and Smith (2011). If the model type being fitted is Poisson modeling 'a' is the same as modeling the mean. For the negative binomial the mean is 'b'(exp('a')-1), 'b' also being as in equation (3) of Faddy and Smith (2011). The 'mean and variance' value fits linear predictor functions to both the mean and the variance. 
}
  \item{model.name}{
If model.type is 'mean only' the model being fitted is one of the three 'Poisson', 'negative binomial', 'Faddy distribution'. If model.type is 'mean and scale-factor' the model being fitted is either 'general' i.e. as equations (4) and (6) of Faddy and Smith (2011), or 'limiting' i.e. as equations (9) and (10) of Faddy and Smith (2011).
}
  \item{link}{
Takes one of one values i.e., 'log'. The default is 'log'. 
}
  \item{initial}{
This is a vector of initial values for the parameters. If this vector is NULL then initial values based on a fitting Poisson models using 'glm' are calculated within the function.
}
  \item{ltvalue}{
Lower truncation value.
}
  \item{utvalue}{
Upper truncation value.
}
  \item{method}{
Optimization method takes one of the two values 'Nelder-Mead' or 'BFGS' these being options for the \code{optim} function.
}
  \item{control}{
'control' is a list of control parameters as used in 'optim' or 'nlm'. If this list is NULL the defaults for 'optim' are set as 'control   <- list(fnscale=-1,trace=0,maxit=1000)' and for 'nlm' are set as 'control <- list(fscale=1,print.level=0,stepmax=1,gradtol=1e-8,steptol=1e-10,iterlim=500)'. For 'optim' the control parameters that can be changed by inputting a variable length list are 'fnscale, trace, maxit, abstol, reltol, alpha, beta, gamma'. For 'nlm' the parameters are 'fscale, print.level, stepmax, gradtol,steptol, iterlim'. Details of 'optim' and 'nlm' and their control parameters are available in the online R help manuals.
}
  \item{fixed.b}{
Set to the value of the parameter b if a fixed.b model is being used.
}
}
\details{
Smith and Faddy (2016) gives further details as well as examples of use.
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
\item{estimates}{Estimates of model parameters}
\item{vnmax}{Vector of maximums of grouped count data vectors in list.counts}
\item{loglikelihood}{Loglikelihood}
}
\references{
Cribari-Neto F, Zeileis A. (2010). Beta Regression in R.
\emph{Journal of Statistical Software}, \bold{34}(2), 1-24. \doi{10.18637/jss.v034.i02}.

Grun B, Kosmidis I, Zeileis A. (2012). Extended Beta Regression in R: Shaken, Stirred, Mixed, and Partitioned.
\emph{Journal of Statistical Software}, \bold{48}(11), 1-25. \doi{10.18637/jss.v048.i11}.

Faddy M, Smith D. (2011). Analysis of count data with covariate dependence in both mean and variance. 
\emph{Journal of Applied Statistics}, \bold{38}, 2683-2694. \doi{10.1002/bimj.201100214}.

Smith D, Faddy M. (2016). Mean and Variance Modeling of Under- and Overdispersed Count Data. 
\emph{Journal of Statistical Software}, \bold{69}(6), 1-23. \doi{10.18637/jss.v069.i06}.

Zeileis A, Croissant Y. (2010). Extended Model Formulas in R: Multiple Parts and Multiple Responses.
\emph{Journal of Statistical Software}, \bold{34}(1), 1-13. \doi{10.18637/jss.v034.i01}.
}
\author{
David M. Smith <smithdm1@us.ibm.com>
}
\examples{
data(herons.group)
initial <- c(0.5623042, 0.4758576, 0.5082486)
names(initial) <- c("Adult mean", "Immature mean", "log(b)")
output.fn <- CountsEPPM(number.attempts ~ 0 + group,
 herons.group, model.type = 'mean only', model = 'negative binomial',
 initial = initial)
print(output.fn)
}
\keyword{ models }
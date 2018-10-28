\name{predict.CountsEPPM}
\alias{predict.CountsEPPM}
\title{
Prediction Method for CountsEPPM Objects
}
\description{
Extract various types of predictions from CountsEPPM regression models.
}
\usage{
\method{predict}{CountsEPPM}(object, newdata = NULL,
 type = c("response", "linear.predictor.mean",
 "linear.predictor.scale.factor", "scale.factor",
 "mean", "variance", "distribution", "distribution.parameters"),
 na.action = na.pass, \dots)
}
\arguments{
  \item{object}{
fitted model object of class "CountsEPPM".
}
  \item{newdata}{
optionally, a data frame in which to look for variables with which to predict. If omitted, the original observations are used.
}
  \item{type}{
character indicating type of predictions: fitted means of responses
("response"), linear predictors ("linear.predictor.mean", "linear.predictor.scale.factor"), fitted value of mean ("mean"),
fitted value of scale-factor ("scale.factor"), fitted value of variance ("variance"), fitted probability distribution ("distribution"), parameters of fitted distributions ("distribution.parameters")
}
  \item{na.action}{
function determining what to do with missing values in \emph{newdata}. The default is to predict NA.
}
  \item{\dots}{
some methods for this generic function require additional arguments.
}
}
\value{
A vector or list of the predicted values from the fitted model object.
}
\references{
Cribari-Neto F, Zeileis A. (2010). Beta Regression in R.
\emph{Journal of Statistical Software}, \bold{34}(2), 1-24. \doi{10.18637/jss.v034.i02}.
}
\author{
David M. Smith <smithdm1@us.ibm.com>
}
\seealso{
\code{\link[betareg]{predict.betareg}} 
}
\examples{
data("herons.group") 
output.fn  <- CountsEPPM(number.attempts ~ 0 + group,
 herons.group, model.type = 'mean only', model.name = 'Poisson')
predict(output.fn, type = "response")
predict(output.fn, type = "linear.predictor.mean")
}
\keyword{ methods }
\keyword{ distribution }

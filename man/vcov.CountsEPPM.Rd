\name{vcov.CountsEPPM}
\alias{vcov.CountsEPPM}
\title{
Variance/Covariance Matrix for Coefficients
}
\description{
Variance/covariance matrix for coefficients of fitted model.
}
\usage{
\method{vcov}{CountsEPPM}(object, model = c("full", "mean", "scale.factor"), \dots)
}
\arguments{
  \item{object}{
fitted model object of class "CountsEPPM".
}
  \item{model}{
character indicating variance/covariance matrix for all coefficients to be output: all coefficients ("full"),
variance/covariance matrix for coefficients of probability of success ("mean"), variance/covariance matrix for coefficients of scale-factor ("scale.factor")
}
  \item{\dots}{
}
}
\value{
The variance/covariance matrix of the parameters of the fitted model object.
}
\references{
Cribari-Neto F, Zeileis A. (2010). Beta Regression in R.
\emph{Journal of Statistical Software}, \bold{34}(2), 1-24. \doi{10.18637/jss.v034.i02}.
}
\author{
David M. Smith <smithdm1@us.ibm.com>
}
\seealso{
\code{\link[betareg]{betareg}} 
}
\examples{
data("herons.group") 
output.fn  <- CountsEPPM(number.attempts ~ 0 + group,
 herons.group, model.type = 'mean only', model.name = 'Poisson')
vcov(output.fn)
}
\keyword{ models }

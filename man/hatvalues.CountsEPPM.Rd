\name{hatvalues.CountsEPPM}
\alias{hatvalues.CountsEPPM}
\title{
Extraction of hat matrix values from CountsEPPM Objects
}
\description{
Extract the values of the hat matrix from models of class "CountsEPMM".
}
\usage{
\method{hatvalues}{CountsEPPM}(model, \dots)
}
\arguments{
  \item{model}{
fitted model object of class "CountsEPPM".
}
  \item{\dots}{
some methods for this generic function require additional arguments.
}
}
\value{
The calculated hat values for the fitted model. These are used to calculate Cook's distances.
}
\references{
Cribari-Neto F, Zeileis A. (2010). Beta Regression in R.
\emph{Journal of Statistical Software}, \bold{34}(2), 1-24. \doi{10.18637/jss.v034.i02}.
}
\author{
David M. Smith <dmccsmith@verizon.net>
}
\seealso{
\code{\link[betareg]{betareg}} 
}
\examples{
data("herons.group") 
output.fn  <- CountsEPPM(number.attempts ~ 0 + group,
 herons.group, model.type = 'mean only', model.name = 'Poisson')
cooks.distance(output.fn)
hatvalues(output.fn)
}
\keyword{ methods }

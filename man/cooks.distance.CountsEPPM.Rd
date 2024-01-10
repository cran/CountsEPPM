\name{cooks.distance.CountsEPPM}
\alias{cooks.distance.CountsEPPM}
\title{
Cook's distance for CountsEPPM Objects
}
\description{
Calculates Cook's distances for CountsEPPM objects.
}
\usage{
\method{cooks.distance}{CountsEPPM}(model, \dots)
}
\arguments{
  \item{model}{
fitted model object of class "CountsEPPM".
}
  \item{\dots}{
some methods for this generic function require additional arguments.
}
}
\details{
Cook's distances as in GLMs.
}
\value{
A vector of Cook's distances.
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
}
\keyword{ methods }

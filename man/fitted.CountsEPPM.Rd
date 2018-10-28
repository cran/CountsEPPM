\name{fitted.CountsEPPM}
\alias{fitted.CountsEPPM}
\title{
Extraction of fitted values from CountsEPPM Objects
}
\description{
This function is generic. Extract the fitted values from models of class "BinaryEPMM".
}
\usage{
\method{fitted}{CountsEPPM}(object, \dots)
}
\arguments{
  \item{object}{
fitted model object of class "CountsEPPM".
}
  \item{\dots}{
currently not used.
}
}
\seealso{
\code{\link[stats:stats-package]{fitted}} 
}
\author{
David M. Smith <smithdm1@us.ibm.com>
}
\examples{
data("herons.group") 
output.fn  <- CountsEPPM(number.attempts ~ 0 + group,
 herons.group, model.type = 'mean only', model.name = 'Poisson')
cooks.distance(output.fn)
fitted(output.fn)
}
\keyword{ methods }

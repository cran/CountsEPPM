\name{coef.CountsEPPM}
\alias{coef.CountsEPPM}
\title{
Extraction of model coefficients for CountsEPPM Objects
}
\description{
Extract the regression model coefficients from models of class "BinaryEPMM".
}
\usage{
\method{coef}{CountsEPPM}(object, prtpar = c("full", "mean", "scale.factor"), \dots)
}
\arguments{
  \item{object}{
fitted model object of class "CountsEPPM".
}
  \item{prtpar}{
character indicating coefficients of the fitted model to be output:
all coefficients ("full"),
coefficients of the model for probability of success ("mean"),
coefficients of the model for scale-factor ("scale.factor")
}
  \item{\dots}{
some methods for this generic function require additional arguments.
}
}
\details{
One of a set of standard extractor functions for fitted model objects of class
"CountsEPPM.
}
\value{
Vector of coefficients of fitted regression model.
}
\author{
David M. Smith <dmccsmith@verizon.net>
}
\seealso{
\code{\link[betareg]{betareg}} 
}
\examples{
data(herons.group)
output.fn  <- CountsEPPM(number.attempts ~ 0 + group,
 herons.group, model.type = 'mean only', model.name = 'Poisson')
coef(output.fn, prtpar = "full")
coef(output.fn, prtpar = "mean")
coef(output.fn, prtpar = "scale.factor")
}
\keyword{ methods }

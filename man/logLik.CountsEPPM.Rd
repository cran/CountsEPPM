\name{logLik.CountsEPPM}
\alias{logLik.CountsEPPM}
\title{
Method for CountsEPPM object
}
\description{
This function is generic and enables the use of functions related to the model fitting involved with \code{lm} and \code{glm} objects such as \code{AIC}.
}
\usage{
\method{logLik}{CountsEPPM}(object, \dots)
}
\arguments{
  \item{object}{
The object output from \code{CountsEPPM}. 
}
  \item{\dots}{
currently not used. 
}
}
\author{
David M. Smith <david.m.smith@truvenhealth.com>
}
\examples{
\dontrun{
output.fn <- CountsEPPM(mean.obs~1,Luningetal.all,model.type,model,initial,
                         ltvalue=4,utvalue=11,optimization.method="nlm")
logLik.CountsEPPM(object=output.fn)
}
}
\keyword{ Methods }

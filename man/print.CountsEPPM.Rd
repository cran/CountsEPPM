\name{print.CountsEPPM}
\alias{print.CountsEPPM}
\title{
Printing of CountsEPPM Objects
}
\description{
Prints objects of class "CountsEPPM".
}
\usage{
\method{print}{CountsEPPM}(x, digits = max(3, getOption("digits") - 3), \dots)
}
\arguments{
  \item{x}{
fitted model object of class "CountsEPPM".
}
  \item{digits}{
digits of printed output.
}
  \item{\dots}{
not currently used.
}
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
print(output.fn)
}
\keyword{ IO }

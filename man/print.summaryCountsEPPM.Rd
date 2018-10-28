\name{print.summaryCountsEPPM}
\alias{print.summaryCountsEPPM}
\title{
Printing of summaryCountsEPPM Objects 
}
\description{
Prints the objects of class "summaryCountsEPPM".
}
\usage{
\method{print}{summaryCountsEPPM}(x, \dots)
}
\arguments{
  \item{x}{
object output by \code{summary.CountsEPPM}.
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
print(summary(output.fn))
}
\keyword{ IO }

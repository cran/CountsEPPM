\name{residuals.CountsEPPM}
\alias{residuals.CountsEPPM}
\title{
Residuals for CountsEPPM Objects
}
\description{
This function is generic. Extract various types of residuals from objects of class "CountsEPPM".
}
\usage{
\method{residuals}{CountsEPPM}(object, type = c("spearson",
 "deviance", "pearson", "response", "likelihood", "sdeviance"),
 \dots)
}
\arguments{
  \item{object}{
Fitted model object of class "CountsEPPM".
}
  \item{type}{
Type of residuals wanted i.e., standardized Pearson "spearson", deviance "deviance", Pearson "pearson", response "response", likelihood "likelihood", standardized deviance "sdeviance".
}
  \item{\dots}{
some methods for this geneic function require additional arguments.
}
}
\details{
Residuals as Cribari-Neto and Zeileis (2010).
}
\seealso{
\code{\link[betareg]{residuals.betareg}} 
}
\references{
Cribari-Neto F, Zeileis A. (2010). Beta Regression in R.
\emph{Journal of Statistical Software}, \bold{34}(2), 1-24. \doi{10.18637/jss.v034.i02}.
}
\author{
David M. Smith <smithdm1@us.ibm.com>
}
\keyword{ models }

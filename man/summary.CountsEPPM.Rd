\name{summary.CountsEPPM}
\alias{summary.CountsEPPM}
\title{
Method for CountsEPPM object
}
\description{
This function is generic and is for printing out a summary of the results of fitting EPPM models to count data.
}
\usage{
\method{summary}{CountsEPPM}(object, \dots)
}
\arguments{
  \item{object}{
The object output from \code{CountsEPPM}. This list includes a vector vnmax of the maximums of the grouped count vectors in list.counts. The vector vnmax can be changed before calling this function in order to give more complete probability vectors i.e.
closer to a total of 1.
}
  \item{\dots}{
currently not used. 
}
}
\author{
David M. Smith <smithdm1@us.ibm.com>
}
\examples{
\dontrun{
output.fn <- CountsEPPM(mean.obs ~ 1, Luningetal.all, model.type,
 model,initial, ltvalue = 4, utvalue = 11, optimization.method =  "nlm")
summary(object=output.fn)
}
}
\keyword{ Methods }

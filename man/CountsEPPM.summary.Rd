\name{CountsEPPM.summary}
\alias{CountsEPPM.summary}
\title{
Print out summary results.
}
\description{
Print out summary of the results of fitting EPPM models to count data.
}
\usage{
CountsEPPM.summary(output.fn)
}
\arguments{
  \item{output.fn}{
The list of objects output from CountsEPPM. This list includes a vector vnmax of the
maximums of the grouped count vectors in list.counts. The vector vnmax can be changed 
before calling this function in order to give more complete probability vectors i.e.
closer to a total of 1.
}
}
\author{
David M. Smith <david.m.smith@truvenhealth.com>
}
\examples{
\dontrun{
output.fn <- CountsEPPM(mean.obs~1,Luningetal.all,model.type,model,initial,
                         ltvalue=4,utvalue=11,optimization.method="nlm")
CountsEPPM.summary(output.fn)
}
}
\keyword{ print }

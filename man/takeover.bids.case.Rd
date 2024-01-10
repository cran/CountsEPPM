\name{takeover.bids.case}
\alias{takeover.bids.case}
\docType{data}
\title{
Takeover bids data.
}
\description{
Data of the number of bids received by 126 U.S. firms that were targets of tender offers from 1978 to 1985.
}
\usage{data("takeover.bids.case")}
\format{
  A data frame with 126 observations on the following 12 variables.
  \describe{
    \item{\code{DOCNO}}{a numeric vector}
    \item{\code{WEEKS}}{a numeric vector}
    \item{\code{NUMBIDS}}{a numeric vector}
    \item{\code{BIDPREM}}{a numeric vector}
    \item{\code{INSTHOLD}}{a numeric vector}
    \item{\code{SIZE}}{a numeric vector}
    \item{\code{LEGLREST}}{a numeric vector}
    \item{\code{REALREST}}{a numeric vector}
    \item{\code{FINREST}}{a numeric vector}
    \item{\code{REGULATN}}{a numeric vector}
    \item{\code{WHTKNGHT}}{a numeric vector}
    \item{\code{SIZESQ}}{a numeric vector}
  }
}
\details{
Data originally from Jaggia and Thosar (1993) and used as an example in Cameron and Trivedi (2013) and Saez-Castillo and Conde-Sanchez (2013). 
}
\source{
Stata data file obtained from A.C. Cameron's webpage 
\url{https://cameron.econ.ucdavis.edu/}.
}
\references{
Cameron, A.C., Trivedi, P.K. (2013). Regression Analysis of Count Data. Cambridge University Press, second edition.

Jaggia, S., Thosar, S. (1993). Multiple Bids as a Consequence of Target Management Resistance.
\emph{Review of Quantitative Finance and Accounting}, 447-457.

Saez-Castillo, A.J., Conde-Sanchez, A. (2013). A hyper-Poisson regression model for overdispersed and underdispersed count data.
\emph{Computational Statistics and Data Analysis}, \bold{61}, 148-157.\doi{10.1016/j.csda.2012.12.009}
}
\examples{
data(takeover.bids.case)
print(takeover.bids.case)
}
\keyword{ datasets }

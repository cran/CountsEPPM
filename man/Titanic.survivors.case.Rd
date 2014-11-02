\name{Titanic.survivors.case}
\alias{Titanic.survivors.case}
\docType{data}
\title{
Titanic survivors data
}
\description{
These data are from the survival log of the Titanic and consist of the number of survivors out of the number of passengers broken down into age, sex and class categories.
}
\usage{data(Titanic.survivors.case)}
\format{
  A data frame with 12 observations on the following 5 variables.
  \describe{
    \item{\code{age}}{a factor with levels \code{ child} \code{ adult}}
    \item{\code{sex}}{a factor with levels \code{ female} \code{ male}}
    \item{\code{class}}{a factor with levels \code{ 1st class} \code{ 2nd class} \code{ 3rd class}}
    \item{\code{cases}}{a numeric vector}
    \item{\code{survive}}{a numeric vector}
  }
}
\details{
Hilbe (2011) first models these data as a logistic model, then finding that they are overdispersed, modeling them as count data (number of survivors, survive) with offset (log of the number of passengers, cases).
}
\source{
Section 9.5, Example 3, pages 263-268, Hilbe, J. (2011). 
}
\references{
Hilbe, J. (2011). Negative Binomial Regression. Cambridge University Press, second edition.
}
\examples{
data(Titanic.survivors.case)
print(Titanic.survivors.case)
}
\keyword{datasets}

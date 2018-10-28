\name{herons.case}
\alias{herons.case}
\docType{data}
\title{
Green-backed herons as two groups
}
\description{
The data are the numbers of attempts at foraging by 20 adult and 20 immature green-backed herons. The data are listed as grouped (adult or immature) count data i.e. number of herons having a particular count value.
}
\usage{data("herons.case")}
\format{
  A data frame with 40 observations on the following 2 variables.
  \describe{
    \item{\code{group}}{a factor with levels \code{ Adult} \code{ Immature}}
    \item{\code{number.attempts}}{a numeric vector}
  }
}
\source{
Zhu J, Eickhoff J, Kaiser M (2003). Modelling the Dependence between Number of Trials and Success Probability in Beta-Binomial-Poisson Mixture Distributions. 
\emph{Biometrics}, \bold{59}, 955-961.
\doi{10.1111/j.0006-341X.2003.00110.x}.
}
\examples{
data(herons.case)
print(herons.case)
}
\keyword{datasets}

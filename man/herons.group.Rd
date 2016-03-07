\name{herons.group}
\alias{herons.group}
\docType{data}
\title{
Green-backed herons as two groups
}
\description{
The data are the numbers of attempts at foraging by 20 adult and 20 immature green-backed herons. The data are listed as grouped (adult or immature) count data i.e. number of herons having a particular count value.
}
\usage{data(herons.group)}
\format{
  The format is:
List of 2
 $ group          : Factor w/ 2 levels " Adult"," Immature": 1 2
 $ number.attempts:List of 2
  ..$ : num [1:25] 0 5 2 1 1 1 0 2 0 1 ...
  ..$ : num [1:26] 0 2 2 1 5 1 2 2 1 1 ...
}
\source{
Zhu J, Eickhoff J, Kaiser M (2003). Modelling the Dependence between Number of Trials and Success Probability in Beta-Binomial-Poisson Mixture Distributions. 
\emph{Biometrics}, \bold{59}, 955-961. \doi{10.1111/j.0006-341x.2003.00110.x}.
}
\examples{
data(herons.group)
print(herons.group)
}
\keyword{datasets}

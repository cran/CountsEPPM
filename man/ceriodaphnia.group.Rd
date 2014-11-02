\name{ceriodaphnia.group}
\alias{ceriodaphnia.group}
\docType{data}
\title{
Ceriodaphnia data
}
\description{
Ceriodaphnia dubia are water fleas used to test the impact of effluents on water quality. The data are counts of young at varying effluent concentrations.
}
\usage{data(ceriodaphnia.group)}
\format{
  The format is:
List of 5
 $ vdose       : num [1:5] 0 1.56 3.12 6.25 12.5
 $ vdose2      : num [1:5] 0 2.44 9.77 39.06 156.25
 $ fdose       : Factor w/ 5 levels "0","1.5625","3.125",..: 1 2 3 4 5
 $ twofdose    : Factor w/ 2 levels "0","1": 1 1 2 2 1
 $ number.young:List of 5
  ..$ : num [1:32] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ : num [1:36] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ : num [1:45] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ : num [1:37] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ : num [1:17] 0 0 1 0 0 0 0 3 0 0 ...
}
\details{
The data is used in Faddy and Smith (2011) as an example. Faddy and Smith (2011) is the main reference for the methods implemented. The data are grouped into number of fleas for each count value.
}
\source{
Bailer, A., Oris, J. (1997). Estimating inhibition concentrations for different response scales using Generalized Linear Models. Environmental Toxicology and Chemistry, 16, 1554-1559.
}
\references{
Faddy, M., Smith, D. (2011). Analysis of count data with covariate dependence in both mean and variance. Journal of Applied Statistics, 38, 2683-2694.
}
\examples{
data(ceriodaphnia.group)
print(ceriodaphnia.group)
}
\keyword{datasets}

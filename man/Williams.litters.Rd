\name{Williams.litters}
\alias{Williams.litters}
\docType{data}
\title{
Number of trials (implantations) of data of Williams (1996).
}
\description{
The data is arranged as a list of grouped counts where the grouping is by dose where dose is included both as a variate (vdose) and as a factor (fdose).
}
\usage{data(Williams.litters)}
\format{
  The format is:
List of 2
 $ covariates :'data.frame':	4 obs. of  2 variables:
  ..$ vdose: num [1:4] 0 0.75 1.5 3
  ..$ fdose: Factor w/ 4 levels "0","0.75","1.5",..: 1 2 3 4
 $ list.counts:List of 4
  ..$ : num [1:18] 0 0 0 0 0 0 0 0 0 2 ...
  ..$ : num [1:20] 0 0 0 0 0 0 0 1 0 1 ...
  ..$ : num [1:17] 0 0 0 0 0 0 0 1 1 2 ...
  ..$ : num [1:17] 0 0 0 0 0 0 0 0 0 3 ...
}
\source{
Williams D (1996). Overdispersion in logistic linear model. In B Morgan (ed.), Statistics in Toxicology, pp. 75-84, Oxford Science Publications.
}
\examples{
data(Williams.litters)
print(Williams.litters)
}
\keyword{datasets}

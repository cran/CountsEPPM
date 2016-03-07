\name{Luningetal.litters}
\alias{Luningetal.litters}
\docType{data}
\title{
Number of trials (implantations) in data of Luning, et al. (1966)
}
\description{
The data are arranged as a list of grouped counts where the grouping is by dose where dose is included both as a variate (vdose) and as a factor (fdose).
}
\usage{data(Luningetal.litters)}
\format{
  The format is:
List of 3
 $ vdose        : num [1:3] 0 300 600
 $ fdose        : Factor w/ 3 levels "0","300","600": 1 2 3
 $ number.trials:List of 3
  ..$ : num [1:11] 0 0 0 0 0 71 156 224 150 70 ...
  ..$ : num [1:11] 0 0 0 0 0 121 170 186 99 24 ...
  ..$ : num [1:11] 0 0 0 0 0 160 153 120 45 7 ...
}
\source{
Luning K, Sheridan W, Ytterborn K, Gullberg U (1966). The relationship between the number of implantations and the rate of intra-uterine death in mice. 
\emph{Mutation Research}, \bold{3}, 444-451. \doi{10.1016/0027-5107(66)90054-6}.
}
\examples{
data(Luningetal.litters)
print(Luningetal.litters)
}
\keyword{datasets}

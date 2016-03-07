\name{distribution.CountsEPPM}
\alias{distribution.CountsEPPM}
\title{
Produce summary distributional information for fitted model.
}
\description{
Constructs a data frame of summary distributional information i.e. means, variances and
(optionally) probabilities from the fitted model information output from CountsEPPM.
}
\usage{
distribution.CountsEPPM(output.fn, output.probabilities = "no", 
output.FDparameters = "no")
}
\arguments{
  \item{output.fn}{
This is the object output by CountsEPPM.
}
  \item{output.probabilities}{
Text input 'no' or 'yes' stating whether or not the vectors of probabilities for the fitted 
model output of output.fn should be included.
}
  \item{output.FDparameters}{
Text input 'no' or 'yes' stating whether or not the vectors of parameters of the underlying 
Faddy distributions for the model output of output.fn should be included.
}
}
\details{
The object output by this function is a data frame consisting of means and variances for the 
fitted model of output.fn together (optionally) with related vectors of probabilities and/or 
the vectors of parameters of the Underlying Faddy distributions. If the model.type is 
'mean only' only the means and variances calculated from the vectors of 
probabilities are output, their variable names being mean.prob and variance.prob respectively. 
If the model.type is 'mean and variance' then the means and variances calculated from estimates 
of the parameters and their associated design matrices are also output, in which case their variable 
names are mean.par and variance.par. It is a feature of this package that the vectors 
of grouped counts in list.counts (see CountsEPPM) need only be input up to the maximum value  
with a non zero count as these are all that is needed to calculate the log likelihoods. 
However, the probability distributions could have non zero probability masses at much larger
values. In order to check the coverages of these distributions, hence also the accuracy 
of mean.prob and variance.prob, the total probabilities are also calculated and included in the
output data frame with the means and variances. If there is a desire for more complete coverages
the values in the vector output.fn\$vnmax can be increased and this function re-run. This process
can be continued until the coverages are considered adequate.
}
\value{
If the model.type is 'mean only' a data frame consisting of the vectors mean.prob, variance.prob
and totalprob is returned. If the model.type is 'mean and variance' these three vectors are
returned in the data frame together with two more mean.par and variance.par. If   
output.probabilities is set to 'yes' a list of probability vectors is added to the data frame.
}
\author{
David M. Smith <david.m.smith@truvenhealth.com>
}
\examples{
\dontrun{
output.distribution <- distribution.CountsEPPM(output.fn)
output.fn$vnmax <- output.fn$vnmax + 60
output.distribution <- distribution.CountsEPPM(output.fn,output.probabilities="yes",
                                               output.FDparameters = "yes")
}
}
\keyword{ distribution }

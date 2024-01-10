\name{plot.CountsEPPM}
\alias{plot.CountsEPPM}
\title{
Diagnostic Plots for CountsEPPM Objects
}
\description{
Various types of standard diagnostic plots can be produced, involving various 
types of residuals, influence measures etc.
}
\usage{
\method{plot}{CountsEPPM}(x, which = 1:4, caption = c("Residuals vs indices of obs.", "Cook's distance plot",
 "Leverage vs predicted values", "Residuals vs linear predictor",
 "Normal Q-Q plot of residuals", "Predicted vs observed values"),
 sub.caption = " ", main = "", ask = prod(par("mfcol"), 1) <
 length(which) && dev.interactive(), \dots, type = "spearson")
}
\arguments{
  \item{x}{
fitted model object of class "CountsEPPM".
}
  \item{which}{
numeric. If a subset of plots is required, specify a subset of the numbers 1:6.
}
  \item{caption}{
character. Captions to appear above the plots.
}
  \item{sub.caption}{
character. Common title-above figures if there are multiple.
}
  \item{main}{
chaacter. Title to each plot in addition to the above caption.
}
  \item{ask}{
logical. If true, the user is asked before each plot.
}
  \item{\dots}{
other parameters to be passed through to plotting functions.
}
  \item{type}{
character indicating type of residual to be used, see residuals.CountsEPPM.
}
}
\details{
The plot method for CountsEPPM objects produces various plots of diagnostic plots similar to those produced by \pkg{betareg}. See Ferrari and Cribari-Neto (2004) for further details of the displays of \pkg{betareg}.
}
\references{
Cribari-Neto F, Zeileis A. (2010). Beta Regression in R.
\emph{Journal of Statistical Software}, \bold{34}(2), 1-24. \doi{10.18637/jss.v034.i02}.
}
\author{
David M. Smith <dmccsmith@verizon.net>
}
\seealso{
\code{\link[betareg]{plot.betareg}} 
}
\examples{
\dontrun{
data("herons.group") 
output.fn  <- CountsEPPM(number.attempts ~ 0 + group,
 herons.group, model.type = 'mean only', model.name = 'Poisson')
cooks.distance(output.fn)
plot(output.fn, which = 1, type= "sdeviance")
}
}
\keyword{ hplot }

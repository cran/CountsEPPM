coef.CountsEPPM <-
function(object, prtpar = c("full", "mean", "scale.factor"), ...) {
      if (missing(prtpar)) { prtpar <- c("full") } 
# Checking for correct prtpar option
      if ((prtpar!="full") & (prtpar!="mean") & (prtpar!="scale.factor")) {
         cat("\n","unknown prtpar option","\n")
         coefficients <- NULL
                                   } else {
         if (prtpar=="full") { 
            if (object$model.name=="Poisson") {
               coefficients <- object$coefficients$mean.est
                                        } else {
               coefficients <- c(object$coefficients$mean.est,
                                 object$coefficients$scalef.est)
                                               } # end of if Poisson
             } else { npar.mean <- ncol(object$covariates.matrix.mean)
               if (prtpar=="mean") { coefficients <- object$coefficients$mean.est
                                }
               if (prtpar=="scale.factor") { 
                  if (is.null(object$coefficients$scalef.est)==TRUE) {
                     coefficients <- NULL
                                           } else {
                     coefficients <- object$coefficients$scalef.est } } # end of if scale-factor
                                    } } # end of if prtpar!= full, etc.
      return(coefficients) }

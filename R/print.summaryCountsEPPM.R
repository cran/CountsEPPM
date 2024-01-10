print.summaryCountsEPPM <-
function(x, ...) {

       if (x$data.type==TRUE) {
          cat("\n","Dependent variable a vector of counts.","\n")
                        } else {
          cat("\n","Dependent variable is a list of frequency distributions for counts","\n")
                               } # end of data.type==

# Checking for truncation i.e. ltvalue and utvalue are not NA
       if (is.na(x$ltvalue)==FALSE) { cat("\n","distribution truncated below at ",x$ltvalue) }
       if (is.na(x$utvalue)==FALSE) { cat("\n","distribution truncated above at ",x$utvalue,"\n") }

       if (is.null(x$converged)==FALSE) {
          cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 
             0.85)), "", sep = "\n")
          cat("Model type        :",x$model.type,"\n")
          cat("Model name        :",x$model.name,"\n")
          if (x$model.type=="mean and scale-factor") { 
                         cat("Link scale-factor : log","\n") }
          offsetid.mean     <- sum(x$offset.mean)
          offsetid.scalef <- sum(x$offset.scalef)
          if ((offsetid.mean!=0) | (offsetid.scalef!=0)) { 
              cat("non zero offsets in linear predictors","\n") }  

          npar <- length(x$optim$par)

          if ((x$model.name=="Faddy distribution") | 
              (x$model.name=="Faddy distribution fixed b")) { 
             if (x$model.name=="Faddy distribution") { wks <- npar - 2
                                              } else { wks <- npar - 1 }
          if (abs(x$optim$par[wks]-1)<1.e-6) {
             cat("Boundary for Faddy distribution c of 1 has been reached","\n") 
             cat("hence its se is set to NA.","\n") }} 

          cat(paste("\n","Coefficients (model for mean with", x$link,"link)\n", sep = " "))
          print(x$coeff.table.mean)    
  
          if (is.null(x$coeff.table.scalef)==FALSE) {
              cat(paste("\n","Coefficients (model for scale-factor with log link)\n")) 
              print(x$coeff.table.scalef) } # end if is.null

          if ((x$model.name=="negative binomial fixed b") | 
              (x$model.name=="Faddy distribution fixed b")) { 
             cat(paste("\n","Value of fixed b", x$fixed.b,"\n", sep = " "))
                                                            } # end of if model.name

          if (is.null(x$weights)==FALSE) {
             cat("\n","Maximum weighted likelihood regression.")
             if (x$data.type==TRUE) {
                cat("\n","Vector of weights used.","\n")
                             } else {
                cat("\n","List of weights used.","\n") }
             if (is.null(attr(x$weights, which="normalize"))==FALSE) {
                if (attr(x$weights, which="normalize")==TRUE) {
                   cat("Normalization to a value of",
                    attr(x$weights, which="norm.to.n"),".\n", sep = " ") }}
                                    } # end of is.null(weights)

          if (is.na(x$loglik)==TRUE) { cat("Log-likelihood is NA","\n")
                                          } else { 
             cat("\n","Type of estimator: ML (maximum likelihood)")
             cat("\n","Log-likelihood:",x$loglik,"on",length(x$optim$par),"Df", sep=" ")
             if (length(x$optim$par)==1) { 
                cat("\n","Single parameter Poisson so no use of optim", sep=" ")
                                  } else {
                optim.method <- x$method
                if (optim.method=="Nelder-Mead") { 
                   cat("\n","Number of iterations:",x$optim$counts[1],"of optim method",optim.method,sep=" ","\n") 
                                         } else { gradient.method <- attr(x$method,which="grad.method")
                   cat("\n","Number of iterations:",x$optim$counts[1],"of optim method",optim.method, 
                       "gradient method",gradient.method,sep=" ","\n")
                   cat("\n final gradients of parameters \n")
                   print(x$gradient) } } # end of if length(x$optim$par)=1
             code <- list(c("successful"),
                          c("iteration limit max has been reached"),
                          c(" "),c(" "),c(" "),c(" "),
                          c(" "),c(" "),c(" "),c(" "),
                          c("degeneracy of the Nelder-Mead"))
             wks <- attr(x$converged, which="code") + 1
             cat("\n","return code",attr(x$converged, which="code"),code[[as.numeric(wks)]],"\n", sep=" ")     
                 } # end of if is.na(x$loglik

                  } else { 
                       cat("\n","Failure of checks on entry arguments to CountsEPPM")
                       cat("\n","or numerical derivative calculations failed.")
                         } # end of if (is.null(converged)

         }

CountsEPPM.summary <-
function(output.fn) {
          cat('\n')
          cat('Model type:',output.fn$model.type,'\n')
          cat('Model     :',output.fn$model,'\n')
          cat('link      :','log','\n')   
          offsetid_mean     <- sum(output.fn$offset.mean)
          offsetid_variance <- sum(output.fn$offset.variance)
          if ((offsetid_mean!=0) | (offsetid_variance!=0)) { 
              cat('non zero offsets in linear predictors','\n') }   
          if (output.fn$scale.factor.model=='yes') { 
              cat('scale factor model fitted','\n') }   
          if (is.na(output.fn$estses[1,3])==TRUE) { 
              cat('Determinant of hessian matrix is zero','\n') 
              cat('or the hessian matrix is ill conditioned.','\n') }
          cat('Parameter estimates and se\'s','\n')
          print.data.frame(output.fn$estses,row.names=FALSE)
          if (is.na(output.fn$loglikelihood)==TRUE) { cat('loglikelihood is NA','\n')
          } else { cat('\n')
                   cat('log likelihood ',output.fn$loglikelihood,'\n')
                   cat('\n')
                   numpar <- length(output.fn$estses[,1]) 
                   AIC <- -2*output.fn$loglikelihood + 2*numpar 
                   cat('AIC',AIC,'\n')
                 } # end of is.na(output.fn$loglikelihood if
                                   }

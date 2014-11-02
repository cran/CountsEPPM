CountsEPPM.summary <-
function(output.fn) { cat('\n')
          cat('Model type:',output.fn$model.type,'\n')
          cat('Model     :',output.fn$model,'\n')
          if ((output.fn$model=='Faddy distribution fixed b') | 
              (output.fn$model=='negative binomial fixed b') | 
              (output.fn$model=='general fixed b')) { 
             cat('b fixed at:',output.fn$fixed.b,'\n') } # end if output$model
          cat('Link for mean         :','log','\n')   
# Introduction of link information for scale factor and variance in Version 1.01
          if (output.fn$model.type=='mean and variance') { 
             if (output.fn$scale.factor.model=='yes') { 
                          cat('Link for scale factor :','log','\n') 
                 } else { cat('Link for variance     :','log','\n') }}   
          offsetid_mean     <- sum(output.fn$offset.mean)
          offsetid_variance <- sum(output.fn$offset.variance)
          if (offsetid_mean!=0) { 
             cat('non zero offsets in linear predictor for mean','\n') }   
          if (offsetid_variance!=0) { 
              cat('non zero offsets in linear predictor for variance','\n') }   
          if (output.fn$model.type=='mean and variance') { 
             if (output.fn$scale.factor.model=='no') { 
                 cat('variance model fitted','\n') 
                    } else {
                 cat('scale factor model model fitted','\n') }}   
          if (is.na(output.fn$estses[1,3])==TRUE) { 
              cat('Determinant of hessian matrix is zero','\n') 
              cat('or the hessian matrix is ill conditioned.','\n') }

          if (is.na(output.fn$loglikelihood)==FALSE) {
             if ((output.fn$model=="Faddy distribution") | 
                 (output.fn$model=="Faddy distribution fixed b")) { 
                  npar   <- nrow(output.fn$estses)
                  nparm1 <- npar - 1
               if (output.fn$model=="Faddy distribution") {
                  wk.c <- round(output.fn$estses[nparm1,2],digits=7)
                  } else {
                  wk.c <- round(output.fn$estses[npar,2],digits=7) }
               if (wk.c==1) {
                 cat('Boundary for c of 1 has been reached','\n') 
                 cat('hence its se is set to NA.','\n') }}} # end is.na 
               
          cat('Parameter estimates and se\'s','\n')
# Introduced in Version 1.01 to deal with the missing name se for 
# the single parameter situation.
          names(output.fn$estses) <- c('name','Estimates','se')
          print.data.frame(output.fn$estses,row.names=FALSE)
          if (is.finite(output.fn$loglikelihood)==FALSE) {
                  cat('loglikelihood is NA or -Inf suggesting an inappropriate model','\n')
          } else { cat('\n')
                   cat('log likelihood ',output.fn$loglikelihood,'\n')
                   AIC  <- -2*output.fn$loglikelihood + 
                            2*length(output.fn$estses[,1])  
                   cat('AIC',AIC,'\n')
                 } # end of is.na(output.fn$loglikelihood if
                                   }

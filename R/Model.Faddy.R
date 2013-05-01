Model.Faddy <-
function(parameter,model.type,model,covariates.matrix.mean,
                offset.mean,vnmax) {
#  data as number of trials & number of successes
#  the parameters log(b) & c are scalars and are the last parameters in the parameter vector   
   nobs <- nrow(covariates.matrix.mean) 
#  covariate matrix for parameter a
   npar.one <- ncol(covariates.matrix.mean)
   npar <- npar.one
   threeparameter <- rep(0,3)  
   va   <- matrix(c(rep(1,nobs)),ncol=1)
   r.parameter <- rep(0,npar.one) 
   r.parameter <- parameter[1:npar.one]
   va <- exp(covariates.matrix.mean%*%r.parameter + offset.mean)
   if (model=="Poisson") { b <- 1
                           c <- 0 } # end of model
   if (model=="negative binomial") { npar <- npar + 1 
                                     b <- exp(parameter[npar])
                                     c <- 1 } # end of model
   if (model=="Faddy distribution") { nparm1 <- npar + 1
                                      npar   <- npar + 2 
                                      b <- exp(parameter[nparm1])
                                      c <- parameter[npar] } # end of model
   probabilities <- rep(list(0),nobs)
   threeparameter[2] <- b
   threeparameter[3] <- c 
   for ( i in 1:nobs) { threeparameter[1] <- va[i] 
                        nmax              <- vnmax[i]
                        nmax1             <- nmax + 1 
                        vid               <- c(0:nmax)
                        if (c<=1) { probability <- Faddyprob.general(threeparameter,nmax)
# c > 1 indicated, all probabilities set to very small value
                                  } else { probability <- rep(1.e-8,nmax1) } # end of if
                        probabilities[[i]] <- probability }
   output <- list(model=model,estimates=parameter,probabilities=probabilities)
   return(output)        }

Model.FaddyJMV.general <-
function(parameter,covariates.matrix.mean,
                    covariates.matrix.variance,offset.mean,offset.variance, 
                    scale.factor.model,vnmax) {
   loglikelihood <- 0
   nobs          <- nrow(covariates.matrix.mean)
   npar.mean     <- ncol(covariates.matrix.mean)
   npar.variance <- ncol(covariates.matrix.variance)
   npar <- npar.mean + npar.variance
   numpar <- length(parameter) 
   r.parameter.mean <- rep(0,npar.mean) 
   r.parameter.mean <- parameter[1:npar.mean] 
   probabilities <- rep(list(0),nobs)
   for ( i in 1:nobs) {
      nmax  <- vnmax[i]
      nmax1 <- nmax + 1
      vnum  <- c(0:nmax)
# log link function for mean
      lp.mean <- covariates.matrix.mean[i,]%*%r.parameter.mean + offset.mean[i]
      cmean <- exp(lp.mean)
# log-linear function for variance 
      r.parameter.variance <- rep(0,npar.variance) 
      wks <- npar.mean + 1
      r.parameter.variance <- parameter[wks:npar] 
      if (scale.factor.model=='no') { lp.variance <- 
             covariates.matrix.variance[i,]%*%r.parameter.variance + 
                                                 offset.variance[i] }
      if (scale.factor.model=='yes') { lp.variance <- 
             covariates.matrix.variance[i,]%*%r.parameter.variance + log(cmean) + 
                                                 offset.variance[i] }
      cvariance   <- exp(lp.variance)
# start of over-dispersed
# extra parameter log(b) 
      wks <- npar + 1
      b <- exp(parameter[wks])
# iteration on c
     aa <- log(cmean/b+1)
     bb=cvariance/(cmean+b)
     if ((exp(aa)-1)>bb) { x <- aa
       xinc <- ((exp(aa)-1-bb)/aa)/((exp(aa)*(aa-1)+1)/aa^2)
       iter <- 0
       while ((iter<51) & (abs(xinc)>=1.e-10)) { x <- x-xinc
          xinc <- ((exp(x)-1)/x-bb/aa)/((exp(x)*(x-1)+1)/x^2)
          iter <- iter + 1
                                               } # end of while 
# c <= 1 solution indicated
      c1 <- (x/aa+1)/2
# c not close to 1
      if (c1<(1-10^(-8))) { a1 <- ((cmean+b)^(1-c1)-b^(1-c1))/(1-c1)
# c close to 1
                          } else { a1 <- log(1+cmean/b) }
      va    <- rep(a1,nmax1)
      vb    <- rep(b,nmax1)
      vc    <- rep(c1,nmax1)
      if ((a1<Inf) & (is.na(a1)==FALSE)) {
         if (c1==0) { vlambda <- va
                    } else { vlambda <- va*(vb+vnum)^vc
                             for (j in 1:nmax1) { 
                             if ((vlambda[j]>10^10) | (vlambda[j]==Inf)) { 
                                                 vlambda[j] <- 10^10 } } } 
         probability <- EPPMprob(vlambda) 
                    } else { probability <- rep(1.e-8,nmax1) }
# end of 1st part of if ((exp(aa)-1)>bb) 
                         } else { probability <- rep(1.e-8,nmax1)
# c > 1 indicated, all probabilities set to very small value
                              } # end of 2nd part of if ((exp(aa)-1)>bb) 
      probabilities[[i]] <- probability } # end of for loop
   model <- "general"
   output <- list(model=model,estimates=parameter,probabilities=probabilities)
   return(output)              }

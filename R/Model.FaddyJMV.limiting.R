Model.FaddyJMV.limiting <-
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
             covariates.matrix.variance[i,]%*%r.parameter.variance + 
                                   log(cmean) + offset.variance[i] }
      cvariance   <- exp(lp.variance)
# variance is > mean, all probabilities set to very small value
      if (cvariance>cmean) { probabilities[[i]] <- rep(1.e-8,nmax1) 
                           } else { 
# start of under-dispersed equations
# solve equation (10) for beta by Newton-Raphson iteration
         x  <- 0
         xinc <- (1 - cvariance/cmean) / (1/2)
         iter <- 0
         while ((iter<51) & (abs(xinc)>=1.e-10)) { x    <- x - xinc 
               xinc <- ((exp(x)-1)/x - cvariance/cmean) / ((x*exp(x)-exp(x)+1)/x^2)
               iter <- iter + 1
                                             } # end of while
         twoparameter <- rep(0,2)
# determine beta
         twoparameter[2] <- x/2/cmean
# determine alpha from equation (9)
         twoparameter[1] <- (1-exp(-twoparameter[2]*cmean)) / twoparameter[2]
         probability <- Faddyprob.limiting(twoparameter,nmax)
         probabilities[[i]] <- probability } # end of if variance > mean
                 } # end of for loop
   model <- "limiting"
   output <- list(model=model,estimates=parameter,probabilities=probabilities)
   return(output)              }
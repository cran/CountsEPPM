Model.FaddyJMV.limiting <-
function(parameter,link,covariates.matrix.mean,
                     covariates.matrix.scalef,offset.mean,offset.scalef, 
                     vnmax) {
   nobs          <- nrow(covariates.matrix.mean)
   npar.mean     <- ncol(covariates.matrix.mean)
   npar.scalef <- ncol(covariates.matrix.scalef)
   npar <- npar.mean + npar.scalef
   numpar <- length(parameter) 
   r.parameter.mean <- rep(0,npar.mean) 
   r.parameter.mean <- parameter[1:npar.mean] 
   r.parameter.scalef <- rep(0,npar.scalef) 
   wks <- npar.mean + 1
   r.parameter.scalef <- parameter[wks:npar] 
   probabilities <- rep(list(0),nobs)

# link function for mean
   lp.mean <- covariates.matrix.mean%*%r.parameter.mean + offset.mean
# log-linear function for scale-factor
   lp.scalef <- covariates.matrix.scalef%*%r.parameter.scalef + offset.scalef 
   lp.mean   <- attr(link, which="mean")$linkinv(lp.mean)
   lp.scalef <- exp(lp.scalef)
# vectors of Faddy distribution parameters
   out.valpha <- rep(0,nobs)
   out.vbeta  <- rep(0,nobs)

   for ( i in 1:nobs) {
      nmax  <- vnmax[i]
      nmax1 <- nmax + 1
      vnum  <- c(0:nmax)
      cmean <- lp.mean[i]
      cvariance <- lp.scalef[i]*cmean
# variance is > mean, all probabilities set to very small value
      if (cvariance>cmean) { probabilities[[i]] <- rep(0,nmax1) 
                           } else { 
# start of under-dispersed equations
# solve equation (10) for beta by Newton-Raphson iteration
         x    <- 0
         xinc <- 2*(1 - cvariance/cmean)
         iter <- 0
         while ((iter<101) & (abs(xinc)>=1.e-10) & (is.finite(xinc)==TRUE)) { 
               x <- x - xinc 
               xinc <- ((exp(x)-1)/x - cvariance/cmean) / ((x*exp(x)-exp(x)+1)/x^2)
               iter <- iter + 1 } # end of while
         twoparameter <- rep(0,2)
# determine beta
         twoparameter[2] <- x/2/cmean
# determine alpha from equation (9)
         if ((cvariance==cmean) | (twoparameter[2]==0)) { twoparameter[1] <- cmean
                  } else {
            twoparameter[1] <- (1-exp(-twoparameter[2]*cmean)) / twoparameter[2] }
         probability <- Faddyprob.limiting(twoparameter,nmax)
         out.valpha[i] <- twoparameter[1] 
         out.vbeta[i]  <- twoparameter[2] 
         probabilities[[i]] <- probability } # end of if variance > mean
                 } # end of for loop 
   output <- list(model.name="limiting",estimates=parameter,probabilities=probabilities,
                  FDparameters=data.frame(out.valpha,out.vbeta))
   return(output)              }

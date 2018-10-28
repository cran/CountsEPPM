LL.Regression.Counts <-
function(parameter,model.type,model.name,link,list.data,
               covariates.matrix.mean,covariates.matrix.scalef,
               offset.mean,offset.scalef,ltvalue,utvalue,
               fixed.b,weights,grad.method) {
   nobs   <- nrow(covariates.matrix.mean) 
   vnmax  <- sapply(list.data,length) - 1
   probabilities <- Model.Counts(parameter,model.type,model.name,link=link,covariates.matrix.mean,
                                 covariates.matrix.scalef,offset.mean,offset.scalef,
                                 fixed.b,vnmax)$probabilities 
# 1.e-323 = exp(-743.7469) is the smallest probability value not to give -Inf
# 1.e-325 has been set as value for -infinity for log likelihoods
      small.loglikelihood    <- log(1.e-323)
      infinity.loglikelihood <- log(1.e-325)

# Calculation of log likelihood
      vlogl <- rep(0,nobs)
      vlogl <- sapply(1:nobs, function(i) {
         probability <- probabilities[[i]]
         count <- list.data[[i]] 
         nmax1 <- vnmax[i] + 1

# With overly large values of the parameters and small counts it is possible to have
# probabilities exactly 0 returned. Also LRTruncation sets the probabilities below
# ltvalue and above utvalue to NA.

# rounding error can cause the sum of the probabilities to <0 or >1
# so rounded to 10 decimal places
        wks <- round(sum(probability),digits=10)
        if ((is.finite(wks)==FALSE) | (wks<=0) | (wks>1)) {
            vlogl[i] <- sum(count)*small.loglikelihood  

                                                    } else {

         log.prob <- rep(0,length(probability))
         log.prob <- sapply(1:length(probability), function(j) {
                       if ((is.na(probability[j])==TRUE) | (probability[j]<0)) { 
                                   log.prob[j] <- -1.e+20
                          } else { 
                           if (probability[j]==0) { log.prob[j] <- 0
                                   } else { log.prob[j] <- log(probability[j]) }} 
                         } ) # end of sapply

           if ((is.na(ltvalue)==FALSE) | (is.na(utvalue)==FALSE)) {
              if (is.na(ltvalue)==FALSE) { wks1 <- ltvalue + 2 
                                         } else { wks1 <- 1 }
              if (is.na(utvalue)==FALSE) { wks2 <- utvalue  
                                         } else { wks2 <- nmax1 }
              rev.probability <- LRTruncation(probability,ltvalue,utvalue)
              wks  <- round(sum(rev.probability[wks1:wks2]),digits=10)
              if (is.finite(wks)==FALSE) { wks <- 0 }
              wkv2 <- count[wks1:wks2]
              if ((is.na(utvalue)==FALSE) & (wks!=1)) { 
                 vlogl[i] <- infinity.loglikelihood 
                                                      } else {
                 wkv1 <- log(rev.probability[wks1:wks2])
                 nlen <- length(wkv1)
                 wkv1 <- sapply(1:nlen, function(j) 
                    if (is.finite(wkv1[j])==FALSE) { wkv1[j] <- small.loglikelihood 
                                            } else { wkv1[j] <- wkv1[j] } ) # end of sapply
                       if (is.null(weights)==TRUE) { vlogl[i] <- t(wkv1)%*%wkv2
                                            } else {
                          if (is.list(weights)==TRUE) { wk.wts <- weights[[i]]
                                            } else {
                             wk.wts <- c(rep(weights[[i]],length(wkv2)))
                       vlogl[i] <- t(wk.wts*wkv1)%*%wkv2 } } # end of is.list(weights)
                                                       } # end wks!=1
                                                                    } else {
                 wkv1 <- log(probability)
                 wkv1 <- sapply(1:nmax1, function(j) 
                    if (is.finite(wkv1[j])==FALSE) { wkv1[j] <- small.loglikelihood 
                                                   } else { wkv1[j] <- wkv1[j] } ) # end of sapply

                 if (is.null(weights)==TRUE) { vlogl[i] <- t(wkv1)%*%count
                                      } else {
                    if (is.list(weights)==TRUE) { wk.wts <- weights[[i]]
                                      } else {
                       wk.wts <- c(rep(weights[[i]],length(count))) } 
                    vlogl[i] <- t(wk.wts*wkv1)%*%count } # end of is.list(weights)
                                       } } } ) # end of first sapply

      if (sum((vlogl==0))==0) { loglikelihood <- sum(vlogl)
                       } else { loglikelihood <- -1.e+20 } 
   return(loglikelihood)                          }

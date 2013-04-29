F.Regression.Counts <-
function(parameter,model.type,model,list.counts,
               covariates.matrix.mean,covariates.matrix.variance,
               offset.mean,offset.variance,ltvalue,utvalue,
               optimization.method,scale.factor.model) {

   nobs <- nrow(covariates.matrix.mean) 
   vnmax <- rep(0,nobs)
   for ( i in 1:nobs) { vnmax[i] <- length(list.counts[[i]]) - 1 }
   output <- Model.Counts(parameter,model.type,model,covariates.matrix.mean,
                   covariates.matrix.variance,offset.mean,offset.variance,
                   scale.factor.model,vnmax) 
# Checking for improper probability distributions
    improper <- rep(0,nobs)
    for ( i in 1:nobs) { probability <- output$probabilities[[i]]
        nmax1 <- vnmax[i] + 1
        wks   <- sum((probability==1.e-8))
        if (wks==nmax1) { improper[i] <- 1 }
        wks <- sum(probability)
        wks <- round(wks,digits=8)
        if ((wks<0) | (wks>1))  { probability <- rep(1.e-8,nmax1)
                                  improper[i] <- 1 }
        wks <- 0
        wks <- wks + sum((probability<0) | (probability>1))
        if (wks>0) { probability <- rep(1.e-8,nmax1)
                     improper[i] <- 1 }
                       } # end of for loop

      loglikelihood <- 0 
      for ( i in 1:nobs) { probability <- output$probabilities[[i]]
         count <- list.counts[[i]] 
         nmax1 <- vnmax[i] + 1
         if ((is.na(ltvalue)==FALSE) | (is.na(utvalue)==FALSE)) { 
             rev.probability <- LRTruncation(probability,ltvalue,utvalue)
             if (is.na(ltvalue)==FALSE) { wks1 <- ltvalue + 2 
                                        } else { wks1 <- 1 }
             if (is.na(utvalue)==FALSE) { wks2 <- utvalue  
                                        } else { wks2 <- nmax1 }

             wkv1 <- rev.probability[wks1:wks2]
             wkv2 <- count[wks1:wks2]
             if (optimization.method=='optim') { 
                loglikelihood  <- loglikelihood + t(log(wkv1))%*%wkv2 }
             if (optimization.method=='nlm')   { 
                loglikelihood  <- loglikelihood - t(log(wkv1))%*%wkv2 } 
                                                                 } else {
# With overly large values of the parameters and small counts it is possible to have
# probabilities exactly 0 returned. Also LRTruncation sets the probabilities below
# ltvalue and above utvalue to 0.
            probability <- probability*(probability>0) + 1.e-8*(probability==0)
            if (optimization.method=='optim') { loglikelihood  <- loglikelihood + 
                                                    t(log(probability))%*%count }
            if (optimization.method=='nlm')   { loglikelihood  <- loglikelihood - 
                                                    t(log(probability))%*%count }
                         } } # end of for loop
   return(loglikelihood)                          }

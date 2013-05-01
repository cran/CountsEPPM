Faddyprob.limiting <-
function(parameter,nmax) { 
   nmax1 <- nmax + 1
   vlambda <- rep(0,nmax1) 
# if parameter[1] is Inf (infinity) na's are produced in all elements of vlambda 
# except the first
   if (parameter[1]==Inf) { probability <- rep(1.e-8,nmax1) 
       } else { if (parameter[2]==0) { vlambda <- rep(parameter[1],nmax1)  
                              } else { vbeta <- rep(parameter[2],nmax1)
                                       vnum  <- c(0:nmax)
                                       vlambda <- parameter[1]*exp(parameter[2]*vnum) 
                                       for (j in 1:nmax1) {                            
                                          if ((vlambda[j]>10^10) | (vlambda[j]==Inf)) {
                                             vlambda[j] <- 10^10 } } } 
         probability <- EPPMprob(vlambda) }
    return(probability)   }

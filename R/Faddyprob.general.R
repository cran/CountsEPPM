Faddyprob.general <-
function(parameter,nmax) {
   nmax1 <- nmax + 1
   vlambda <- rep(0,nmax1) 
   va   <- rep(parameter[[1]],nmax1)
   vb   <- rep(parameter[[2]],nmax1)
    c   <- parameter[[3]]
   vc   <- rep(c,nmax1)
   vnum <- c(0:nmax)
#  case of c > 1 dealt with in Model.Faddy
   if (c==0) { vlambda <- va } 
   if ((c<(1-10^(-8))) & (c!=0)) { vlambda <- va*(vb+vnum)^vc
      for (j in 1:nmax1) { 
         if ((vlambda[j]>10^10) | (vlambda[j]==Inf)) {
            vlambda[j] <- 10^10 } } }
# c very close to 1 i.e. negative binomial
               if (c>=(1-10^(-8))) { vlambda <- va*(vb+vnum) }
               probability <- EPPMprob(vlambda)
   return(probability)   }

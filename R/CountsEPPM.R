CountsEPPM <-
function(formula,data,model.type='mean and variance',model='general',
                  initial=c(NA),ltvalue=NA,utvalue=NA,optimization.method='optim',
                  scale.factor.model='no') {

   nlist       <- length(data)
   nlistm1     <- nlist - 1
   list.counts <- data$list.counts
   nobs        <- length(list.counts)
   mean.obs       <- rep(0,nobs)
   variance.obs   <- rep(0,nobs)
   for ( i in 1:nobs ) { 
      count        <- list.counts[[i]]
      nmax1        <- length(count)
      nmax         <- nmax1 - 1
      cnum         <- 0:nmax
      ncount       <- sum(count)  
      mean.obs[i]     <- sum(cnum*count)/ncount 
      variance.obs[i] <- (sum(cnum*cnum*count) - ncount*mean.obs[i]*mean.obs[i]) / (ncount - 1)
                 } # end of for loop
    if (scale.factor.model=='yes') { scalef.obs <- variance.obs/mean.obs 
        wkdata <- data.frame(mean.obs,variance.obs,scalef.obs,data$covariates) 
                                   } else {
        wkdata <- data.frame(mean.obs,variance.obs,data$covariates) }

# Checking arguments of function
   inerror <- 1

# Checking that formula has either 1 lhs and 1 rhs, or 2 ihs and 2 rhs,
# corresponding to just a formula for mean, and a formula for both mean
# and variance 
   FBoth  <- Formula(formula) 
   mfBoth <- model.frame(FBoth,data=wkdata)
   lenFB <- length(FBoth)
   if ((lenFB[1]==1) & (lenFB[2]==1)) {
      inerror <- 0
      covariates.matrix.mean     <- model.matrix(FBoth,data=mfBoth,lhs=1,rhs=1)
      covariates.matrix.variance <- matrix(c(rep(1,nrow(covariates.matrix.mean))),ncol=1)
      resp.mean <- model.part(FBoth,data=mfBoth,lhs=1)
      if (attributes(resp.mean)$names!='mean.obs') {
         cat('\n','WARNING: mean response variate is not named mean.obs','\n') 
         cat('as the name is not used the function continues','\n') }   
                                      } # mean formula only
   if ((lenFB[1]==2) & (lenFB[2]==2)) {
      inerror <- 0
      covariates.matrix.mean     <- model.matrix(FBoth,data=mfBoth,lhs=1,rhs=1)
      covariates.matrix.variance <- model.matrix(FBoth,data=mfBoth,lhs=2,rhs=2)
      resp.mean     <- model.part(FBoth,data=mfBoth,lhs=1)
      resp.variance <- model.part(FBoth,data=mfBoth,lhs=2)
      if (attributes(resp.mean)$names!='mean.obs') {
         cat('\n','WARNING: mean response variate is not named mean.obs','\n') 
         cat('as the name is not used the function continues','\n') }
         if (scale.factor.model=='yes') { 
            if (attributes(resp.variance)$names!='scalef.obs') {
         cat('\n','WARNING: variance response variate is not named scalef.obs','\n') 
         cat('as the name is not used the function continues','\n') }
                                        } else {      
            if (attributes(resp.variance)$names!='variance.obs') {
         cat('\n','WARNING: variance response variate is not named variance.obs','\n') 
         cat('as the name is not used the function continues','\n') } }
                                      } # mean and variance formula 
   if (inerror==1) {
      cat('\n','formula error, lhs & rhs not both the same value of 1 or 2','\n') }

# Checking for correct combinations of model.type and model
   if ((model.type!='mean only') & (model.type!='mean and variance')) {
      inerror <- 1
      cat('\n','unknown model.type','\n') }
   if (model.type=='mean only') {
      if ((model!='Poisson') & (model!='negative binomial') & (model!='Faddy distribution')) {
         inerror <- 1
         cat('\n','unknown model for this model.type','\n') } }
   if (model.type=='mean and variance') {
      if ((model!='general') & (model!='limiting')) {
         inerror <- 1
         cat('\n','unknown model for this model.type','\n') }
# Checking that the formula has two parts 
      if ((lenFB[1]==1) & (lenFB[2]==1)) {
         cat('\n','variance model set to default of an intercept only','\n') } }

# Checking that list.counts is a list  
if (inerror==0) {   
   if (is.list(list.counts)==FALSE) { inerror <- 1
      cat('\n','list.counts is not a list','\n') }

# Checking for offsets i.e. offset.mean, offset.variance are not NA
    if (is.null(data$offset.mean)==TRUE) {
       offset.mean=c(rep(0,nobs)) } else {
          nobswk <- length(data$offset.mean)
          if (nobswk!=nobs) { offset.mean=c(rep(0,nobs)) 
             cat('\n','WARNING: length of offset.mean not same as number of grouped count vectors','\n')
             cat('\n','calculations proceed but with the default for offset.mean used','\n')
                            } else {
             offset.mean=data$offset.mean } } # end if offset.mean
    if (is.null(data$offset.variance)==TRUE) { 
       offset.variance=c(rep(0,nobs)) } else {
          nobswk <- length(data$offset.variance)
          if (nobswk!=nobs) { offset.variance=c(rep(0,nobs)) 
             cat('\n','WARNING: length of offset.variance not same as number of grouped count vectors','\n')
             cat('\n','calculations proceed but with the default for offset.variance used','\n')
                            } else {
             offset.variance=data$offset.variance } } # end if offset.variance
    wkdata <- data.frame(wkdata,offset.mean,offset.variance)

# Setting up initial estimates of parameters if not input 
   if (is.na(initial[1])==TRUE) { 
# Using glm to obtain initial estimates 
      glm.Results <- glm(formula(FBoth,lhs=1,rhs=1),family=gaussian(link="log"),
                         data=wkdata,offset=offset.mean) 
      initial.mean <- coefficients(glm.Results)
      if (model.type=='mean and variance') {
         if ((lenFB[1]==2) & (lenFB[2]==2)) {
            glm.Results <- glm(formula(FBoth,lhs=2,rhs=2),family=gaussian(link="log"),
                               data=wkdata,offset=offset.variance) 
                                            } else {
            glm.Results <- glm(formula(variance.obs~1),family=gaussian(link="log"),data=wkdata) 
                                                   } 
            initial.variance <- coefficients(glm.Results)
                                           } # end of if mean and variance
      if (model.type=='mean only') {
          if (model=='Poisson') { parameter <- initial.mean
                           names(parameter) <- names(initial.mean) }
          if (model=="negative binomial")  { parameter <- c(initial.mean,0) 
                           names(parameter) <- c(names(initial.mean),"log(b)") }
          if (model=="Faddy distribution") { parameter <- c(initial.mean,0,0) 
                           names(parameter) <- c(names(initial.mean),"log(b)","c") }
                                   } # of if model.type mean only
      if (model.type=='mean and variance') {
         if (model=='general')  { parameter <- c(initial.mean,initial.variance,0) 
            names(parameter) <- c(names(initial.mean),names(initial.variance),"log(beta)") }
         if (model=='limiting') { parameter <- c(initial.mean,initial.variance) 
             names(parameter) <- c(names(initial.mean),names(initial.variance)) }
                                   } # of if model.type mean and variance
                                   } else {
# Checking length of input initial against model value
      npar.one <- ncol(covariates.matrix.mean)
      numpar <- length(initial) 
      if (model.type=="mean only") { npar <- npar.one 
              if (model=="negative binomial")  { npar <- npar + 1 }
              if (model=="Faddy distribution") { npar <- npar + 2 } }
      if (model.type=="mean and variance") { 
              npar.two <- ncol(covariates.matrix.variance)
              npar <- npar.one + npar.two 
              if (model=="general")  { npar <- npar + 1 } }
      if (numpar!=npar) { inerror <- 1  
            cat('\n','number of parameters error','\n') }
      parameter <- initial
      if (is.null(names(initial))==TRUE) { 
         cat('\n','WARNING: initial has no associated names','\n')
         names(parameter) <- 1:numpar } # end if is.null(names)
                           } # end if is.null(initial[1])
                                  } # end of inerror=0

# Checking that scale.factor.model is only yes if model is under- 
# or over-dispersed means & variances  
if (inerror==0) {   
    if (scale.factor.model=='yes') { if (model.type!='mean and variance') { 
       inerror <- 1 
       cat('\n','scale factor model should be set to no','\n') } }
                } # end of inerror=0

# Checking for truncation i.e. ltvalue and utvalue are not NA
    if (is.na(ltvalue)==FALSE) { cat('\n','distribution truncated below at ',ltvalue) }
    if (is.na(utvalue)==FALSE) { cat('\n','distribution truncated above at ',utvalue) }

# start of if inerror if else
if (inerror==0) {   

# Seting up vector vnmax and calculating likelihood for initial estimates
    vnmax <- rep(0,nobs)                 
    for ( i in 1:nobs) { vnmax[i] <- length(list.counts[[i]]) - 1 } # end of for loop

    output <- Model.Counts(parameter,model.type,model,covariates.matrix.mean,
              covariates.matrix.variance,offset.mean,offset.variance,scale.factor.model,
              vnmax) 

# Checking for improper probability distributions for initial estimates
    improper <- rep(0,nobs)
    for ( i in 1:nobs) { probability <- output$probabilities[[i]]
        nmax1 <- vnmax[i] + 1
        wks   <- sum((probability==1.e-8))
        if (wks==nmax1) { improper[i] <- 1 }
        wks <- sum(probability)
# using round to handle totals of probabilities very slightly greater than 1 
# due to rounding errors
        wks <- round(wks,digits=8)
        if ((wks<0) | (wks>1))  { probability <- rep(1.e-8,nmax1)
                                  improper[i] <- 1 }
        wks <- sum((probability<0) | (probability>1))
        if (wks>0) { probability <- rep(1.e-8,nmax1)
                     improper[i] <- 1 }
                       } # end of for loop

# Fitting model using optim or nlm

          if (optimization.method=='optim') {
             Results  <- optim(parameter,fn=F.Regression.Counts,gr=NULL,model.type,model,list.counts,
                               covariates.matrix.mean,covariates.matrix.variance,
                               offset.mean,offset.variance,ltvalue,utvalue,optimization.method,
                               scale.factor.model,method="Nelder-Mead",hessian=TRUE,
                               control=list(fnscale=-1,trace=0,maxit=1000)) 
             convergence <- c("successful","iteration limit max has been reached",
                              rep(" ",8),"degeneracy of the Nelder-Mead simplex")
             cat('\n','optimization method optim:','\n')
             wks <- Results$convergence + 1
             cat('convergence ',Results$convergence,convergence[wks],'\n')
             estimates <- Results$par } # end optimization.method=optim
          if (optimization.method=='nlm') {
             Results <- nlm(F.Regression.Counts,parameter,model.type,model,list.counts,
                            covariates.matrix.mean,covariates.matrix.variance,
                            offset.mean,offset.variance,ltvalue,utvalue,optimization.method,
                            scale.factor.model,hessian=TRUE,fscale=1,print.level=0,stepmax=1,
                            gradtol=1e-8,steptol=1e-10,iterlim=500)
             code <- c("relative gradient is close to zero, current iterate is probably solution",
                       "succesive iterates within tolerance, current iterate is probably solution",
                       "last global step failed to locate a point lower than estimate. 
                        Either estimate is an approximate local minimum of the function 
                        or steptol is too small","iteration limit max has been reached",
                       "Maximum step size stepmax exceeded five consecutive times. 
                        Either the function is unbounded below, becomes asymptotic to a 
                        finite value from above in some direction, or stepmax is too small")
             cat('\n', 'optimization method nlm:','\n')
             cat('code       ',Results$code,code[Results$code],'\n')
             cat('iterations ',Results$iterations,'\n')            
             estimates <- Results$estimate 
             names(estimates) <- names(parameter) } # end optimization.method=nlm 
          output <- Model.Counts(estimates,model.type,model,covariates.matrix.mean,
                   covariates.matrix.variance,offset.mean,offset.variance,scale.factor.model,
                   vnmax) 
          npar  <- length(estimates)
          deter <- det(Results$hessian)
          if (deter==0) { se <- rep(NA,npar) 
                        } else { listlen <- nrow(covariates.matrix.mean)
                                 if ((npar==1) & (listlen==1)) { 
                if (optimization.method=='optim') { se <- sqrt(-1/Results$hessian) }
                if (optimization.method=='nlm')   { se <- sqrt( 1/Results$hessian) }
                                  } else { condition <- rcond(Results$hessian)
# function rcond is from the package Matrix and gives the (reciprocal) condition number
# near 0 is ill-conditioned, near 1 is well conditioned
                                           if (condition>1e-8) { 
                if (optimization.method=='optim') { varcov <- - solve(Results$hessian) }
                if (optimization.method=='nlm')   { varcov <-   solve(Results$hessian) }
                                                                 se     <- sqrt(diag(varcov)) 
                                                               } else { se <- rep(NA,npar) } } } 
          nobs <- nrow(covariates.matrix.mean) 
# Calculation of log likelihood
      loglikelihood <- 0 
      for ( i in 1:nobs) { probability <-  output$probabilities[[i]]     
         count <- list.counts[[i]] 
         nmax1 <- vnmax[i] + 1
         if ((is.na(ltvalue)==FALSE) | (is.na(utvalue)==FALSE)) { 
             rev.probability <- LRTruncation(probability,ltvalue,utvalue)
             if (is.na(ltvalue)==FALSE) { wks1 <- ltvalue + 2 
                                        } else { wks1 <-1 }
             if (is.na(utvalue)==FALSE) { wks2 <- utvalue 
                                        } else { wks2 <- nmax1 }
             wkv1 <- rev.probability[wks1:wks2]
             wkv2 <- count[wks1:wks2]
             loglikelihood  <- loglikelihood + t(log(wkv1))%*%wkv2 
                                                                 } else {
# With overly large values of the parameters and small counts it is possible to have
# probabilities exactly 0 returned. Also LRTruncation sets the probabilities below
# ltvalue and above utvalue to 0.
             probability <- probability*(probability>0) + 1.e-8*(probability==0)
             loglikelihood  <- loglikelihood + t(log(probability))%*%count }
                      } # end of for loop
     name <- names(estimates) 
     output.fn <- list(model.type=model.type,model=model,
                       covariates.matrix.mean=covariates.matrix.mean,
                       covariates.matrix.variance=covariates.matrix.variance,
                       offset.mean=offset.mean,offset.variance=offset.variance,
                       ltvalue=ltvalue,utvalue=utvalue,
                       scale.factor.model=scale.factor.model,
                       estses=data.frame(name,estimates,se),
                       vnmax=vnmax,loglikelihood=loglikelihood,
                       mean.obs=mean.obs,variance.obs=variance.obs)
                 } else {output.fn <- list(model.type=model.type,model=model,
                           covariates.matrix.mean=covariates.matrix.mean,
                           covariates.matrix.variance=covariates.matrix.variance,
                           offset.mean=offset.mean,offset.variance=offset.variance,
                           ltvalue=ltvalue,utvalue=utvalue,
                           scale.factor.model=scale.factor.model,
                           estses=data.frame(NA,NA,NA),
                           vnmax=NA,loglikelihood=NA,
                           mean.obs=mean.obs,variance.obs=variance.obs) } # inerror if else
          return(output.fn)     }

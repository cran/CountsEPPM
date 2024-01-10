CountsEPPM <-
function(formula,data,subset=NULL,na.action=NULL,weights=NULL,
                  model.type="mean and scale-factor",model.name="general",
                  link="log",initial=NULL,ltvalue=NA,utvalue=NA,
                  method="Nelder-Mead",control=NULL,
                  fixed.b=NA) {

# Checking for correct combinations of model.type and model.name
   if (model.type!="mean only") {
      if ((model.name!="general") & (model.name!="general fixed b") & (model.name!="limiting")) {
          cat("\n","unknown model.name for this model.type","\n")
          return(object=NULL) }
                     } else {
      if ((model.name!="Poisson") & (model.name!="negative binomial") &
          (model.name!="negative binomial fixed b") & (model.name!="Faddy distribution") &
          (model.name!="Faddy distribution fixed b")) {
          cat("\n","unknown model.name for this model.type","\n")
          return(object=NULL) } } # end of if mean only

# Checking method 
   if ((method!="Nelder-Mead") & (method!="BFGS")) {
       cat("\n","unknown function optim method","\n")
       return(object=NULL) }

# Checking that data is data.frame or list.
   if ((is.data.frame(data)==FALSE) & (is.list(data)==FALSE)) { 
      cat("\n","Input data is neither data frame nor list.","\n")
      return(object=NULL) } # end of check data.frame or list

    cl <- match.call()

# Checking for correct link functions
   if ((link!="log")) {
       cat("\n","unknown link function","\n")
       return(object=NULL) 
           } else {
       attr(link, which="mean") <- make.link(link) } #end of if link

# as in betareg package handling subset=option 
    if (missing(data)) { data <- environment(formula) }
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE 

    FBoth <- Formula(formula) 
    lenFB <- length(FBoth) 
    mf[[1L]] <- as.name("model.frame")

# name of response variable
    wk.name <- attr(FBoth, which="lhs") 
    if (length(wk.name)>1) {
       cat("\n","more than one variable name on lhs of the formula","\n")
       return(object=NULL) }

# data.frame or list

   if (is.data.frame(data)==TRUE) { 

# indicator of whether a data frame (TRUE) or a list (FALSE)
# data frame input
      data.type <- TRUE
      mfBoth    <- model.frame(FBoth,data=data)
# Checking formula and data
      resp.var  <- model.part(FBoth,data=mfBoth,lhs=1)
      nvar      <- length(data)
      nobs      <- length(data[[1]])
      if (nvar==1) { covariates  <- NULL
                   } else { covariates <- data } 
      list.data  <- lapply(1:nobs, function(i) 
                          c(rep(0,resp.var[[1]][i]),1) ) # end of lapply
      mean.obs     <- resp.var[[1]]
      variance.obs <- rep(0,nobs)
      scalef.obs   <- rep(1,nobs)

      wkdata <- data.frame(mean.obs,scalef.obs,data)

                                   } else { 

# list of frequency distributions input
      data.type <- FALSE

# Checking for name of dependent variable in list of data
      wk.name <- attr(FBoth, which="lhs") 
      nvar <- length(data)
      nobs <- length(data[[1]])
      if (nvar==1) { covariates  <- NULL }

      if ((nvar==1) & (is.list(data[[1]])==TRUE)) { 
         if ((wk.name==names(data)[1])==TRUE) {
              list.set <- TRUE
              list.data <- data[[1]]
                                       } else {          
              cat("\n","single list with list of data is not named ",wk.name,"\n")
              return(object=NULL) } # end of if wk.name==names(data)[1]
          } else { list.set <- FALSE
                   for ( i in 1:nvar ) { 
            if (is.list(data[[i]])==TRUE) { 
               if ((wk.name==names(data)[i])==TRUE) { 
                  if (list.set==TRUE)  { 
                     cat("\n","More than one list named ",wk.name," within list of data so","\n")
                     cat("not clear which is the dependent variable.","\n")
                                         return(object=NULL)
                                } else { list.data <- data[[i]] 
                                         list.set <- TRUE } # end of list.set==TRUE
                                                    } # end of wk.name==names(data)[i]

                                       } else { 
                            if (i==1) { covariates <- data.frame(data[1])
                                        names(covariates[1]) <- names(data[1])                   
                                             } else { 
                               if ((i==2) & ((list.set)==TRUE)) { 
                                               covariates   <- data.frame(data[2])
                                               names(covariates[1]) <- names(data[2]) 
                                                         } else { 
                                               wks <- length(covariates) + 1
                                               covariates <- data.frame(covariates,data[i])
                                               names(covariates[wks]) <- names(data[i])  
                                                           }}} } # end for loop
                                                 } # end if nvar==1 & is.list

      if (list.set==FALSE) { 
         cat("\n","No list named ",wk.name," within list of data.","\n")
         return(object=NULL) } # end of if list.set==FALSE

      mean.obs     <- rep(0,nobs)
      variance.obs <- rep(0,nobs)
      scalef.obs   <- rep(1,nobs)
      for ( i in 1:nobs ) { 
         count  <- list.data[[i]]
         ncount <- sum(count)  
         nmax1  <- length(count)
         nmax   <- nmax1 - 1
         cnum   <- 0:nmax
         if (ncount==1) { mean.obs[i] <- t(cnum)%*%count
                        } else {
            mean.obs[i]     <- t(cnum)%*%count/ncount 
            variance.obs[i] <- (t(cnum*cnum)%*%count - 
                            ncount*mean.obs[i]*mean.obs[i]) / (ncount - 1)
            if (mean.obs[i]>0) {
               scalef.obs[i]   <- variance.obs[i]/mean.obs[i] }
# Checking for NA for variance.obs and scalef.obs introduced Version 1.01
            if (is.na(variance.obs[i])==TRUE) { variance.obs[i] <- 0 }
            if (is.na(scalef.obs[i])==TRUE)   { scalef.obs[i] <- 1 }
                              } } # end of for loop

# Checking for no covariates 
      if (is.null(covariates)==TRUE) { wkdata <- data.frame(mean.obs,scalef.obs) 
                              } else { wkdata <- data.frame(mean.obs,scalef.obs,covariates) } 

                            }  # end of if is.data.frame

# weights is not null  
   if (is.null(weights)==FALSE) {
      if (is.null(attributes(weights))==TRUE) {
         attr(weights, which="normalize") <- FALSE }
# normalizing weights 
      if (attr(weights, which="normalize")==TRUE)  {
         wkv <- c(rep(0, nobs))
         wkv <- sapply(1:nobs, function(i) { 
                   wkv[i] <- sum(list.data[[i]]) } ) # end of sapply 
         total.n <- sum(wkv)
         if (is.null(attr(weights, which="norm.to.n"))==FALSE) {
               total.n <- as.numeric(attr(weights, which="norm.to.n")) }
         if (data.type==TRUE) {
            weights <- total.n*weights/sum(weights)
                     } else {
            weights <- lapply(1:nobs, function(i) { weights[[i]] <- 
                         (list.data[[i]]>0)*weights[[i]] } ) # end of lapply
            wkv <- sapply(1:nobs, function(i) { wkv[i] <- 
                         sum(as.numeric(list.data[[i]]*weights[[i]])) } ) # end of sapply
            weights <- lapply(1:nobs, function(i) {
               weights[[i]] <- total.n*weights[[i]] / sum(wkv) } ) # end of lapply
                            } } # end if(attr(weights, which="normalize")==TRUE)
				  } # end is.null(weights) 

# Setting up vector vnmax 
    vnmax  <- sapply(list.data,length) - 1

   mf$data <- wkdata
   mf$resp.var <- mean.obs
   mf$formula  <- update(formula(FBoth,lhs=NULL,rhs=1), mean.obs ~ . )
# Setting up weights
   if (data.type==FALSE) { mf$weights=vnmax } # end if data.type

# This statement evaluates mf as a data.frame
   wkdata <- eval(mf, parent.frame())

#   evaluating scale-factor model if used 
    if (model.type=="mean and scale-factor") { 
          mf.scalef <- mf
          mf.scalef$resp.var <- scalef.obs
       if (lenFB[2]==2) {
             mf.scalef$formula <- update(formula(FBoth,lhs=NULL,rhs=2), scalef.obs ~ . )
                        } else {
             mf.scalef$formula <- update(formula(FBoth,lhs=NULL,rhs=NULL), scalef.obs ~ 1 )
                               } # end if (lenFB[2]==2)

# This statement evaluates mf.scalef as a data.frame
      temp.wkdata <- eval(mf.scalef, parent.frame())
# removing from temp.wkdata the variables common to both data frames, i.e., it & wkdata
      intersection.var <- intersect(names(wkdata), names(temp.wkdata))
      dup.var <- duplicated(c(intersection.var, names(temp.wkdata)))
      wks <- length(intersection.var) + 1
      wke <- length(dup.var)
      dup.var <- dup.var[wks:wke]
      temp.wkdata <- subset(temp.wkdata, select=(dup.var==FALSE))

# merging two data sets i.e., that for mean.obs and that for scalef.obs
      wkdata <- merge(wkdata, temp.wkdata, by=0, incomparables=TRUE, sort=FALSE)
                    } # end of if mean and scale-factor

# removing the variable (resp.var) and
# removing the variable Row.names added in by the merge operation
# when model type is mean and scale-factor 
   wkdata <- subset(wkdata, select=((names(wkdata)!="(resp.var)") &
                                    (names(wkdata)!="Row.names")))

# subsetting mean.obs, variance.obs, scalef.obs if subset in operation
   if (is.null(subset)==FALSE) { 
      mean.obs     <- mean.obs[subset]
      variance.obs <- variance.obs[subset] 
      scalef.obs   <- scalef.obs[subset] 
      weights      <- weights[subset] 
# list.data also needs to be subsetted 
      list.data <- list.data[subset]
# following inserted 7th July 2017 to handle an error resulting when list input is used
      if (data.type==TRUE) { resp.var     <- resp.var[subset] 
         n.var        <- n.var[subset] } } # end of is.null(subset)

# if data frame adding in the variable resp.var to the data frame
# for use calculating initial estimates with glm
      if (data.type==TRUE) { wkdata <- data.frame(wkdata, resp.var) }

# Resetting nobs if subset is being used
      if (is.null(subset)==FALSE) { nobs <- length(wkdata[[1]]) }

# attaching wkdata to search path 
#      attach(wkdata, warn.conflicts=FALSE)

# Resetting nobs if subset is being used
      if (is.null(subset)==FALSE) { nobs <- length(wkdata[[1]])
                                    list.data <- list.data[subset] } # end of is.null(subset)

# Checking if offset.var involved and if it is adding it to workdata so that 
# both variables offset(****) and **** are in the data frame in order to use offsets in glm 
# N.B. DO NOT USE THE NAME offset.*** as an offset name in the data sets input
    offset.var <- NULL
    wks <- length(wkdata)
    for ( i in 1:wks) { nchar <- nchar(names(wkdata[i]))
       if ((nchar>7) & (((substr(names(wkdata[i]), 1, 7)=="offset(")==TRUE) | 
                         (substr(names(wkdata[i]), 1, 7)=="offset.")==TRUE)) { 
             offset.var <- wkdata[i]
             names(offset.var) <- substr(names(wkdata[i]), 8, (nchar-1))
             wkdata <- data.frame(wkdata, offset.var)
                       } } # end of for i

# Changing offsets names from offset.xxx. to offset(xxx)
    wks <- length(wkdata)
    names(wkdata) <- sapply(1:wks, function(i) {
        nchar <- nchar(names(wkdata[i]))
        if ((nchar>7) & ((substr(names(wkdata[i]), 1, 7)=="offset.")==TRUE)) { 
              new.label <- paste("offset(", substr(names(wkdata[i]), 8, (nchar-1)), ")", sep="") 
              names(wkdata[i]) <- new.label
                    } else { 
              names(wkdata[i]) <- names(wkdata[i]) }
                        } ) # end of sapply

# removing wkdata on exit 
    on.exit(rm(wkdata))

# Checking arguments of function

# Checking for correct model.types 
   if ((model.type!="mean only") & (model.type!="mean and scale-factor")) {
      cat("\n","unknown model.type","\n")
      return(object=NULL) }
   if (model.type=="mean only") {
      if (lenFB[2]==2) { 
         cat("\n","model.type is mean only but rhs of formula has two parts to it")
         cat("\n","2nd part of rhs of formula is ignored","\n") } # end if lenFB[2]
         FBoth <- update(FBoth,mean.obs ~ .) } # end of mean only
   if (model.type=="mean and scale-factor") {
      if (lenFB[2]==1) { 
            cat("\n","Model for scale-factor set to intercept only.","\n")
            FBoth <- update(FBoth,mean.obs | scalef.obs ~ . | 1 )
                                  } else {
            FBoth <- update(FBoth,mean.obs | scalef.obs ~ . | . )
                       } } # end of if mean and scale-factor 

   terms.mean <- terms(formula(FBoth,lhs=1,rhs=1))
   if (model.type=="mean only") { 
      terms.scalef <- terms(formula( 0 ~ 1 )) 
                             } else { 
      if (lenFB[2]==1) { terms.scalef <- terms(formula( 0 ~ 1 ))
                } else { terms.scalef <- terms(formula(FBoth,lhs=1,rhs=2)) } } # end of if model.type
   terms.full   <- terms(formula(FBoth))

# Checking that list.data is a list  
   if (is.list(list.data)==FALSE) { cat("\n","list.data is not a list","\n")
      return(object=NULL) }

# extracting offsets and model matrices 
   mf.mean <- model.frame(formula(FBoth,lhs=1,rhs=1), data=wkdata)

# The following command seems to have problems with handling a variable
# named class. This was the original name of the variable in the Titanic
# data. When changed to pclass there were no problems.
   covariates.matrix.mean <- model.matrix(mf.mean, data=wkdata) 
   offset.mean <- model.offset(mf.mean)
   covariates.matrix.scalef <- matrix(c(rep(1,nrow(covariates.matrix.mean))),ncol=1)
   offset.scalef <- NULL
   if (model.type=="mean and scale-factor") {
       if (lenFB[2]==2) {
         mf.scalef <- model.frame(formula(FBoth,lhs=2,rhs=2), data=wkdata)
         covariates.matrix.scalef <- model.matrix(mf.scalef, data=wkdata)
         offset.scalef <- model.offset(mf.scalef) } } # end if mean and scale-factor 

   if (is.null(offset.mean)==TRUE) { offset.mean <- c(rep(0,nobs)) }
   if (is.null(offset.scalef)==TRUE) { offset.scalef <- c(rep(0,nobs)) } 

# Setting up initial estimates of parameters if not input 
   if (is.null(initial)==TRUE) {
#      nmax1   <- length(mean.obs)
# data.frame input or case data input as list, 
# switching off warnings from glm usually caused by non integer values for the counts
      options(warn=-1)
# setting up new formula in standard form for data.frame input
   if (data.type==TRUE) { 
# changing to standard form of arguments to glm for Poisson for data as a data frame
      if (is.null(weights)==TRUE) { 
         glm.Results <- glm(formula(FBoth,lhs=1,rhs=1),family=poisson(link=link),
                            data=wkdata) 
                           } else { 
         wkdata <- data.frame(wkdata, weights) 
         glm.Results <- glm(formula(FBoth,lhs=1,rhs=1),family=poisson(link=link),
                            data=wkdata, weights=weights) }
                     } else {
         weights.mean <- c(rep(1,nobs))
         if (is.null(weights)==FALSE) { 
            weights.mean <- as.vector(sapply(1:nobs, function(i) { 
                            weights.mean[i] <- t(weights[[i]])%*%list.data[[i]] } ) ) } # end of is.null(weights)
         weights.mean <- weights.mean*vnmax
         wkdata <- data.frame(wkdata, weights.mean) 
         glm.Results <- glm(formula(FBoth,lhs=1,rhs=1),family=poisson(link=link),
                            data=wkdata, weights=weights.mean)  
                     } # end if data.type
# switching warnings from glm back on
      options(warn=0)                      

      initial.mean        <- coefficients(glm.Results)
      names(initial.mean) <- names(coefficients(glm.Results))

      if (model.type=="mean and scale-factor") {
         glm.Results <- glm(formula(FBoth,lhs=2,rhs=2), family=gaussian(link="log"),
                            subset=(scalef.obs>0), data=wkdata)
         initial.scalef <- coefficients(glm.Results)
         names(initial.scalef) <- names(coefficients(glm.Results))
         if (model.name=="general")  { parameter <- c(initial.mean,initial.scalef,0) 
            names(parameter) <- c(names(initial.mean),names(initial.scalef),"log(b)") }
         if ((model.name=="general fixed b") | (model.name=="limiting")) { 
            parameter <- c(initial.mean,initial.scalef) 
            names(parameter) <- c(names(initial.mean),names(initial.scalef)) }
# Calculating log likelihood for initial estimates for regression of log(scale-factor) (mean)
# for mean and scale-factor model 
         loglikelihood <- LL.Regression.Counts(parameter,model.type,model.name,
               link=link,list.data,covariates.matrix.mean,covariates.matrix.scalef,
               offset.mean,offset.scalef,ltvalue,utvalue,fixed.b,weights,grad.method)
# Calculating log likelihood for initial estimates for log(scale-factor) of 0 for
# the variance model 
         wks <- length(initial.scalef)
         if (model.name=="general")  { wk.parameter <- c(initial.mean,rep(0,wks),0) 
            names(wk.parameter) <- c(names(initial.mean),names(initial.scalef),"log(b)") }
         if ((model.name=="general fixed b") | (model.name=="limiting")) { 
            wk.parameter <- c(initial.mean,rep(0,wks)) 
            names(wk.parameter) <- c(names(initial.mean),names(initial.scalef)) }
         wk.loglikelihood <- LL.Regression.Counts(wk.parameter,model.type,model.name,
               link=link,list.data,covariates.matrix.mean,covariates.matrix.scalef,
               offset.mean,offset.scalef,ltvalue,utvalue,fixed.b,weights,grad.method)
# Use the parameter estimates from the initial one with largest log likelihood
         if (wk.loglikelihood>loglikelihood) { parameter <- wk.parameter }
                                           } # end of if mean and scale-factor

      if (model.type=="mean only") { 
          if ((model.name=="Poisson") | (model.name=="negative binomial fixed b")) {
                                 parameter        <- initial.mean
                                 names(parameter) <- names(initial.mean) }
          if (model.name=="negative binomial") { parameter <- c(initial.mean,0) 
                           names(parameter) <- c(names(initial.mean),"log(b)") }
          if (model.name=="Faddy distribution") { parameter <- c(initial.mean,0,0) 
                           names(parameter) <- c(names(initial.mean),"c","log(b)") }
          if (model.name=="Faddy distribution fixed b") { parameter <- c(initial.mean,0) 
                           names(parameter) <- c(names(initial.mean),"c") }
                                   } # of if model.type mean only

                                   } else { # else of initial

# Checking length of input initial against model value
      npar.one <- ncol(covariates.matrix.mean)
      numpar <- length(initial) 
      if (model.type=="mean only") { npar <- npar.one 
              loc.c <- npar + 1 
              if ((model.name=="negative binomial") | 
                  (model.name=="Faddy distribution fixed b")) { npar <- loc.c }
              if (model.name=="Faddy distribution") { npar <- npar + 2 } }
      if (model.type=="mean and scale-factor") { 
              npar.two <- ncol(covariates.matrix.scalef)
              npar <- npar.one + npar.two 
              if (model.name=="general")  { npar <- npar + 1 } }
      if (numpar!=npar) { cat("\n","number of parameters error","\n")
         return(object=NULL) }
      parameter <- initial
      if (is.null(names(initial))==TRUE) { 
         cat("\n","WARNING: initial has no associated names","\n")
         names(parameter) <- 1:numpar } # end if is.null(names)
                           } # end if is.null(initial)
   start <- parameter

# Checking nobs is the same value as length(list.data)
   wks <- length(list.data)
   if (wks!=nobs) { 
         cat("\n","number of rows of covariates.matrix.mean not equal to length of list.data","\n")
         return(object=NULL) }

# Checking for improper probability distributions for initial estimates
# Setting up vector vnmax and calculating likelihood for initial estimates
#    vnmax  <- sapply(list.data,length) - 1
    output <- Model.Counts(parameter,model.type,model.name,link=link,covariates.matrix.mean,
              covariates.matrix.scalef,offset.mean,offset.scalef,
              fixed.b,vnmax)
    for ( i in 1:nobs) { probability <- output$probabilities[[i]]
       nmax1 <- vnmax[i] + 1
       wks   <- sum(probability)
       if (is.finite(wks)==TRUE) {
          wks <- round(sum(probability),digits=10)
# rounding error can cause the sum of the probabilities to <0 or >1
# so rounded to 10 decimal places                                        
          if ((wks<0) | (wks>1))  { cat("\n","improper distribution produced by initial estimates","\n")
             return(object=NULL) }
          wks <- sum((round(probability,digits=10)<0) | 
                     (round(probability,digits=10)>1))
          if (wks>0) { cat("\n","improper distribution produced by initial estimates","\n")
              return(object=NULL) }
                                 } else { cat("\n","improper distribution produced by initial estimates","\n")
                                          return(object=NULL) } # end of if is.finite(wks)
                    } # end of for loop

# Checking that fixed.b has a positive value if that model is used
    if (((model.name=="general fixed b") | 
         (model.name=="negative binomial fixed b distribution fixed b") |
         (model.name=="Faddy distribution fixed b")) & 
        ((is.na(fixed.b)==TRUE) | (fixed.b<=0))) { cat("\n","value of fixed.b is NA or <=0","\n")
              return(object=NULL) }

# Checking that the input initial estimate for c is <1 if a 
# Faddy distribution model is being fitted.
    if ((is.null(initial)==FALSE) & ((model.name=="Faddy distribution") | 
        (model.name=="Faddy distribution fixed b"))) {
       if (round(parameter[loc.c],digits=4)>=1) { cat("\n","initial c>=1 is not allowed for the Faddy distribution","\n")
              return(object=NULL) }}

##########################################################################
# start of main calculations 
##########################################################################

# link function for mean
      npar.mean   <- ncol(covariates.matrix.mean)
      npar <- npar.mean
#      r.parameter.mean <- rep(0,npar.mean) 
      r.parameter.mean <- parameter[1:npar.mean] 
      lp.mean <- covariates.matrix.mean%*%r.parameter.mean + offset.mean
# link function for scale-factor 
      if (model.type=="mean and scale-factor") { 
         npar.scalef <- ncol(covariates.matrix.scalef)
         npar <- npar.mean + npar.scalef
#         r.parameter.scalef <- rep(0,npar.scalef) 
         wks <- npar.mean + 1
         r.parameter.scalef <- parameter[wks:npar] 
         lp.scalef <- covariates.matrix.scalef%*%r.parameter.scalef + offset.scalef } # end of if statement

# Checking grad.method if method is BFGS
      if (method=="BFGS") { 
         if (is.null(attr(method,which="grad.method"))==TRUE) { 
               attr(method,which="grad.method") <- "simple" }
         if ((attr(method,which="grad.method")!="simple") & 
             (attr(method,which="grad.method")!="Richardson")) {
            cat("\n","unknown gradient method with method BFGS so reset to simple","\n")
               attr(method,which="grad.method") <- "simple" }
               grad.method <- attr(method,which="grad.method")
                           } else {
               grad.method <- NULL
            } # end if method=BFGS

# Setting defaults and checking control parameters for optim
   if (is.null(control)==TRUE) { 
#      control=list(fnscale=-1,trace=0,maxit=1000,abstol=1e-8,reltol=1e-8,
#                   alpha=1.0,beta=0.5,gamma=2.0,REPORT=10) 
      control=list(fnscale=-1, trace=0) 
                               } else { 
# Control parameters related to Nelder-Mead except for REPORT which
# is for BFGS. Those related only other methods of optim are considered to 
# be unknown and ignored.
     ncontrol <- length(control)
#     wk.control <- list(fnscale=-1,trace=0,maxit=1000,abstol=1e-8,reltol=1e-8,
#                        alpha=1.0,beta=0.5,gamma=2.0,REPORT=10) 
     wk.control=list(fnscale=-1, trace=0) 
     for ( i in 1:ncontrol) { 
         if ((names(control[i])!="fnscale") & (names(control[i])!="trace") & 
             (names(control[i])!="maxit") & (names(control[i])!="abstol") & 
             (names(control[i])!="reltol") & (names(control[i])!="alpha") & 
             (names(control[i])!="beta") & (names(control[i])!="gamma") & 
             (names(control[i])!="REPORT")) { 
            wkname <- names(control[i])
            cat("\n","WARNING: Argument in control is unknown so ignored.","\n")
                                  } # end of if
         if (names(control[i])=="fnscale") { wk.control$fnscale <- control[[i]] }
         if (names(control[i])=="trace")   { wk.control$trace   <- control[[i]] }
         if (names(control[i])=="maxit")   { wk.control$maxit   <- control[[i]] }
         if (names(control[i])=="abstol")  { wk.control$abstol  <- control[[i]] }
         if (names(control[i])=="reltol")  { wk.control$reltol  <- control[[i]] }
         if (names(control[i])=="alpha")   { wk.control$alpha   <- control[[i]] }
         if (names(control[i])=="beta")    { wk.control$beta    <- control[[i]] }
         if (names(control[i])=="gamma")   { wk.control$gamma   <- control[[i]] }
         if (names(control[i])=="REPORT")  { wk.control$REPORT  <- control[[i]] }
                             } # end of for loop
     control <- wk.control
                                   } # end of is.null(control)

# Fitting model using optim, options BFGS with gr=NULL, BFGS with gr=gradient function

       if (length(parameter)==1) {
# Using glm to obtain estimate of mean for a Poisson model
# switching off warnings from glm usually caused by non integer values for the counts
          options(warn=-1)
# setting up new formula in standard form for data.frame input
          if (data.type==TRUE) { 

# changing to standard form of arguments to glm for Poisson for data as a data frame
         if (is.null(weights)==TRUE) { 
            FBoth_one <- update(FBoth, mean.obs ~ .) 
            glm.Results <- glm(formula(FBoth_one,lhs=1,rhs=1),family=poisson(link=link),
                               data=wkdata)
                           } else { 
            wkdata <- data.frame(wkdata, weights) 
            FBoth_one <- update(FBoth, mean.obs ~ .) 
            glm.Results <- glm(formula(FBoth_one,lhs=1,rhs=1),family=poisson(link=link),
                               data=wkdata, weights=weights) }
                        } else {
            weights.mean <- c(rep(1,nobs))
            if (is.null(weights)==FALSE) { 
               weights.mean <- as.vector(sapply(1:nobs, function(i) { 
                  weights.mean[i] <- t(weights[[i]])%*%list.data[[i]] } ) ) } # end of is.null(weights)
            weights.mean <- weights.mean*vnmax
            wkdata <- data.frame(wkdata, weights.mean) 
            FBoth_one <- update(FBoth, mean.obs ~ . ) 
            glm.Results <- glm(formula(FBoth_one,lhs=1,rhs=1),family=poisson(link=link),
                               data=wkdata, weights=weights.mean) 
                        } # end if data.type
          estimates <- coefficients(glm.Results)
# switching warnings from glm back on
          options(warn=0)
          wks <- LL.Regression.Counts(parameter,model.type,model.name,link,list.data,
                   covariates.matrix.mean,covariates.matrix.scalef,offset.mean,offset.scalef,
                   ltvalue,utvalue,fixed.b,weights,grad.method)
          converged   <- FALSE 
          attr(converged, which="code") <- 0        
          iterations <- glm.Results$iter
          wk.optim <- list(par=glm.Results$coefficients,value=wks,counts=c(glm.Results$iter, NA),
                           convergence=0,message=NULL)
                     } else {
          if (method=="Nelder-Mead") { 
              wk.method <- "Nelder-Mead" } else {
              wk.method <- "BFGS" } # end of if method       
          wk.optim <- optim(parameter,fn=LL.Regression.Counts,gr=LL.gradient,
                           model.type,model.name,link=link,list.data,
                           covariates.matrix.mean,covariates.matrix.scalef,
                           offset.mean,offset.scalef,ltvalue,utvalue,
                           fixed.b,weights,grad.method,
                           method=wk.method,hessian=FALSE,control=control) 
          if (wk.optim$convergence==0) { converged <- TRUE 
                               } else { converged <- FALSE } 
          attr(converged, which="code") <- wk.optim$convergence                          
          if (is.null(wk.optim$message)==FALSE) {
             cat(" message        ",wk.optim$message,"\n") }
          iterations <- wk.optim$counts[1]
          estimates  <- wk.optim$par
                                     } # end of if length(parameter)=1 
          names(estimates) <- names(parameter) 

          npar <- length(wk.optim$par)
          nobs <- nrow(covariates.matrix.mean) 
          mean.par       <- rep(0,nobs)
          variance.par   <- rep(1,nobs)
          variance.limit <- rep(0,nobs)
# Calculation of p and means from parameter estimates and design matrices
          npar.mean     <- ncol(covariates.matrix.mean)
          r.parameter.mean <- rep(0,npar.mean) 
          r.parameter.mean <- wk.optim$par[1:npar.mean] 

# inverse of link function
          mean.par <- attr(link, which="mean")$linkinv(lp.mean)
          if (model.type=="mean only") { 
             if ((model.name=="Poisson") | (model.name=="negative binomial fixed b")) { wkv.coefficients <- 
                list(mean.est=wk.optim$par[1:npar.mean], scalef.est=NULL)
                                 } else { wks <- npar.mean + 1
                                          wkv.coefficients <- 
                list(mean.est=wk.optim$par[1:npar.mean], scalef.est=wk.optim$par[wks:npar]) }
                                } else {
# Calculation of scale-factors and variances from parameter estimates and design matrices
             npar.scalef <- ncol(covariates.matrix.scalef)
             wks <- npar.mean + 1
             wkv.coefficients <- list(mean.est=wk.optim$par[1:npar.mean], scalef.est=wk.optim$par[wks:npar]) 
# inverse link for scale-factor 
             scalef.par <- exp(lp.scalef)
                                              } # if model.type
# model.hessian from hessian from numDeriv method Richardson
          model.hessian <- hessian(LL.Regression.Counts,x=wk.optim$par,
                        method="Richardson",method.args=list(r=6,eps=1.e-4),
                        model.type=model.type,model.name=model.name,link=link,
                        list.data=list.data,
                        covariates.matrix.mean=covariates.matrix.mean,
                        covariates.matrix.scalef=covariates.matrix.scalef,
                        offset.mean=offset.mean,offset.scalef=offset.scalef,
                        ltvalue=ltvalue,utvalue=utvalue,fixed.b=fixed.b,
                        weights=weights,grad.method=grad.method)

# Deleting appropriate row and column of matrix for Faddy distribution if c=1
          if ((model.name=="Faddy distribution") | 
              (model.name=="Faddy distribution fixed b")) { nparm1 <- npar - 1
                       nparm1sq <- nparm1*nparm1
                       nparm2   <- nparm1 - 1
                       wk.hessian <- matrix(c(rep(0,nparm1sq)),ncol=nparm1)
              if (model.name=="Faddy distribution") { wk.c <- round(estimates[nparm1],digits=7)
                                             } else { wk.c <- round(estimates[npar],digits=7) }
              if ((model.name=="Faddy distribution") & (wk.c==1)) { 
                       wk.hessian[1:nparm2,1:nparm2] <- model.hessian[1:nparm2,1:nparm2]
                       wk.hessian[nparm1,] <- c(model.hessian[npar,1:nparm2],model.hessian[npar,npar])
                       wk.hessian[,nparm1] <- c(model.hessian[1:nparm2,npar],model.hessian[npar,npar])
                       model.hessian       <- wk.hessian
                                          } # if Faddy & c=1
              if ((model.name=="Faddy distribution fixed b") & (wk.c==1)) { 
                       wk.hessian[1:nparm1,1:nparm1] <- model.hessian[1:nparm1,1:nparm1]
                       model.hessian                 <- wk.hessian
                                          } # if Faddy fixed b & c=1
                                                     } # end first if

# checking condition number of vcov and inverting plus square rooting to get variance/covariance matrix
          deter <- det(model.hessian)
          wk.npar <- nrow(model.hessian)
          if ((is.finite(deter)==FALSE) | (deter==0)) { vcov <- matrix(c(rep(NA,(wk.npar*wk.npar))),ncol=wk.npar)
                        } else { if (wk.npar==1) { vcov <- -1/model.hessian
                                          } else { condition <- rcond(model.hessian)
# function rcond is from the package Matrix and gives the (reciprocal) condition number
# near 0 is ill-conditioned, near 1 is well conditioned
                                    if (condition>1e-16) {
                                        vcov <- - solve(model.hessian)
                                           } else { vcov <- matrix(c(rep(NA,(wk.npar*wk.npar))),ncol=wk.npar) }}} 
# Setting up row and column names for vcov
          colnames(vcov) <- rownames(vcov) <- names(wk.optim$par[1:wk.npar]) 
          output.model <- Model.Counts(wk.optim$par,model.type,model.name,link,covariates.matrix.mean,
                                 covariates.matrix.scalef,offset.mean,offset.scalef,
                                 fixed.b,vnmax)

# Calculate the fitted value vector for the fitted model parameters
          mean.prob <- rep(0, nobs)
          vone <- rep(1,nobs)

          if (model.name=="limiting") {
             valpha <- output.model$FDparameters$out.valpha
             vbeta  <- output.model$FDparameters$out.vbeta
                                } else {
             va <- output.model$FDparameters$out.va
             vb <- output.model$FDparameters$out.vb
             vc <- output.model$FDparameters$out.vc
             v1mc <- vone - output.model$FDparameters$out.vc } # end if "limiting"
          if ((model.name=="Poisson") | (model.name=="negative binomial") | 
              (model.name=="negative binomial fixed b")) {
             if (model.name=="Poisson") { mean.prob <- va 
                                 } else { 
                 mean.prob <- vb*(attr(link,which="mean")$linkinv(va) - vone) }
                                                  } else { 
                 if (model.name=="limiting") {
                     mean.prob <- sapply(1:nobs, function(j) 
                        mean.prob[j] <- - log(1 - valpha[j]*vbeta[j]) / vbeta[j] )
                                } else {
                     mean.prob <- sapply(1:nobs, function(j) 
                        if ((abs(v1mc[j])<1.e-6)==TRUE) { 
                           mean.prob[j] <- vb[j]*(attr(link,which="mean")$linkinv(va[j]) - 1) 
                                                 } else {  
                           mean.prob[j] <- vb[j]*((1+va[j]*v1mc[j]/(vb[j]**v1mc[j]))**(1/v1mc[j])-1) } )
                                                          } } # end if "Poisson" | "negative binomial"

# For the list data type (frequency distribution of data) the degrees of freedom
# for null and residual need to be set to those of the data frame data type

          if (data.type==TRUE) { total.ninlist <- nobs
                               } else {
             vninlist <- c(rep(0,length(list.data)))
             vninlist <- sapply(1:length(list.data), function(ilist) 
                               vninlist[ilist] <- sum(list.data[[ilist]]) )
             total.ninlist <- sum(vninlist)
                                      } # end of is.data.frame

          object <- list(data.type=data.type, list.data=list.data, call=cl,
                       formula=formula, model.type=model.type, model.name=model.name,
                       link=link,covariates.matrix.mean=covariates.matrix.mean,
                       covariates.matrix.scalef=covariates.matrix.scalef,
                       offset.mean=offset.mean,offset.scalef=offset.scalef,
                       ltvalue=ltvalue,utvalue=utvalue,fixed.b=fixed.b,
                       coefficients=wkv.coefficients,loglik=wk.optim$value,vcov=vcov,
                       n=nobs, nobs=nobs, df.null=total.ninlist, df.residual=(total.ninlist-length(wk.optim$par)),
                       vnmax=vnmax, weights=weights,converged=converged, method=method,
                       start=start, optim=wk.optim, control=control,
                       fitted.values=mean.prob, y=mean.obs, 
                       terms=list(mean=terms.mean,scale.factor=terms.scalef,full=terms.full)) 
     attr(object, "class") <- c("CountsEPPM")
     return(object)     }

waldtest.CountsEPPM <-
function(object, ..., vcov = NULL, test = c("Chisq", "F"))
{
  ## methods needed:
  ## - terms()
  ## - update()
  ## - formula()
  ## - nobs() or residuals() -> only for determining number of observations
  ## - df.residual() or logLik
  ## - coef() -> needs to be named, matching names in terms() and vcov()
  ## - vcov(), potentially user-supplied

  ## use S4 methods if loaded
  coef0   <- function(x, ...) {
#    coef1 <- if("stats4" %in% loadedNamespaces()) stats4::coef else coef
    coef1 <-  coef
    na.omit(coef1(x, ...))
  } # end of use S4 methods if loaded


#  logLik0 <- if("stats4" %in% loadedNamespaces()) stats4::logLik else logLik
#  update0 <- if("stats4" %in% loadedNamespaces()) stats4::update else update
  logLik0 <- logLik
  update0 <- update
  nobs0   <- function(x, ...) {
#    nobs1 <- if("stats4" %in% loadedNamespaces()) stats4::nobs else nobs
    nobs1 <- nobs
    nobs2 <- function(x, ...) NROW(residuals(x, ...))
    rval <- try(nobs1(x, ...), silent = TRUE)
    if(inherits(rval, "try-error") | is.null(rval)) rval <- nobs2(x, ...)
    return(rval)
  }

  vcov0   <- if(!is.null(vcov)) vcov else stats::vcov
#  vcov0   <- if(!is.null(vcov)) vcov else {
#    if("stats4" %in% loadedNamespaces()) stats4::vcov else stats::vcov
#  }
  df.residual0 <- function(x) {
    df <- try(df.residual(x), silent = TRUE)
    if(inherits(df, "try-error") | is.null(df)) df <- try(nobs0(x) - attr(logLik0(x), "df"), silent = TRUE)
    if(inherits(df, "try-error") | is.null(df)) df <- try(nobs0(x) - length(as.vector(coef0(x))), silent = TRUE)
    if(inherits(df, "try-error")) df <- NULL
    return(df)
  }

  ## model class
  cls <- class(object)[1]

  ## convenience functions:
  ## 1. extracts term labels
  tlab <- function(x) {
    tt <- try(terms(x), silent = TRUE)
    if(inherits(tt, "try-error")) "" else attr(tt, "term.labels")
  } # end of extracts term labels

  ## 2. extracts model name
#  if(is.null(name)) name <- function(x) {
    name <- function(x) {
    rval <- try(formula(x), silent = TRUE)
    if(inherits(rval, "try-error") | is.null(rval)) rval <- try(x$call, silent = TRUE)
    if(inherits(rval, "try-error") | is.null(rval)) return(NULL) else return(paste(deparse(rval), collapse="\n"))
  } # end of extracts model name

  ## 3. compute an updated model object
  modelUpdate <- function(fm, update) {
    ## if `update' is numeric or character, then assume that the 
    ## corresponding variables (i.e., terms) are redundant (i.e., should be omitted)
    if(is.numeric(update)) {
      ## sanity checking of numeric update specification
      if(any(update < 1)) {
        warning("for numeric model specifications all values have to be >=1")
	update <- abs(update)[abs(update) > 0]
      }
      if(any(update > length(tlab(fm)))) {
        warning(paste("more terms specified than existent in the model:",
	        paste(as.character(update[update > length(tlab(fm))]), collapse = ", ")))
	update <- update[update <= length(tlab(fm))]
      }
      ## finally turn numeric into character update specification
      update <- tlab(fm)[update]
    }
    if(is.character(update)) {
      ## sanity checking of character update specification
      if(!all(update %in% tlab(fm))) {
        warning(paste("terms specified that are not in the model:",
	        paste(dQuote(update[!(update %in% tlab(fm))]), collapse = ", ")))
        update <- update[update %in% tlab(fm)]
      }
      if(length(update) < 1) stop("empty model specification")  
      ## finally turn character into formula update specification       
      update <- as.formula(paste(". ~ . -", paste(update, collapse = " - ")))
    }
    if(inherits(update, "formula")) update <- update0(fm, update)
    if(!inherits(update, cls)) stop(paste("original model was of class \"", cls,
      "\", updated model is of class \"", class(update)[1], "\"", sep = ""))
    return(update)
  } # end of modelUpdate 

  ## 4. compare two fitted model objects
  modelCompare <- function(fm, fm.up, vfun = NULL) {
    q <- length(coef0(fm)) - length(coef0(fm.up))

    if(q > 0) {
      fm0 <- fm.up
      fm1 <- fm
    } else {
      fm0 <- fm
      fm1 <- fm.up
    }

    k <- length(coef0(fm1))
    n <- nobs0(fm1)

    llcoeffm0 <- coef0(fm0)
    llcoeffm1 <- coef0(fm1)

# Changed section to handle having two sets of coefficient names, one for each 
# of p and scale factor, which could include the same names of coefficients.
    if(!all(tlab(fm0) %in% tlab(fm1))) stop("models are not nested")

    if ((fm0$model.type=="p only") & (fm1$model.type=="p only")) {

    ## determine omitted variables
    ovar <- which(!(names(llcoeffm1) %in% names(llcoeffm0)))

                                          } else {

    ## determine omitted variables

    fm0.p.names <- names(fm0$coefficients$p.est)
    fm0.scalef.names <- names(fm0$coefficients$scalef.est)
    fm1.p.names <- names(fm1$coefficients$p.est)
    fm1.scalef.names <- names(fm1$coefficients$scalef.est)
    wk.fm0.names <- fm0.p.names
    for (i in 1:length(fm0.p.names)) {
       wk.fm0.names[i] <- paste("p", fm0.p.names[i], sep=" ") } # end of for

    if (is.null(fm0.scalef.names)==FALSE) {
       wk.fm0.scalef.names <- fm0.scalef.names
       for(i in 1:length(fm0.scalef.names)) {
          wk.fm0.scalef.names[i] <- paste("scalef", fm0.scalef.names[i], sep=" ") } # end of for 
       wk.fm0.names <- c(wk.fm0.names,wk.fm0.scalef.names) } # end of is.null

    wk.fm1.names <- fm1.p.names
    for (i in 1:length(fm1.p.names)) {
       wk.fm1.names[i] <- paste("p", fm1.p.names[i], sep=" ") } # end of for 

    if (is.null(fm1.scalef.names)==FALSE) {
       wk.fm1.scalef.names <- fm1.scalef.names
       for(i in 1:length(fm1.scalef.names)) {
          wk.fm1.scalef.names[i] <- paste("scalef", fm1.scalef.names[i], sep=" ") } # end of for 
        wk.fm1.names <- c(wk.fm1.names,wk.fm1.scalef.names) } # end of is.null

    ovar <- which(!(wk.fm1.names %in% wk.fm0.names)) } # end if ((fm0$model.type

# end of changed section

    ## get covariance matrix estimate
    vc <- if(is.null(vfun)) vcov(fm1)
          else if(is.function(vfun)) vfun(fm1)
	  else vfun

    ## compute Chisq statistic
    stat <- t(llcoeffm1[ovar]) %*% solve(vc[ovar,ovar]) %*% llcoeffm1[ovar]
    return(c(-q, stat))
  } # end of modelCompare

  ## recursively fit all objects (if necessary)
  objects <- list(object, ...)
  nmodels <- length(objects)
  if(nmodels < 2) {
    objects <- c(objects, . ~ 1)
    nmodels <- 2
  }
  
  # remember which models are already fitted and which are described
  # by an update mechanism
  no.update <- sapply(objects, function(obj) inherits(obj, cls))
  
  ## updating
  for(i in 2:nmodels) objects[[i]] <- modelUpdate(objects[[i-1]], objects[[i]])

  ## start of a section changed from betareg code to accommodate having possibly two
  ## sets of responses, one for the p, the other for the scale factor

  ## check responses
  p.getresponse <- function(x) {
    tt <- try(terms(x)$p, silent = TRUE)
    if(inherits(tt, "try-error")) "" else deparse(tt[[2]])
  } # end of p.getresponse
  scalef.getresponse <- function(x) {
    tt <- try(terms(x)$scalef, silent = TRUE)
    if(inherits(tt, "try-error")) "" else deparse(tt[[2]])
  } # end of scalef.getresponse

  responses <- as.character(lapply(objects, p.getresponse))
  sameresp <- responses == responses[1]
  responses <- as.character(lapply(objects, scalef.getresponse))
  sameresp <- c(sameresp, responses == responses[1])

  ## end of changed section

  if(!all(sameresp)) {
    objects <- objects[sameresp]
    warning("models with response ", deparse(responses[!sameresp]),
	    " removed because response differs from ", "model 1")
  }

  ## check sample sizes
  ns <- sapply(objects, nobs0)
  if(any(ns != ns[1])) {
    for(i in 2:nmodels) {
      if(ns[1] != ns[i]) {
        if(no.update[i]) stop("models were not all fitted to the same size of dataset")
	  else {
	    commonobs <- row.names(model.frame(objects[[i]])) %in% row.names(model.frame(objects[[i-1]]))
	    objects[[i]] <- eval(substitute(update(objects[[i]], subset = commonobs),
	      list(commonobs = commonobs)))
	    if(nobs0(objects[[i]]) != ns[1]) stop("models could not be fitted to the same size of dataset")
	  }
      }
    }
  }

  ## check vcov0
  if(nmodels > 2 && !is.null(vcov0) && !is.function(vcov0))
    stop("to compare more than 2 models `vcov' needs to be a function")

  ## setup ANOVA matrix
  test <- match.arg(test)
  rval <- matrix(rep(NA, 4 * nmodels), ncol = 4)
  colnames(rval) <- c("Res.Df", "Df", test, paste("Pr(>", test, ")", sep = ""))
  rownames(rval) <- 1:nmodels
  rval[,1] <- as.numeric(sapply(objects, df.residual0))
  for(i in 2:nmodels) rval[i, 2:3] <- modelCompare(objects[[i-1]], objects[[i]], vfun = vcov0)
  if(test == "Chisq") {
    rval[,4] <- pchisq(rval[,3], round(abs(rval[,2])), lower.tail = FALSE)
  } else {
    df <- rval[,1]
    for(i in 2:nmodels) if(rval[i,2] < 0) df[i] <- rval[i-1,1]
    rval[,3] <- rval[,3]/abs(rval[,2])
    rval[,4] <- pf(rval[,3], abs(rval[,2]), df, lower.tail = FALSE)
  }

  variables <- lapply(objects, name)
  if(any(sapply(variables, is.null))) variables <- lapply(match.call()[-1L], deparse)[1L:nmodels]
  title <- "Wald test\n"
  topnote <- paste("Model ", format(1:nmodels),": ", variables, sep="", collapse="\n")

  structure(as.data.frame(rval), heading = c(title, topnote),
	    class = c("anova", "data.frame"))
}

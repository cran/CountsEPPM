hatvalues.CountsEPPM <-
function (model, ...) {
    x <- model$covariates.matrix.mean
    if (NCOL(x) < 1L) 
        return(structure(rep.int(0, NROW(x)), .Names = rownames(x)))
    if (is.null(model$offset$mean))   { offset.mean <- rep(0, NROW(x)) }
    wts <- weights(model$model)
    if (is.null(wts)) { wts <- rep(1, NROW(x)) }
    beta <- model$coefficients$mean.est[1:ncol(x)]   
    eta <- as.vector(x %*% beta + offset.mean)
    vone <- rep(1,length(eta))
# inverse of link function
    mean <- attr(model$link, which="mean")$linkinv(eta)

# modeling scalefactor

    if (model$model.type=="mean only") { scale.factor <- matrix(rep.int(1, NROW(x)), ncol=1)
                                } else {
       z <- model$covariates.matrix.scalef
       if (is.null(model$offset$scalef)) { offset.scalef <- rep(0, NROW(z)) }
       gamma <- model$coefficients$scalef.est[1:ncol(z)]
       scale.factor_eta <- as.vector(z %*% gamma + offset.scalef) 
       scale.factor <- exp(scale.factor_eta) } # end if model.type

     vvariance <- mean*scale.factor
# first differential of mean w.r.t. eta (linear predictor)
     dmu <- attr(model$link, which="mean")$mu.eta(eta)
# constructing a vector of cases from the list of groups
     if (model$data.type==TRUE) { 
          wk.mean    <- mean
          wk.wts  <- wts
          wk.vvariance <- vvariance
          wk.dmu <- dmu
          wk.x    <- x
                            } else {
       wks <- length(x[1,])
       nend <- 0
       for (ilist in 1:length(model$list.data)) { 
          for ( i in 1:length(model$list.data[[ilist]])) { 
             nt <- model$list.data[[ilist]][i] 
             if (nt>0) { 
                if (nend==0) {
                    wk.x <- as.vector(c(rep(x[ilist,],nt)))
                    wk.mean <- c(rep(mean[ilist],nt))
                    wk.wts <- c(rep(wts[ilist],nt))
                    wk.vvariance <- c(rep(vvariance[ilist],nt))
                    wk.dmu <- c(rep(dmu[ilist],nt))
                             } else {
                    wk.x <- append(wk.x, as.vector(c(rep(x[ilist,],nt))), after=(nend*wks))
                    wk.mean <- append(wk.mean, c(rep(mean[ilist],nt)), after=nend)
                    wk.wts <- append(wk.wts, c(rep(wts[ilist],nt)), after=nend)
                    wk.vvariance <- append(wk.vvariance, c(rep(vvariance[ilist],nt)), after=nend)
                    wk.dmu <- append(wk.dmu, c(rep(dmu[ilist],nt)), after=nend) }
                nend <- nend + nt } } # end of for i loop 
                                 } # end of for ilist loop 
       wk.x <- t(matrix(wk.x, nrow=wks))
                                } # end of if model$data.type
    wkm <- sqrt(as.vector(wk.wts * wk.dmu * wk.dmu / wk.vvariance)) * wk.x
    xwx1 <- chol2inv(qr.R(qr(wkm)))
    as.vector(diag(wkm %*% xwx1 %*% t(wkm))) }

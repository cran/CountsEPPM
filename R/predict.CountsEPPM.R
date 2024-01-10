predict.CountsEPPM <-
function (object, newdata = NULL, type = c("response", "linear.predictor.mean",
                      "linear.predictor.scale.factor", "scale.factor", "mean", "variance", 
                      "distribution", "distribution.parameters"), na.action = na.pass, ...) {

    type <- match.arg(type)

    if (missing(newdata)) {
        nobs <- nrow(object$covariates.matrix.mean) 
        ntrials <- list(rep(c(0),nobs))
        ntrials <- lapply(1:nobs, function(i) {
           nmax <- object$vnmax[i]
           ntrials[[i]] <- c(rep(nmax,nmax+1)) } )
        output.model  <- Model.Counts(parameter=object$optim$par,
                           model.type=object$model.type,model.name=object$model.name,link=object$link,
                           covariates.matrix.mean=object$covariates.matrix.mean,
                           covariates.matrix.scalef=object$covariates.matrix.scalef,
                           offset.mean=object$offset.mean,offset.scalef=object$offset.scalef,
                           fixed.b=object$fixed.b,vnmax=object$vnmax) 
     } else {

        mf <- model.frame(delete.response(object$terms$mean), 
                          data = newdata, na.action = na.action, xlev = object$levels[["mean"]])
        newdata <- newdata[rownames(mf), , drop = FALSE]
        offset <- list(p = rep.int(0, nrow(mf)), scale.factor = rep.int(0, nrow(mf)))
        X <- model.matrix(delete.response(object$terms$mean), 
                   mf, contrasts = object$contrasts$mean)
        if (!is.null(object$call$offset)) 
            offset[[1L]] <- offset[[1L]] + eval(object$call$offset, newdata)
        if (!is.null(off.num <- attr(object$terms$mean, "offset"))) {
            for (j in off.num) offset[[1L]] <- offset[[1L]] + 
              eval(attr(object$terms$mean, "variables")[[j + 1L]], newdata) }
        Z <- matrix(c(rep(1,nrow(newdata))),ncol=1)
        if (object$model.type=="mean and scale-factor") {
            mf <- model.frame(delete.response(object$terms$scale.factor), 
                      data = newdata, na.action = na.action, xlev = object$levels[["scale.factor"]])
            Z <- model.matrix(delete.response(object$terms$scale.factor),
                       mf, contrasts = object$contrasts[["scale.factor"]])
            if (!is.null(off.num <- attr(object$terms$scale.factor, "offset"))) {
                for (j in off.num) offset[[2L]] <- offset[[2L]] + 
                  eval(attr(object$terms$scale.factor, "variables")[[j + 1L]], newdata) }
                                                     } # end of if p and scale.factor
          nobs <- nrow(newdata) 

          output.model  <- Model.Counts(parameter=object$optim$par,
                             model.type=object$model.type,model.name=object$model.name,link=object$link,
                             covariates.matrix.mean=X,covariates.matrix.scalef=Z,
                             offset.mean=offset[[1L]],offset.scalef=offset[[2L]],
                             fixed.b=object$fixed.b,vnmax=newdata$vnmax) } # end of if missing(newdata)
    vone          <- c(rep(1,nobs))
    mean.prob     <- c(rep(0, nobs))
    variance.prob <- c(rep(0, nobs))
    scalef.prob   <- vone
        rval <- switch(type, response = {
           if (missing(newdata)) { mean.prob <- object$fitted.values
                                 } else {
             va <- output.model$FDparameters$out.va
             vb <- output.model$FDparameters$out.vb
             vc <- output.model$FDparameters$out.vc
             v1mc <- vone - output.model$FDparameters$out.vc
             mean.prob <- rep(0, nobs)
             vone <- rep(1,nobs)
             va <- output.model$FDparameters$out.va
             vb <- output.model$FDparameters$out.vb
             vc <- output.model$FDparameters$out.vc
             v1mc <- vone - output.model$FDparameters$out.vc
             if ((object$model.name=="Poisson") | (object$model.name=="negative binomial") | 
                 (object$model.name=="negative binomial fixed b")) {
                if (object$model.name=="Poisson") { mean.prob <- va 
                                    } else { 
                    mean.prob <- vb*(attr(object$link,which="mean")$linkinv(va) - vone) }
                                                     } else { 
                mean.prob <- sapply(1:nobs, function(j) 
                   if ((abs(v1mc[j])<1.e-6)==TRUE) { 
                      mean.prob[j] <- vb[j]*(attr(object$link,which="mean")$linkinv(va[j]) - 1) 
                                         } else {  
                      mean.prob[j] <- vb[j]*((1+va[j]*v1mc[j]/(vb[j]**v1mc[j]))**(1/v1mc[j])-1) } ) }
                                            } # end of if missing(newdata)
             mean.prob
        }, linear.predictor.mean = {
           if (missing(newdata)) {
              vlp <- as.vector(object$covariates.matrix.mean %*% 
                      object$coefficients$mean.est[1:ncol(object$covariates.matrix.mean)] +
                      object$offset.mean)   
                              } else {
              vlp <- drop(X %*% object$coefficients$mean.est[1:ncol(X)] + offset[[1L]])   
                                     } # end of if missing(newdata)    
             vlp
        }, linear.predictor.scale.factor = {
             if (object$model.type=="mean and scale-factor") { 
                if (missing(newdata)) {
                   vlp  <- as.vector(object$covariates.matrix.scalef %*% 
                          object$coefficients$scalef.est[1:ncol(object$covariates.matrix.scalef)] +
                          object$offset.scalef)
                                   } else {
                   vlp  <- drop(Z %*% object$coefficients$scalef.est[1:ncol(Z)] + offset[[2L]])
                                          } # end of if missing(newdata)
                                      } else { vlp <- NULL}
             vlp
        }, scale.factor = { 
             if (object$model.name=="limiting") {
                valpha <- output.model$FDparameters$out.valpha
                vbeta  <- output.model$FDparameters$out.vbeta
                                   } else {
                va <- output.model$FDparameters$out.va
                vb <- output.model$FDparameters$out.vb
                vc <- output.model$FDparameters$out.vc
                v1mc <- vone - output.model$FDparameters$out.vc } # end if "limiting"
             if ((object$model.name=="Poisson") | (object$model.name=="negative binomial") | 
                 (object$model.name=="negative binomial fixed b")) {
                if (object$model.name!="Poisson") { 
                    scalef.prob <- (attr(object$link,which="mean")$linkinv(va) - vone) + vone } # end model.name!=Poisson
                                                     } else { 
                    if (object$model.name=="limiting") {
                        scalef.prob <- sapply(1:nobs, function(j) {
                           wks <- - 2*log(1 - valpha[j]*vbeta[j]) 
                           scalef.prob[j] <- (exp(wks) - 1) / wks } )
                                   } else {
                        mean.prob <- sapply(1:nobs, function(j) 
                          if ((abs(v1mc[j])<1.e-6)==TRUE) { 
                             mean.prob[j] <- vb[j]*(attr(object$link,which="mean")$linkinv(va[j]) - 1) 
                                                   } else {  
                             mean.prob[j] <- vb[j]*((1+va[j]*v1mc[j]/(vb[j]**v1mc[j]))**(1/v1mc[j])-1) } )
                        variance.prob <- sapply(1:nobs, function(j) 
                          if ((abs(v1mc[j])<1.e-6)==TRUE) { 
                             variance.prob[j] <- mean.prob[j]*(mean.prob[j]/vb[j] + 1)
                                                   } else {  
                             tcm1 <- 2*vc[j] - 1
                             mdb <- mean.prob[j] / vb[j] + 1
                             variance.prob[j] <- - vb[j]*mdb*(1 - mdb**tcm1) / tcm1 } )
                        scalef.prob <- variance.prob / mean.prob  
                                                             } } # end if "Poisson" | "negative binomial"
             scalef.prob
        }, mean = { 
             mean.prob <- rep(0, nobs)
             vone <- rep(1,nobs)
             if (object$model.name=="limiting") {
                valpha <- output.model$FDparameters$out.valpha
                vbeta  <- output.model$FDparameters$out.vbeta
                                   } else {
                va <- output.model$FDparameters$out.va
                vb <- output.model$FDparameters$out.vb
                vc <- output.model$FDparameters$out.vc
                v1mc <- vone - output.model$FDparameters$out.vc } # end if "limiting"
             if ((object$model.name=="Poisson") | (object$model.name=="negative binomial") | 
                 (object$model.name=="negative binomial fixed b")) {
                if (object$model.name=="Poisson") { mean.prob <- va 
                                    } else { 
                    mean.prob <- vb*(attr(object$link,which="mean")$linkinv(va) - vone) }
                                                     } else { 
                    if (object$model.name=="limiting") {
                        mean.prob <- sapply(1:nobs, function(j) 
                           mean.prob[j] <- - log(1 - valpha[j]*vbeta[j]) / vbeta[j] )
                                   } else {
                        mean.prob <- sapply(1:nobs, function(j) 
                           if ((abs(v1mc[j])<1.e-6)==TRUE) { 
                              mean.prob[j] <- vb[j]*(attr(object$link,which="mean")$linkinv(va[j]) - 1) 
                                                    } else {  
                              mean.prob[j] <- vb[j]*((1+va[j]*v1mc[j]/(vb[j]**v1mc[j]))**(1/v1mc[j])-1) } )
                                                             } } # end if "Poisson" | "negative binomial"
             mean.prob
        }, variance = { 
             mean.prob <- rep(0, nobs)
             variance.prob <- rep(0, nobs)
             vone <- rep(1,nobs)
             if (object$model.name=="limiting") {
                valpha <- output.model$FDparameters$out.valpha
                vbeta  <- output.model$FDparameters$out.vbeta
                                   } else {
                va <- output.model$FDparameters$out.va
                vb <- output.model$FDparameters$out.vb
                vc <- output.model$FDparameters$out.vc
                v1mc <- vone - output.model$FDparameters$out.vc } # end if "limiting"
             if ((object$model.name=="Poisson") | (object$model.name=="negative binomial") | 
                 (object$model.name=="negative binomial fixed b")) {
                if (object$model.name=="Poisson") { variance.prob <- va  
                                    } else { 
                    mean.prob <- vb*(attr(object$link,which="mean")$linkinv(va) - vone) 
                    variance.prob <- mean.prob*(mean.prob/vb + vone) }
                                                     } else { 
                    if (object$model.name=="limiting") {
                        mean.prob <- sapply(1:nobs, function(j) 
                           mean.prob[j] <- - log(1 - valpha[j]*vbeta[j]) / vbeta[j] )
                        variance.prob <- sapply(1:nobs, function(j) 
                           variance.prob[j] <- (exp(2*vbeta[j]*mean.prob[j]) - 1) / (2*vbeta[j]) )
                                   } else {
                        mean.prob <- sapply(1:nobs, function(j) 
                          if ((abs(v1mc[j])<1.e-6)==TRUE) { 
                             mean.prob[j] <- vb[j]*(attr(object$link,which="mean")$linkinv(va[j]) - 1) 
                                                   } else {  
                             mean.prob[j] <- vb[j]*((1+va[j]*v1mc[j]/(vb[j]**v1mc[j]))**(1/v1mc[j])-1) } )
                        variance.prob <- sapply(1:nobs, function(j) 
                          if ((abs(v1mc[j])<1.e-6)==TRUE) { 
                             variance.prob[j] <- mean.prob[j]*(mean.prob[j]/vb[j] + 1)
                                                   } else {  
                             tcm1 <- 2*vc[j] - 1
                             mdb <- mean.prob[j] / vb[j] + 1
                             variance.prob[j] <- - vb[j]*mdb*(1 - mdb**tcm1) / tcm1 } )
                                                             } } # end if "Poisson" | "negative binomial"
             variance.prob  
        }, distribution = { probabilities <- output.model$probabilities
# Truncation
              for ( i in 1:nobs) { probability <- output.model$probabilities[[i]]     
                 nmax1 <- object$vnmax[i] + 1
                 if ((is.na(object$ltvalue)==FALSE) | (is.na(object$utvalue)==FALSE)) { 
                    rev.probability <- LRTruncation(probability,object$ltvalue,object$utvalue)
                    if (is.na(object$ltvalue)==FALSE) { wks1 <- object$ltvalue + 2 
                                        } else { wks1 <- 1 }
                    if (is.na(object$utvalue)==FALSE) { wks2 <- object$utvalue  
                                               } else { wks2 <- nmax1 }
                 probabilities[[i]] <- rep(0,nmax1)
                 probabilities[[i]][wks1:wks2] <- rev.probability[wks1:wks2] }
                              } # end of for loop
             probabilities  
        }, distribution.parameters = { 


for ( i in 1:nobs) {
   nmax <- object$vnmax[i]
   nmax1 <- nmax + 1
   vnum <- c(0:nmax)
if (object$model.name!="limiting") {
   a   <- rep(output.model$FDparameters$out.va[i],nmax1)
   b   <- rep(output.model$FDparameters$out.vb[i],nmax1) 
   c   <- rep(output.model$FDparameters$out.vc[i],nmax1)
   if ((output.model$FDparameters$out.vc[i]!=1) & 
       (output.model$FDparameters$out.vc[i]!=0)) { 
      vlambda <- exp(log(a)+c*log(b+vnum))
                                             } else {
      if (output.model$FDparameters$out.vc[i]==0) { vlambda <- a }
      if (output.model$FDparameters$out.vc[i]==1) { vlambda <- a*(b+vnum) }
                } # end if c < 1 & = 0                                               
                                   } else {
      alpha <- rep(output.model$FDparameters$out.valpha[i],nmax1)
      beta  <- rep(output.model$FDparameters$out.beta[i],nmax1)
      vlambda <- alpha*exp(beta*vnum)
                                          } # end of model.name==limiting
# limiting value for lambda
      lambda.limit <- 745
      vlambda <- sapply(1:nmax1, function(j) 
         if ((is.finite(vlambda[j])==FALSE) | (vlambda[j]>lambda.limit)) { 
            vlambda[j] <- lambda.limit  
               } else { vlambda[j] <- vlambda[j] } )
                   } # end of for loop i

             output.model$FDparameters } )
        return(rval) }

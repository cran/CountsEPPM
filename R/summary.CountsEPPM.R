summary.CountsEPPM <-
function(object, ...) {

       nobs <- nrow(object$covariates.matrix.mean) 
       mean.par     <- rep(0,nobs)
       variance.par <- rep(0,nobs)
       scalef.par   <- rep(1,nobs)
       vone         <- rep(1,nobs)
       variance.limit <- rep(0,nobs)
# Calculation of means from parameter estimates and design matrices
       lp.mean <- object$covariates.matrix.mean%*%object$coefficients$mean.est + 
                    object$offset.mean
# inverse of link function
       mean.par <- attr(object$link, which="mean")$linkinv(lp.mean)
       wk.object <- object
       wk.object$coefficients <- as.vector(object$coefficients[1])
       wk.object$vcov <- vcov(object,model="mean")
       coeff.table.mean <- coeftest(wk.object, vcov.=vcov(object,model="mean"))

       if (is.null(object$coefficients$scalef.est)==FALSE) {
          wk.object <- object
          wk.object$coefficients <- as.vector(object$coefficients[2])
          wk.object$vcov <- vcov(object,model="scale.factor")
          coeff.table.scalef <- coeftest(wk.object, vcov.=vcov(object,model="scale.factor"))
                                                      } else { 
          coeff.table.scalef <- NULL } # end of is.null

# gradient 

      if (object$method=="BFGS") { 
          wk.gradient <- LL.gradient(parameter=object$optim$par,object$model.type,object$model.name,object$link,
                                     object$list.data,object$covariates.matrix.mean,object$covariates.matrix.scalef,
                                     object$offset.mean,object$offset.scalef,object$ltvalue,object$utvalue,
                                     object$fixed.b,object$weights,grad.method=attr(object$method,which="grad.method")) 
                          } else { wk.gradient <- NULL } # end of if object$method=="BFGS"

       object <- list(data.type=object$data.type, call=object$call, formula=object$formula, 
                    model.type=object$model.type, model.name=object$model.name, 
                    link=object$link, offset.mean=object$offset.mean,offset.scalef=object$offset.scalef,
                    coeff.table.mean=coeff.table.mean, coeff.table.scalef=coeff.table.scalef, loglik=object$loglik, 
                    n=object$nobs, nobs=object$nobs, df.null=object$df.null, df.residual=object$df.residual,
                    vnmax=object$vnmax, weights=object$weights, converged=object$converged, 
                    method=object$method,optim=object$optim,control=object$control,
                    fitted.values=object$fitted.values, y=object$y, terms=object$terms, 
                    ltvalue=object$ltvalue,utvalue=object$utvalue,
                    fixed.b=object$fixed.b, gradient=wk.gradient) 
      class(object) <- "summaryCountsEPPM"
      object }

logLik.CountsEPPM <-
function (object, ...) {
      structure(object$loglik, df = length(object$optim$par),
                nobs = nrow(object$covariates.matrix.mean),
                class = "logLik") }

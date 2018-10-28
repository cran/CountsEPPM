vcov.CountsEPPM <-
function (object, model = c("full", "mean", "scale.factor"), ...) {
    vc <- object$vcov
    k <- length(object$coefficients$mean.est)
    m <- length(object$coefficients$scalef.est)
    match.arg(model)
    if (missing(model)) { model <- c("full") } 
# Checking for correct model option
      if ((model!="full") & (model!="mean") & (model!="scale.factor")) {
         cat("\n","unknown model option","\n")
         vc <- NULL
                                   } else {
         switch(model, full = { vc 
                  }, mean = {
           vc <- vc[seq.int(from = 1, to = k, by = 1), 
              seq.int(from = 1, to = k, by = 1), drop = FALSE]
                  }, "scale.factor" = {
           if (m==0) { cat("\n","the scale-factor model has no elements","\n")
              vc <- NULL
                     } else { 
              vc <- vc[seq.int(length.out = m) + k, seq.int(length.out = m) + 
                  k, drop = FALSE]
              colnames(vc) <- rownames(vc) <- names(object$coefficients$scalef.est)
              vc }} ) } # end of if model!="full", etc.,
    return(vc) }

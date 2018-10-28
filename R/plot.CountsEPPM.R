plot.CountsEPPM <-
function (x, which = 1:4, caption = c("Residuals vs indices of obs.", 
    "Cook's distance plot", "Leverage vs predicted values", 
    "Residuals vs linear predictor", "Normal Q-Q plot of residuals", 
    "Predicted vs observed values"), sub.caption = " ", main = "", 
    ask = prod(par("mfcol"), 1) < length(which) && dev.interactive(), 
    ..., type = "spearson") 
  {
    if (!is.numeric(which) || any(which < 1) || any(which > 6)) 
        stop("'which' must be in 1:6")
    types <- c("pearson", "deviance", "response", "likelihood", 
        "sdeviance", "spearson")
    Types <- c("Pearson residuals", "Deviance residuals", "Raw response residuals", 
        "Likelihood residuals", "Standardized deviance residuals", 
        "Standardized Pearson residuals")
    type <- match.arg(type, types)
    Type <- Types[type == types]
    res <- residuals(x, type = type)
    n <- length(res)
    k <- length(x$coefficients$mean)
    show <- rep(FALSE, 6)
    show[which] <- TRUE
    Main <- rep("", 6)
    Main[which] <- rep(main, length.out = sum(show))
    one.fig <- prod(par("mfcol")) == 1
    if (ask) { op <- par(ask = TRUE)
               on.exit(par(op)) }
    if (show[1]) {
        plot(1:n, res, xlab = "Obs. number", ylab = Type, main = Main[1], 
            ...)
        if (one.fig) { title(sub = sub.caption, ...) }
        mtext(caption[1], 3, 0.25)
        abline(h = 0, lty = 3, col = "gray")
    }
    if (show[2]) {
        plot(1:n, cooks.distance(x), xlab = "Obs. number", ylab = "Cook's distance", 
            type = "h", main = Main[2])
        if (one.fig) { title(sub = sub.caption, ...) }
        mtext(caption[2], 3, 0.25)
    }
    if (show[3]) {
       if (x$data.type==TRUE) {
          fitted.values <- fitted.CountsEPPM(x)
                            } else {
          fitted.values <- c(rep(0,n))
          nstart <- 1
          nend <- 0
          for (ilist in 1:length(fitted.CountsEPPM(x))) { 
             ninlist <- sum(x$list.data[[ilist]])
             for ( i in 1:length(x$list.data[[ilist]])) { 
                nt <- x$list.data[[ilist]][i] 
                if (nt>0) {
                   nend <- nend + nt
                   fitted.values[nstart:nend] <- fitted.CountsEPPM(x)[ilist]
                   nstart <- nstart + nt } } # end of for i loop
                                   } # end of for ilist loop 
                                      } # end of if x$data.type  
        plot(fitted.values, hatvalues.CountsEPPM(x), xlab = "Predicted values", 
            ylab = "hatvalues as leverage", main = Main[3], ...)
        if (one.fig) { title(sub = sub.caption, ...) }
        mtext(caption[3], 3, 0.25)
    }
    if (show[4]) {
       if (x$data.type==TRUE) {
          linear.predictor <- predict(x, type = "linear.predictor.mean")
                            } else {
          linear.predictor <- c(rep(0,n))
          nstart <- 1
          nend <- 0
          for (ilist in 1:length(predict(x, type = "linear.predictor.mean"))) { 
             ninlist <- sum(x$list.data[[ilist]])
             for ( i in 1:length(x$list.data[[ilist]])) { 
                nt <- x$list.data[[ilist]][i] 
                if (nt>0) {
                   nend <- nend + nt
                   linear.predictor[nstart:nend] <- predict(x, type = "linear.predictor.mean")[ilist]
                   nstart <- nstart + nt } } # end of for i loop
                                   } # end of for ilist loop 
                                      } # end of if x$data.type  
        plot(linear.predictor, res, xlab = "Linear predictor mean", 
            ylab = Type, main = Main[4], ...)
        if (one.fig) { title(sub = sub.caption, ...) }
        mtext(caption[4], 3, 0.25)
        abline(h = 0, lty = 3, col = "gray")
    }
    if (show[5]) {
        qqnorm(y = as.vector(residuals.CountsEPPM(x, type)),
               main = Main[5], xlab = "Normal quantiles", ylab = Type)
        qqline(y = as.vector(residuals.CountsEPPM(x, type)),
               distribution = qnorm)
        if (one.fig) { title(sub = sub.caption, ...) }
        mtext(caption[5], 3, 0.25)
    }
    if (show[6]) {
        y <- if (is.null(x$y)) 
            model.response(model.frame(x))
        else x$y
        plot(y, fitted(x), xlab = "Observed values", ylab = "Predicted values", 
            main = Main[6], ...)
        if (one.fig) { title(sub = sub.caption, ...) }
        mtext(caption[6], 3, 0.25)
        abline(0, 1, lty = 2, col = "gray")
    }
    if (!one.fig && par("oma")[3] >= 1) 
        mtext(sub.caption, outer = TRUE, cex = 1.25)
    invisible() }

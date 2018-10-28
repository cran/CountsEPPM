residuals.CountsEPPM <-
function (object, type = c("spearson", "deviance", "pearson", 
    "response", "likelihood", "sdeviance"), ...) {
    type <- match.arg(type)
    mean <- object$fitted.values
    scale.factor <- predict(object, type = "scale.factor")
    nobs <- length(mean) 
    vone <- rep(1, nobs) 
    vnmax <- as.vector(object$vnmax) 
    vnmax1 <- vnmax + vone
    wts <- weights(object$model)
    if (object$data.type==TRUE) { 
       y <- if (is.null(object$y)) { model.response(model.frame(object)) 
                                   } else { object$y }
                            } else { y <- object$list.data } # end of data.type=TRUE
    if (is.null(wts)) { wts <- vone }

    if (object$data.type==TRUE) { 
       wk.y     <- y
       wk.mean     <- mean
       wk.scale.factor <- scale.factor
       wk.vnmax <- vnmax
       wk.wts   <- wts
       wk.vone  <- vone
       wk.resid <- wk.y - wk.mean
       total.ninlist <- nobs

                            } else { 

       vninlist <- c(rep(0,length(object$list.data)))
       vninlist <- sapply(1:length(object$list.data), function(ilist) 
                        vninlist[ilist] <- sum(object$list.data[[ilist]]) )
       total.ninlist <- sum(vninlist)

       wk.y     <- c(rep(0,total.ninlist))
       wk.mean     <- c(rep(0,total.ninlist))
       wk.scale.factor <- c(rep(1,total.ninlist))
       wk.vnmax <- c(rep(0,total.ninlist))
       wk.wts   <- c(rep(1,total.ninlist))
       wk.vone  <- c(rep(1,total.ninlist))
       wk.resid <- c(rep(0,total.ninlist))
       nstart <- 1
       nend <- 0
       for (ilist in 1:length(object$list.data)) { 
          ninlist <- sum(object$list.data[[ilist]])
          for ( i in 1:length(object$list.data[[ilist]])) { 
             nt <- object$list.data[[ilist]][i] 
             if (nt>0) {
                nend <- nend + nt
                wk.mean[nstart:nend]     <- mean[ilist]
                wk.scale.factor[nstart:nend] <- scale.factor[ilist]
                wk.vnmax[nstart:nend] <- vnmax[ilist]
                wk.y[nstart:nend]     <- i - 1 
                wk.wts[nstart:nend]   <- as.vector(wts[ilist])
                wk.resid[nstart:nend] <- wk.y[nstart:nend] - wk.mean[nstart:nend]
                nstart <- nstart + nt } } # end of for i loop
                                } # end of for ilist loop 
                                   } # end of if object$data.type  

    wk.resp.resid <- wk.resid
    if ((type=="pearson") | (type=="spearson") | (type=="likelihood")) {
       wk.resid <- wk.resid * sqrt( (wk.wts) / ( wk.mean * wk.scale.factor )) }
    if ((type=="deviance") | (type=="sdeviance") | (type=="likelihood")) {

        distribution <- predict(object, type = "distribution")

        ll.obs <- c(rep(0,total.ninlist))
        ll.obs <- sapply(1:total.ninlist, function(i) { 
           nmax1 <- object$vnmax[i] + 1
           vnum <- c(0:wk.vnmax[i])
           probability <- dpois(x=vnum,lambda=wk.y[i])
# Truncation
           if ((is.na(object$ltvalue)==FALSE) | (is.na(object$utvalue)==FALSE)) { 
              rev.probability <- LRTruncation(probability,object$ltvalue,object$utvalue)
              if (is.na(object$ltvalue)==FALSE) { wks1 <- object$ltvalue + 2 
                                  } else { wks1 <- 1 }
              if (is.na(object$utvalue)==FALSE) { wks2 <- object$utvalue  
                                  } else { wks2 <- nmax1 }
           probability[wks1:wks2] <- rev.probability[wks1:wks2] }
           wks <- probability[(wk.y[i]+1)]
           if (wks>0) { ll.obs[i] <- log(wks)
               } else { ll.obs[i] <- 0 } } ) # end of sapply

        ll.fit <- c(rep(0,total.ninlist))

        if (object$data.type==TRUE) { 

           ll.fit <- sapply(1:total.ninlist, function(i) { 
           probability <- distribution[[i]]
# Truncation
           if ((is.na(object$ltvalue)==FALSE) | (is.na(object$utvalue)==FALSE)) { 
              nmax1 <- object$vnmax[i] + 1
              rev.probability <- LRTruncation(probability,object$ltvalue,object$utvalue)
              if (is.na(object$ltvalue)==FALSE) { wks1 <- object$ltvalue + 2 
                                  } else { wks1 <- 1 }
              if (is.na(object$utvalue)==FALSE) { wks2 <- object$utvalue  
                                  } else { wks2 <- nmax1 }
           probability[wks1:wks2] <- rev.probability[wks1:wks2] }
           wks <- probability[(wk.y[i]+1)]
           if (wks>0) { ll.fit[i] <- log(wks)
               } else { ll.fit[i] <- 0 } } ) # end of sapply 

                             } else {

           wk.ind   <- c(rep(0,total.ninlist))
           wk.dist   <- c(rep(0,total.ninlist))
           nstart <- 1
           nend <- 0
           for (ilist in 1:length(object$list.data)) { 
              ninlist <- sum(object$list.data[[ilist]])
              for ( i in 1:length(object$list.data[[ilist]])) { 
                 nt <- object$list.data[[ilist]][i] 
                 if (nt>0) {
                    nend <- nend + nt
                    wk.ind[nstart:nend] <- i
                    wk.dist[nstart:nend] <- ilist
                    nstart <- nstart + nt } } # end of for i loop
                                     } # end of for ilist loop 

        ll.fit <- sapply(1:total.ninlist, function(i) { 
           probability <- distribution[[wk.dist[i]]]
# Truncation
           if ((is.na(object$ltvalue)==FALSE) | (is.na(object$utvalue)==FALSE)) { 
              nmax1 <- object$vnmax[i] + 1
              rev.probability <- LRTruncation(probability,object$ltvalue,object$utvalue)
              if (is.na(object$ltvalue)==FALSE) { wks1 <- object$ltvalue + 2 
                                  } else { wks1 <- 1 }
              if (is.na(object$utvalue)==FALSE) { wks2 <- object$utvalue  
                                  } else { wks2 <- nmax1 }
           probability[wks1:wks2] <- rev.probability[wks1:wks2] }
           ll.fit[i] <- log(probability[wk.ind[i]]) } ) # end of sapply ll.fit
                                } # end of if object$data.type
        wk.resid.dev <- sqrt(wk.wts)*sign(wk.resid)*sqrt(2*abs(ll.obs - ll.fit)) 
                                   } # end of if deviance or sdeviance or likelihood

    res <- switch(type, pearson = {
        wk.resid
    }, response = {
        wk.resp.resid
    }, deviance = {
        wk.resid.dev
    }, likelihood = {
# likelihood residuals
          sqrt(wk.wts)*sign(wk.resid.dev)*sqrt( hatvalues(object)*wk.resid*wk.resid  
                    + wk.resid.dev*wk.resid.dev ) 
    }, sdeviance = {
          wkv <- rep(0,length(wk.resid.dev)) 
# deviance residuals standardized to have an asymptotic variance of 1
          wkv <- sapply(1:length(wk.resid.dev), function(i) { 
              if (wk.scale.factor[i]<=0) { wkv[i] <- wkv[i]
                              } else { wkv[i] <-
                  wk.resid.dev[i] / sqrt( wk.scale.factor[i] * (1 - hatvalues(object)[i]) )
                                     } } ) 
          wkv
    }, spearson = {
# Pearson residuals standardized to have an asymptotic variance of 1
          wk.resid / sqrt( wk.vone - hatvalues(object) )
    }) # end of switch
    return(res) }

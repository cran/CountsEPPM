LRTruncation <-
function (probability,ltvalue,utvalue) {
    nmax1 <- length(probability)
    wks   <- sum((probability==1.e-8))
    rev.probability <- rep(1.e-8,nmax1)
# Checking for an improper distribution initial distribution
    if (wks<nmax1) { x <- 0 
       if (is.na(utvalue)==FALSE) { x <- 1 - sum(probability) 
               if (nmax1>utvalue) { x <- x + sum(probability[utvalue:nmax1]) } }
       if (is.na(ltvalue)==FALSE) { wks <- ltvalue + 1
                                      x <- x + sum(probability[1:wks]) }
       if ((x>0) & (x<1)) { rev.probability <- probability / (1 - x)
                          } else { if (x<=0) { rev.probability <- probability }
                                 } } # end of if (wks<nmax1) 
    if (is.na(utvalue)==FALSE) { wks <- utvalue + 1
            if (nmax1>utvalue) { rev.probability[wks:nmax1] <- NA } }
    if (is.na(ltvalue)==FALSE) { wks <- ltvalue + 1
                                 rev.probability[1:wks] <- NA } 
return(rev.probability) }

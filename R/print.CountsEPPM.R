print.CountsEPPM <-
function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 
        0.85)), "", sep = "\n")
    if ((x$converged)==FALSE) { cat("model did not converge\n") 
                              } else { 
        if (length(x$coefficients$mean.est)) {
            cat(paste("Coefficients (model for mean with ", x$link, 
                " link):\n", sep = ""))
            print.default(format(x$coefficients$mean.est, digits = digits), 
                print.gap = 2, quote = FALSE)
            cat("\n")
        }
        else cat("No coefficients (in mean model)\n\n")
            if (length(x$coefficients$scalef.est)) {
                cat(paste("Coefficients (model for scale-factor with log link):\n", sep = ""))
                print.default(format(x$coefficients$scalef.est, 
                  digits = digits), print.gap = 2, quote = FALSE)
                cat("\n")
            }
            else cat("No coefficients (in precision model)\n\n")
        } 
    invisible(x) }

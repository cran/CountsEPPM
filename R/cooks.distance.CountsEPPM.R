cooks.distance.CountsEPPM <-
function (model, ...)  {
    h <- hatvalues.CountsEPPM(model)
    k <- length(model$coefficients$mean.est)
    res <- residuals.CountsEPPM(model, type = "pearson")
    h * (res^2)/(k * (1 - h)^2) }

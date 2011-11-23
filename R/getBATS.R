# TODO: Add comment
# 
# Author: Slava
###############################################################################


residuals.bats<-function(object, ...) {
	residuals<-msts(object$errors, seasonal.periods=object$seasonal.periods, start=object$start.time, ...)
	return(residuals)
}

fitted.bats<-function(object, ...) {
	residuals<-msts(object$fitted.values, seasonal.periods=object$seasonal.periods, start=object$start.time, ...)
	return(residuals)
}

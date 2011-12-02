# TODO: Add comment
# 
# Author: srazbash
###############################################################################
#setwd("/Volumes/NO\ NAME/SlavaBATS/")
#setwd("E:/SlavaBATS")
#source("fitBATS.R", verbose=TRUE)
#source("forecastBATS.R", verbose=TRUE)
#source("makeMatrices.R", verbose=TRUE)
#source("checkAdmissibility.R", verbose=TRUE)
#source("makeParamVector.R", verbose=TRUE)
#source("adjustSeasonalSeeds.R", verbose=TRUE)
#source("getBATS.R", verbose=TRUE)

filterSpecifics<-function(y, box.cox, trend, damping, seasonal.periods, use.arma.errors, force.seasonality=FALSE, ...) {
	if((trend == FALSE) & (damping == TRUE)) {
		return(list(AIC=Inf))
	}
	#printCASE(box.cox, trend, damping, seasonal.periods, NULL, NULL, 0, 0)
	first.model<-fitSpecificBATS(y, use.box.cox=box.cox, use.beta=trend, use.damping=damping, seasonal.periods=seasonal.periods)
	if((!is.null(seasonal.periods)) & (!force.seasonailty)) {
		#printCASE(box.cox, trend, damping, NULL, NULL, NULL, 0, 0)
		non.seasonal.model<-fitSpecificBATS(y, use.box.cox=box.cox, use.beta=trend, use.damping=damping, seasonal.periods=NULL)
		if(first.model$AIC > non.seasonal.model$AIC) {
			seasonal.periods<-NULL
			first.model<-non.seasonal.model
		}
	}
	if(use.arma.errors) { 
		##Turn off warnings
		old.warning.level <- options()$warn
		options(warn=-1)
		arma<-auto.arima(as.numeric(first.model$e), d=0, ...)
		###Re-enable warnings
		options(warn=old.warning.level)
		p<-arma$arma[1]
		q<-arma$arma[2]
		if((p != 0) | (q != 0)) { #Did auto.arima() find any AR() or MA() coefficients?
			if(p != 0) {
				ar.coefs<-numeric(p)
			} else {
				ar.coefs<-NULL
			}
			if(q != 0) {
				ma.coefs<-numeric(q)
			} else {
				ma.coefs<-NULL
			}
			starting.params<-first.model$parameters
			#printCASE(box.cox, trend, damping, seasonal.periods, ar.coefs, ma.coefs, p, q)
			second.model<-fitSpecificBATS(y, use.box.cox=box.cox, use.beta=trend, use.damping=damping, seasonal.periods=seasonal.periods, ar.coefs=ar.coefs, ma.coefs=ma.coefs)
			if(second.model$AIC < first.model$AIC) {
				return(second.model)
			} else {
				return(first.model)
			}
		} else { #Else auto.arima() did not find any AR() or MA()coefficients
			return(first.model)
		}
	} else {
		return(first.model)
	}
}


bats<-function(y, use.box.cox=NULL, use.trend=NULL, use.damped.trend=NULL, seasonal.periods=NULL, use.arma.errors=TRUE, ...) {
	if(any((y <= 0))) {
		stop("BATS requires positive data")
	}
	if(any(class(y) == "msts")) {
		start.time<-start(y)
		seasonal.periods<-attr(y,"msts")
		y<-as.numeric(y)
	} else if(class(y) == "ts") {
		start.time<-start(y)
		seasonal.periods<-frequency(y)
		y<-as.numeric(y)
	}  else {
		start.time<-1
		y<-as.numeric(y)
	}
	best.aic<-NULL
	if(is.null(use.box.cox)) {
		use.box.cox<-c(FALSE, TRUE)
	} 
	if(is.null(use.trend)) {
		use.trend<-c(FALSE, TRUE)
	} else if(use.trend == FALSE) {
		use.damped.trend<-FALSE
	}
	if(is.null(use.damped.trend)) {
		use.damped.trend<-c(FALSE, TRUE)
	}
	for(box.cox in use.box.cox) {
		for(trend in use.trend) {
			for(damping in use.damped.trend) {
				current.model<-filterSpecifics(y, box.cox=box.cox, trend=trend, damping=damping, seasonal.periods=seasonal.periods, use.arma.errors=use.arma.errors, ...)
				if(!is.null(best.aic)) {
					if(current.model$AIC < best.aic) {
						best.aic<-current.model$AIC
						best.model<-current.model
					}
				} else {
					best.model<-current.model
					best.aic<-best.model$AIC
				}
			}
		}
	}
	best.model$call<-match.call()
	best.model$start.time<-start.time
	return(best.model)
}

print.bats<-function(object) {
	cat("\n")
	cat("BATS( {")
	if(!is.null(object$lambda)) {
		cat(object$lambda)
	} else {
		cat("1")
	}
	cat("}, {")
	if(!is.null(object$ar.coefficients)) {
		cat(length(object$ar.coefficients))
	} else {
		cat("0")
	}
	cat(", ")
	if(!is.null(object$ma.coefficients)) {
		cat(length(object$ma.coefficients))
	} else {
		cat("0")
	}
	cat("}, {")
	if(!is.null(object$damping.parameter)) {
		cat(object$damping.parameter)
	} else {
		cat("0")
	}
	
	if(!is.null(object$seasonal.periods)) {
		cat("}, { ")
		for(i in object$seasonal.periods) {
			cat(i)
			if(i != object$seasonal.periods[length(object$seasonal.periods)]) {
				cat(", ")
			} else {
				cat("})")
			}
		}
	} else {
		cat("})\n\n")	
	}
	cat("\nCall: ")
	print(object$call)
	cat("\nParameters:\n")
	cat("\nBox-Cox Parameter: ")
	cat(object$lambda)
	cat("\nAlpha: ")
	cat(object$alpha)
	cat("\nBeta: ")
	cat(object$beta)
	cat("\nDamping Parameter: ")
	cat(object$damping.parameter)
	cat("\nGamma Values: ")
	cat(object$gamma.values)
	cat("\nAR() Coefficients: ")
	cat(object$ar.coefficients)
	cat("\nMA() Coefficients: ")
	cat(object$ma.coefficients)
	cat("\n\n")
	cat("\nSeed States:\n")
	print(object$seed.states)
	
	cat("\nSigma: ")
	cat(sqrt(object$variance))
	
	cat("\n\nAIC: ")
	cat(object$AIC)
	cat("\n\n")
	
}

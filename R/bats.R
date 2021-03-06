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
	

	first.model <- fitSpecificBATS(y, use.box.cox=box.cox, use.beta=trend, use.damping=damping, seasonal.periods=seasonal.periods)
	if((!is.null(seasonal.periods)) & (!force.seasonality)) {
		non.seasonal.model <- fitSpecificBATS(y, use.box.cox=box.cox, use.beta=trend, use.damping=damping, seasonal.periods=NULL)
		if(first.model$AIC > non.seasonal.model$AIC) {
			seasonal.periods <- NULL
			first.model <- non.seasonal.model
		}
	}
	if(use.arma.errors) { 
		##Turn off warnings
		old.warning.level <- options()$warn
		options(warn=-1)
		arma <- auto.arima(as.numeric(first.model$errors), d=0, ...)
		###Re-enable warnings
		options(warn=old.warning.level)
		p <- arma$arma[1]
		q <- arma$arma[2]
		if((p != 0) | (q != 0)) { #Did auto.arima() find any AR() or MA() coefficients?
			if(p != 0) {
				ar.coefs <- numeric(p)
			} else {
				ar.coefs <- NULL
			}
			if(q != 0) {
				ma.coefs <- numeric(q)
			} else {
				ma.coefs <- NULL
			}
			starting.params <- first.model$parameters
			#printCASE(box.cox, trend, damping, seasonal.periods, ar.coefs, ma.coefs, p, q)
			second.model <- fitSpecificBATS(y, use.box.cox=box.cox, use.beta=trend, use.damping=damping, seasonal.periods=seasonal.periods, ar.coefs=ar.coefs, ma.coefs=ma.coefs)
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

parFilterSpecifics<-function(control.number, control.array, y, seasonal.periods, use.arma.errors, force.seasonality=FALSE, ...) {
	box.cox <- control.array[control.number, 1] 
	trend <- control.array[control.number, 2]
	damping <- control.array[control.number, 3]
	
	
	if((trend == FALSE) & (damping == TRUE)) {
		return(list(AIC=Inf))
	}
	
	
	first.model <- fitSpecificBATS(y, use.box.cox=box.cox, use.beta=trend, use.damping=damping, seasonal.periods=seasonal.periods)
	if((!is.null(seasonal.periods)) & (!force.seasonality)) {
		non.seasonal.model <- fitSpecificBATS(y, use.box.cox=box.cox, use.beta=trend, use.damping=damping, seasonal.periods=NULL)
		if(first.model$AIC > non.seasonal.model$AIC) {
			seasonal.periods <- NULL
			first.model <- non.seasonal.model
		}
	}
	if(use.arma.errors) { 
		##Turn off warnings
		old.warning.level <- options()$warn
		options(warn=-1)
		arma <- auto.arima(as.numeric(first.model$errors), d=0, ...)
		###Re-enable warnings
		options(warn=old.warning.level)
		p <- arma$arma[1]
		q <- arma$arma[2]
		if((p != 0) | (q != 0)) { #Did auto.arima() find any AR() or MA() coefficients?
			if(p != 0) {
				ar.coefs <- numeric(p)
			} else {
				ar.coefs <- NULL
			}
			if(q != 0) {
				ma.coefs <- numeric(q)
			} else {
				ma.coefs <- NULL
			}
			starting.params <- first.model$parameters
			#printCASE(box.cox, trend, damping, seasonal.periods, ar.coefs, ma.coefs, p, q)
			second.model <- fitSpecificBATS(y, use.box.cox=box.cox, use.beta=trend, use.damping=damping, seasonal.periods=seasonal.periods, ar.coefs=ar.coefs, ma.coefs=ma.coefs)
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

bats <- function(y, use.box.cox=NULL, use.trend=NULL, use.damped.trend=NULL, seasonal.periods=NULL, use.arma.errors=TRUE, use.parallel=TRUE, num.cores=NULL, ...) {
	if(any((y <= 0))) {
		stop("BATS requires positive data")
	}
  origy <- y
	if(any(class(y) == "msts")) {
		seasonal.periods <- attr(y,"msts")
	} else if(class(y) == "ts") {
		seasonal.periods <- frequency(y)
	}
	if(all((seasonal.periods == 1))) {
		seasonal.periods <- NULL	
	}
	if(!is.null(seasonal.periods)) {
		seasonal.mask <- (seasonal.periods == 1)
		seasonal.periods <- seasonal.periods[!seasonal.mask]
	}
  y <- as.numeric(y)
	best.aic <- NULL
	if(is.null(use.box.cox)) {
		use.box.cox <- c(FALSE, TRUE)
	} 
	if(is.null(use.trend)) {
		use.trend <- c(FALSE, TRUE)
	} else if(use.trend == FALSE) {
		use.damped.trend <- FALSE
	}
	if(is.null(use.damped.trend)) {
		use.damped.trend <- c(FALSE, TRUE)
	}
	if(use.parallel) {
		#Set up the control array
		control.array <- NULL
		for(box.cox in use.box.cox) {
			for(trend in use.trend) {
				for(damping in use.damped.trend) {
					if((trend == FALSE) & (damping == TRUE)) {
						next
					}
					control.line <- c(box.cox, trend, damping)
					if(!is.null(control.array)) {
						control.array <- rbind(control.array, control.line)
					} else {
						control.array <- control.line
					}
				}
			}
		}
		##Fit the models
		if(is.null(num.cores)) {
			num.cores<-detectCores(all.tests = FALSE, logical = TRUE)
		}
		clus <- makeCluster(num.cores)
		models.list <- clusterApplyLB(clus, c(1:nrow(control.array)), parFilterSpecifics, y=y, control.array=control.array, seasonal.periods=seasonal.periods, use.arma.errors=use.arma.errors)
		stopCluster(clus)
		##Choose the best model
		####Get the AICs
		aics <- numeric(nrow(control.array))
		for(i in 1:nrow(control.array)) {
			aics[i] <- models.list[[i]]$AIC
		}
		best.number <- which.min(aics)
		best.model <- models.list[[best.number]]

	} else {
		for(box.cox in use.box.cox) {
			for(trend in use.trend) {
				for(damping in use.damped.trend) {
					current.model <- filterSpecifics(y, box.cox=box.cox, trend=trend, damping=damping, seasonal.periods=seasonal.periods, use.arma.errors=use.arma.errors, ...)
					if(!is.null(best.aic)) {
						if(current.model$AIC < best.aic) {
							best.aic <- current.model$AIC
							best.model <- current.model
						}
					} else {
						best.model <- current.model
						best.aic <- best.model$AIC
					}
				}
			}
		}
	}
	best.model$call <- match.call()
	if(best.model$optim.return.code != 0) {
		warning("optim() did not converge.")
	}
  
  # Add ts attributes
  if(!any(class(origy) == "ts"))
  {
    if(is.null(seasonal.periods))
      origy <- ts(origy,s=1,f=1)
    else 
      origy <- msts(origy,seasonal.periods)
  }
  attributes(best.model$fitted.values) <- attributes(best.model$errors) <- attributes(origy)
  best.model$y <- origy
  
	return(best.model)
}

print.bats <- function(x,...) {
	cat(makeText(x))
	cat("\n")
	cat("\nCall: ")
	print(x$call)
	cat("\nParameters")
  if(!is.null(x$lambda))
  {
    cat("\n  Lambda: ")
    cat(round(x$lambda,6))
  }
	cat("\n  Alpha: ")
	cat(x$alpha)
  if(!is.null(x$beta))
  {
    cat("\n  Beta: ")
    cat(x$beta)
    cat("\n  Damping Parameter: ")
    cat(round(x$damping.parameter,6))
	}
  if(!is.null(x$gamma.values))
   {
    cat("\n  Gamma Values: ")
    cat(x$gamma.values)
  }
  if(!is.null(x$ar.coefficients))
  {
    cat("\n  AR coefficients: ")
    cat(round(x$ar.coefficients,6))
	}
  if(!is.null(x$ma.coefficients))
  {
    cat("\n  MA coefficients: ")
    cat(round(x$ma.coefficients,6))
  }
	cat("\n")
	cat("\nSeed States:\n")
	print(x$seed.states)
	
	cat("\nSigma: ")
	cat(sqrt(x$variance))
	
	cat("\nAIC: ")
	cat(x$AIC)
	cat("\n")	
}
residuals.bats <- function(object, ...) {
	object$errors
}


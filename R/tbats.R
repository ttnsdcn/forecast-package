tbats<-function(y, use.box.cox=NULL, use.trend=NULL, use.damped.trend=NULL, seasonal.periods=NULL, use.arma.errors=TRUE, ...) {
	if(any((y <= 0))) {
		stop("TBATS requires positive data")
	}
	non.seasonal.model<-bats(as.numeric(y), use.box.cox=use.box.cox, use.trend=use.trend, use.damped.trend=use.damped.trend, use.arma.errors=use.arma.errors)
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
		if(is.null(seasonal.periods)) {
			return(best.model)
		}
	}
	if(is.null(seasonal.periods)) {
		return(not.seasonal.model)
	}
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
	#Set a vector of model params for later comparison
	model.params<-logical(length=3)
	model.params[1]<-any(use.box.cox)
	model.params[2]<-any(use.trend)
	model.params[3]<-any(use.damped.trend)
	
	###The OLS setup
	#Get seasonal states
	bats.states<-bats(y, model.params[1], model.params[2], model.params[3], seasonal.periods=seasonal.periods, force.seasonality=TRUE)$x
	#
	if(model.params[2]) {
		adj.beta<-1
	} else {
		adj.beta<-2
	}
	seasonals<-numeric(length(y)*length(seasonal.periods))
	dim(seasonals)<-c(length(y), length(seasonal.periods))
	previous.season<-0
	k.vector<-rep(1, length(seasonal.periods))
	n<-length(y)
	for(i in 1:length(seasonal.periods)) {
		seasonals[,i]<-as.numeric(bats.states[(1+adj.beta+previous.season+seasonal.periods[i])])
		p.val<-0
		fourier.terms<-makeSingleFourier(1, seasonal.periods[i], n)
		previous.sse<-sum(residuals(lm(seasonals[,i] ~ fourier.terms -1))^2)
		#while((p.val < .001) & ((2*k.vector[i]) < (seasonal.periods[i]-1))) {
		#	
		#}
		repeat {
			if((2*(k.vector[i]+1)) >= (seasonal.periods[i]-1)) {
				break
			}
			new.fourier.terms<-makeSingleFourier((k.vector[i]+1), seasonal.periods[i], n)
			new.sse<-sum(residuals(lm(seasonals[,i] ~ fourier.terms + new.fourier.terms -1))^2)
			p.val<-calcFTest(previous.sse, new.sse, 2, (2 + ncol(fourier.terms)), n)
			if(p.val > .001) {
				break
			} else {
				k.vector[i]<-k.vector[i]+1
				four.terms<-cbind(four.terms, new.four.terms)
			}
		}
	}
	best.model<-fitSpecificTBATS(y, model.params[1], model.params[2], model.params[3], seasonal.periods, k.vector)
	for(i in 1:length(seasonal.periods)) {
		max.k<-floor(((seasonal.periods[i]-1)/2))
		#repeat {
			if(k.vector[i] == max.k) {
				next
			}
			if(max.k <= 6) {
				old.k<-k.vector[i]
				k.vector[i]<-max.k
				repeat {
					new.model<-fitSpecificTBATS(y, model.params[1], model.params[2], model.params[3], seasonal.periods, k.vector)
					if(new.model$AIC >= best.model$AIC) {
						k.vector[i]<-old.k
						break
					} else {
						old.k<-k.vector[i]
						k.vector[i]<-k.vector[i]-1
						best.model<-new.model
					}
				}
				print("here-db")
				next
			}
			if(k.vector[i] >= 6) {
				k.vector[i]<-k.vector[i]+1
				new.model<-fitSpecificTBATS(y, model.params[1], model.params[2], model.params[3], seasonal.periods, k.vector)
				if(new.model$AIC >= best.model$AIC) {
					k.vector[i]<-k.vector[i]-1
					break
				} else {
					best.model<-new.model
				}
			} else if(max.k > 6) {
				#step.up.k<-k.vector
				#step.down.k<-k.vector
				#step.up.k[i]<-7
				#step.down.k[i]<-5
				
				#up.model<-fitSpecificTBATS(y, model.params[1], model.params[2], model.params[3], seasonal.periods, step.up.k)
				#down.model<-fitSpecificTBATS(y, model.params[1], model.params[2], model.params[3], seasonal.periods, step.down.k)
				
				#if(up.model$AIC < down.model$AIC) {
					
				#} else {
				#	if(down.model$AIC < best.model$AIC)
				#}
				
				
			}
			
		#}
	}
	aux.model<-best.model
	if(non.seasonal.model$AIC < best.model$AIC) {
		#print(best.model)
		#return(non.seasonal.model)
		best.model<-non.seasonal.model
	}
	#print(best.model)
	for(box.cox in use.box.cox) {
		for(trend in use.trend) {
			for(damping in use.damped.trend) {
				if(all((model.params == c(box.cox, trend, damping)))) {
					new.model<-filterTBATSSpecifics(y, box.cox, trend, damping, seasonal.periods, k.vector, use.arma.errors, aux.model=aux.model, ...)
					#print(new.model)
					#print("####-a")
					#print(new.model$AIC)
					#print("$$$$-a")
					if(new.model$AIC < best.model$AIC) {
						best.model<-new.model	
					}
				} else if(!((trend == FALSE) & (damping == TRUE))) {
					new.model<-filterTBATSSpecifics(y, box.cox, trend, damping, seasonal.periods, k.vector, use.arma.errors, ...)
					#print(new.model)
					#print("####")
					#print(new.model$AIC)
					#print("$$$$")
					if(new.model$AIC < best.model$AIC) {
						best.model<-new.model	
					}
				}
				
				
			}
		}
	}
	
	
	best.model$call<-match.call()
	best.model$start.time<-start.time
	return(best.model)
}

filterTBATSSpecifics<-function(y, box.cox, trend, damping, seasonal.periods, k.vector, use.arma.errors, aux.model=NULL, ...) {
	if(is.null(aux.model)) {
		first.model<-fitSpecificTBATS(y, use.box.cox=box.cox, use.beta=trend, use.damping=damping, seasonal.periods=seasonal.periods, k.vector=k.vector)
	} else {
		first.model<-aux.model	
	}
	if(use.arma.errors) { 
		##Turn off warnings
		old.warning.level <- options()$warn
		options(warn=-1)
		arma<-try(auto.arima(as.numeric(first.model$e), d=0, ...), silent=TRUE)
		###Re-enable warnings
		options(warn=old.warning.level)
		if(class(arma) != "try-error") {
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
	
				second.model<-fitSpecificTBATS(y, use.box.cox=box.cox, use.beta=trend, use.damping=damping, seasonal.periods=seasonal.periods, k.vector=k.vector, ar.coefs=ar.coefs, ma.coefs=ma.coefs)
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
	} else {
		return(first.model)
	}
}


makeSingleFourier<-function(j, m, T) {
	frier<-matrix(0, nrow=T, ncol=2)
	for(t in 1:T) {
		frier[t,1]<-cos((2*pi*j)/m)
		frier[t,2]<-sin((2*pi*j)/m)
	}
	return(frier)
} 

calcFTest<-function(r.sse, ur.sse, num.restrictions, num.u.params, num.observations) {
	f.stat<-((r.sse - ur.sse)/num.restrictions)/(r.sse/(num.observations - num.u.params))
	p.value<-pf(f.stat, num.restrictions, (num.observations - num.u.params),lower.tail=FALSE )
	return(p.value)
}


print.tbats<-function(object) {
	cat("\n")
	cat("TBATS( {")
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
		k.pos<-1
		for(i in object$seasonal.periods) {
			cat("< ")
			cat(object$k.vector[k.pos])
			cat(", ")
			cat(i)
			cat(" >")
			k.pos<-k.pos+1
			if(i != object$seasonal.periods[length(object$seasonal.periods)]) {
				cat(", ")
			} else {
				cat("})")
			}
		}
	} else {
		cat("{)\n\n")	
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
	cat("\nGamma-1 Values: ")
	cat(object$gamma.one.values)
	cat("\nGamma-2 Values: ")
	cat(object$gamma.two.values)
	cat("\nK-vector: ")
	cat(object$k.vector)
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
# TODO: Add comment
# 
# Author: srazbash
###############################################################################


forecast.tbats<-function(object, h=10, level=c(80,95), fan=FALSE, ts.frequency=(if(!is.null(object$seasonal.periods)) { max(object$seasonal.periods)} else { 1})) {
	if(h<=0) {
		stop("Forecast horizon out of bounds")
	}
	if(fan) {
		level <- seq(51,99,by=3)
	}
	
	if(!is.null(object$k.vector)) {
		tau<-2*object$k.vector
	} else {
		tau<-0
	}
	x<-matrix(0,nrow=nrow(object$x), ncol=h)
	y.forecast<-numeric(h)
	if(!is.null(object$beta)) {
		adj.beta<-1
	} else {
		adj.beta<-0
	}
	
	w<-.Call("makeTBATSWMatrix", smallPhi_s = object$damping.parameter, kVector_s=as.integer(object$k.vector), arCoefs_s = object$ar.coefficients, maCoefs_s = object$ma.coefficients, tau_s=as.integer(tau), PACKAGE = "forecast")
	
	if(!is.null(object$seasonal.periods)) {
		gamma.bold<-matrix(0,nrow=1,ncol=tau)
		.Call("updateTBATSGammaBold", gammaBold_s=gamma.bold, kVector_s=as.integer(object$k.vector), gammaOne_s=object$gamma.one.v, gammaTwo_s=object$gamma.two.v, PACKAGE = "forecast")
	} else {
		gamma.bold<-NULL	
	}
	g<-matrix(0, nrow=(tau+1+adj.beta+object$p+object$q), ncol=1)
	if(object$p != 0) {
		g[(1+adj.beta+tau+1),1]<-1
	}
	if(object$q != 0) {
		g[(1+adj.beta+tau+object$p+1),1]<-1
	}
	.Call("updateTBATSGMatrix", g_s=g, gammaBold_s=gamma.bold, alpha_s=object$alpha, beta_s=object$beta.v, PACKAGE = "forecast")
	
	#print(g)
	
	F<-makeTBATSFMatrix(alpha=object$alpha, beta=object$beta, small.phi=object$damping.parameter, seasonal.periods=object$seasonal.periods, k.vector=as.integer(object$k.vector), gamma.bold.matrix=gamma.bold, ar.coefs=object$ar.coefficients, ma.coefs=object$ma.coefficients)
	
	y.forecast[1]<-w$w.transpose %*% object$x[,ncol(object$x)]
	x[,1]<- F %*% object$x[,ncol(object$x)] + g %*% object$e[length(object$e)]
	
	for(t in 2:h) {
		x[,t]<-F %*% x[,(t-1)]
		y.forecast[t]<-w$w.transpose %*% x[,(t-1)]
	}
	
	##Make prediction intervals here
	lower.bounds <- upper.bounds <- matrix(NA,ncol=length(level),nrow=h)
	variance.multiplier<-numeric(h)
	variance.multiplier[1]<-1
	if(h > 1) {
		for(j in 1:(h-1)) {
			if(j == 1) {
				f.running<-diag(ncol(F))
			} else {
				f.running<-f.running %*% F
			}				
			c.j<-w$w.transpose %*% f.running %*% g 
			variance.multiplier[(j+1)]<-variance.multiplier[j]+ c.j^2
		}
	}
	
	variance<-object$variance * variance.multiplier
	#print(variance)
	st.dev<-sqrt(variance)
	for(i in 1:length(level)) {
		marg.error <- st.dev * abs(qnorm((100-level[i])/200))
		lower.bounds[,i] <- y.forecast - marg.error
		upper.bounds[,i] <- y.forecast + marg.error
		
	}
	#Inv Box Cox transform if required
	if(!is.null(object$lambda))
	{
		y.forecast <- InvBoxCox(y.forecast,object$lambda)
		lower.bounds <- InvBoxCox(lower.bounds,object$lambda)
		upper.bounds <- InvBoxCox(upper.bounds,object$lambda)
	}
	##Calc a start time for the forecast
	y<-object$y
	y[(length(y)+1)]<-0
	y<-ts(y, start=object$start.time, frequency=ts.frequency)
	fcast.start.time<-end(y)
	#Make msts object for x and mean
	x<-msts(object$y, seasonal.periods=(if(!is.null(object$seasonal.periods)) { object$seasonal.periods} else { 1}), ts.frequency=ts.frequency, start=object$start.time)
	fitted.values<-msts(object$fitted.values, seasonal.periods=(if(!is.null(object$seasonal.periods)) { object$seasonal.periods} else { 1}), start=object$start.time)
	y.forecast<-msts(y.forecast, seasonal.periods=(if(!is.null(object$seasonal.periods)) { object$seasonal.periods} else { 1}), start=fcast.start.time)
		
	forecast.object<-list(model=object, mean=y.forecast, level=level, x=x, upper=upper.bounds, lower=lower.bounds, fitted=fitted.values, method=makeTextTBATS(object), residuals=object$e)
	class(forecast.object)<-"forecast"
	return(forecast.object)
}


makeTextTBATS<-function(object) {
	name<-"TBATS( {"
	if(!is.null(object$lambda)) {
		name<-paste(name, round(object$lambda, digits=6), sep="")
	} else {
		name<-paste(name, "1", sep="")
	}
	name<-paste(name, "}, {", sep="")
	if(!is.null(object$ar.coefficients)) {
		name<-paste(name, length(object$ar.coefficients), sep="")
	} else {
		name<-paste(name, "0", sep="")
	}
	name<-paste(name, ", ", sep="")
	if(!is.null(object$ma.coefficients)) {
		name<-paste(name, length(object$ma.coefficients), sep="")
	} else {
		name<-paste(name, "0", sep="")
	}
	name<-paste(name, "}, {", sep="")
	if(!is.null(object$damping.parameter)) {
		name<-paste(name, round(object$damping.parameter, digits=6), sep="")
	} else {
		name<-paste(name, "0", sep="")
	}
	
	if(!is.null(object$seasonal.periods)) {
		name<-paste(name, "}, { ", sep="")
		k.pos<-1
		for(i in object$seasonal.periods) {
			name<-paste(name, "< ", sep="")
			name<-paste(name, object$k.vector[k.pos], sep="")
			name<-paste(name, ", ", sep="")
			name<-paste(name, i, sep="")
			name<-paste(name, " >", sep="")
			k.pos<-k.pos+1
			if(i != object$seasonal.periods[length(object$seasonal.periods)]) {
				name<-paste(name, ", ", sep="")
			} else {
				name<-paste(name, "})", sep="")
			}
		}
	} else {
		name<-paste(name, "})", sep="")	
	}
	return(name)
}


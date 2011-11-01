# TODO: Add comment

msts<-function(data, seasonal.periods, ts.frequency=floor(max(seasonal.periods)), ...) {
	object<-ts(data=data, frequency=ts.frequency, ...)
	class(object)<-c("msts", "ts")
	attr(object, "msts")<-seasonal.periods
	return(object)
}

print.msts<-function(object) {
	if(length(attr(object, "msts")) == 1) {
		attr(object, "msts")<-NULL
		print.ts(object)
	} else {
		cat("Multi-Seasonal Time Series:\n")
		cat("Start: ")
		cat(start(object))
		#cat("\nEnd: ")
		#cat(object$end)
		cat("\nSeasonal Periods: ")
		cat(attr(object,"msts"))
		cat("\nData:\n")
		print(as.numeric(object))
		#print(matrix(object, ncol=length(object)), nrow=1)
		cat("\n")
	}
}

window.msts<-function(object, ...) {
	seasonal.periods<-attr(object,"msts")
	class(object)<-c("ts")
	object<-window(object, ...)
	class(object)<-c("msts", "ts")
	attr(object, "msts")<-seasonal.periods
	return(object)
}


### Functions to handle seasonality

monthdays <- function(x)
{
    if(!is.ts(x))
        stop("Not a time series")
    f <- frequency(x)
    if(f==12)
        days <- c(31,28,31,30,31,30,31,31,30,31,30,31)
    else if(f==4)
        days <- c(90,91,92,92)
    else
        stop("Not monthly or quarterly data")
    nyears <- round(length(x)/f+1)+1
    years <- (1:nyears) + (start(x)[1] - 1)
    leap.years <- ((years %% 4 == 0) & !(years %% 100 ==0 & years %% 400 != 0))[1:nyears]
    dummy <- t(matrix(rep(days,nyears),nrow=f))
    if(f==12)
        dummy[leap.years,2] <- 29
    else
        dummy[leap.years,1] <- 91
    xx <- c(t(dummy))[start(x)[2]-1+(1:length(x))]
    return(ts(xx,start=start(x),f=f))
}

sindexf <- function(object,h)
{
    if(class(object)=="stl")
    {
        ss <- object$time.series[,1]
        m <- frequency(ss)
        ss <- ss[length(ss)-(m:1)+1]
        tsp.x <- tsp(object$time.series)
    }
    else if(class(object)=="decomposed.ts")
    {
        ss <- object$figure
        m <- frequency(object$seasonal)
        n <- length(object$trend)
        ss <- rep(ss,n/m+1)[1:n]
        ss <- ss[n-(m:1)+1]
        tsp.x <- tsp(object$seasonal)
    }
    else
        stop("Object of unknown class")
    out <- ts(rep(ss,h/m+1)[1:h], f=m, start=tsp.x[2]+1/m)

    return(out)
}

seasadj <- function(object)
{
    if(class(object)=="stl")
        return(object$time.series[,2]+object$time.series[,3])
    else if(class(object)=="decomposed.ts")
    {
        if(object$type=="additive")
            return(object$x-object$seasonal)
        else
            return(object$x/object$seasonal)
    }
    else
        stop("Object of unknown class")
}

seasonaldummy <- function(x)
{
    if(!is.ts(x))
        stop("Not a time series")
    else
        fr.x <- frequency(x)
    if(fr.x==1)
        stop("Non-seasonal time series")
    dummy <- as.factor(cycle(x))
    dummy.mat <- matrix(0,ncol=frequency(x)-1,nrow=length(x))
    nrow <- 1:length(x)
    for(i in 1:(frequency(x)-1))
        dummy.mat[dummy==paste(i),i] =1
    colnames(dummy.mat) <- if (fr.x == 12)
                month.abb[1:11]
            else if(fr.x == 4)
                c("Q1", "Q2", "Q3")
            else paste("S",1:(fr.x-1),sep="")

    return(dummy.mat)
}

seasonaldummyf <- function(x, h)
{
    if(!is.ts(x))
        stop("Not a time series")
    f=frequency(x)
    return(seasonaldummy(ts(rep(0,h),start=tsp(x)[2]+1/f,freq=f)))
}

forecast.stl <- function(object, method=c("ets","arima"), etsmodel="ZZN",
     h = frequency(object$time.series)*2, level = c(80, 95), fan = FALSE, lambda=NULL, ...)
{
  method <- match.arg(method)
  m <- frequency(object$time.series)
  n <- nrow(object$time.series)
  lastseas <- rep(object$time.series[n-(m:1)+1,"seasonal"],trunc(1+(h-1)/m))[1:h]
  # De-seasonalize
  x.sa <- seasadj(object)
  # Forecast
  if(method=="ets")
  {
    # Ensure non-seasonal model
    if(substr(etsmodel,3,3) != "N")
    {
      warning("The ETS model must be non-seasonal. I'm ignoring the seasonal component specified.")
      substr(etsmodel,3,3) <- "N"
    }
    fit <- ets(x.sa,model=etsmodel,...)
  }
  else
    fit <- auto.arima(x.sa,D=0,max.P=0,max.Q=0,...)
  fcast <- forecast(fit,h=h,level=level,fan=fan)
  # Reseasonalize
  fcast$mean <- fcast$mean + lastseas
  fcast$upper <- fcast$upper + lastseas
  fcast$lower <- fcast$lower + lastseas
  fcast$x <- ts(rowSums(object$time.series))
  tsp(fcast$x) <- tsp(object$time.series)
  fcast$method <- paste("STL + ",fcast$method)
  fcast$seasonal <- ts(lastseas[1:m],f=m,start=tsp(object$time.series)[2]-1+1/m)
  fcast$fitted <- fitted(fcast)+object$time.series[,1]
  fcast$residuals <- fcast$x - fcast$fitted
  
	if (!is.null(lambda)) 
	{
		fcast$x <- InvBoxCox(fcast$x,lambda)
		fcast$fitted <- InvBoxCox(fcast$fitted, lambda)
		fcast$mean <- InvBoxCox(fcast$mean, lambda)
		fcast$lower <- InvBoxCox(fcast$lower, lambda)
		fcast$upper <- InvBoxCox(fcast$upper, lambda)
		fcast$lambda <- lambda
	}
  
   return(fcast)
}

stlf <- function(x ,h=frequency(x)*2, s.window=7, method=c("ets","arima"), etsmodel="ZZN", level = c(80, 95), fan = FALSE, lambda=NULL, ...)
{
	if (!is.null(lambda)) 
	{
		origx <- x
		x <- BoxCox(x, lambda)
	}

	fit <- stl(x,s.window=s.window)
	fcast <- forecast(fit,h=h,method=method,etsmodel=etsmodel, level=level,fan=fan,...)

	if (!is.null(lambda)) 
	{
		fcast$x <- origx
		fcast$fitted <- InvBoxCox(fcast$fitted, lambda)
		fcast$mean <- InvBoxCox(fcast$mean, lambda)
		fcast$lower <- InvBoxCox(fcast$lower, lambda)
		fcast$upper <- InvBoxCox(fcast$upper, lambda)
		fcast$lambda <- lambda
	}

	return(fcast)
}

fourier <- function(x, K)
{
    n <- length(x)
    period <- frequency(x)
    X <- matrix(,nrow=n,ncol=2*K)
    for(i in 1:K)
    {
        X[,2*i-1] <- sin(2*pi*i*(1:n)/period)
        X[,2*i] <- cos(2*pi*i*(1:n)/period)
    }
    colnames(X) <- paste(c("S","C"),rep(1:K,rep(2,K)),sep="")
    return(X)
}

fourierf <- function(x, K, h)
{
    n <- length(x)
    period <- frequency(x)
    X <- matrix(,nrow=h,ncol=2*K)
    for(i in 1:K)
    {
        X[,2*i-1] <- sin(2*pi*i*((n+1):(n+h))/period)
        X[,2*i] <- cos(2*pi*i*((n+1):(n+h))/period)
    }
    colnames(X) <- paste(c("S","C"),rep(1:K,rep(2,K)),sep="")
    return(X)
}

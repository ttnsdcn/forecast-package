ets <- function(y, model="ZZZ", damped=NULL,
    alpha=NULL, beta=NULL, gamma=NULL, phi=NULL, additive.only=FALSE, lambda=NULL, 
    lower=c(rep(0.0001,3), 0.8), upper=c(rep(0.9999,3),0.98),
    opt.crit=c("lik","amse","mse","sigma","mae"), nmse=3, bounds=c("both","usual","admissible"),
    ic=c("aic","aicc","bic"),restrict=TRUE)
{
    #dataname <- substitute(y)
    opt.crit <- match.arg(opt.crit)
    bounds <- match.arg(bounds)
    ic <- match.arg(ic)
      
    #if(max(y,na.rm=TRUE) > 1e6)
    #    warning("Very large numbers which may cause numerical problems. Try scaling the data first")

    if(class(y)=="data.frame" | class(y)=="list" | class(y)=="matrix" | is.element("mts",class(y)))
        stop("y should be a univariate time series")
    y <- as.ts(y)
    # Remove missing values near ends
    ny <- length(y)
    y <- na.contiguous(y)
    if(ny != length(y))
        warning("Missing values encountered. Using longest contiguous portion of time series")

  orig.y <- y
  if(!is.null(lambda))
  {
    y <- BoxCox(y,lambda)
    additive.only=TRUE
  }
    
    if(nmse < 1 | nmse > 10)
        stop("nmse out of range")
    m <- frequency(y)

    if(sum((upper-lower)>0)<4)
        stop("Lower limits must be less than upper limits")

    # Define model
    if(class(model)=="ets")
    {
        alpha=model$par["alpha"]
        beta=model$par["beta"]
        if(is.na(beta))
            beta <- NULL
        gamma=model$par["gamma"]
        if(is.na(gamma))
            gamma <- NULL
        phi=model$par["phi"]
        if(is.na(phi))
            phi <- NULL
        damped=(model$components[4]=="TRUE")
        model=paste(model$components[1],model$components[2],model$components[3],sep="")
    }

    errortype  <- substr(model,1,1)
    trendtype  <- substr(model,2,2)
    seasontype <- substr(model,3,3)

    if(!is.element(errortype,c("M","A","Z")))
        stop("Invalid error type")
    if(!is.element(trendtype,c("N","A","M","Z")))
        stop("Invalid trend type")
    if(!is.element(seasontype,c("N","A","M","Z")))
        stop("Invalid season type")

    if(m < 1)
    {
      warning("I can't handle data with frequency less than 1. Seasonality will be ignored.")
      seasontype=="N"
    }
    if(m == 1)
    {
      if(seasontype=="A" | seasontype=="M")
        stop("Nonseasonal data")
      else
        substr(model,3,3) <- seasontype <- "N"
    }
    if(m > 24)
    {
      if(is.element(seasontype,c("A","M")))
        stop("Frequency too high")
      else if(seasontype=="Z")
      {
        warning("I can't handle data with frequency greater than 24. Seasonality will be ignored. Try stlf() if you need seasonal forecasts.")
        substr(model,3,3) <- seasontype <- "N"
        #m <- 1
      }
    }

    # Check inputs
    if(restrict)
    {
        if((errortype=="A" & (trendtype=="M" | seasontype=="M")) |
            (errortype=="M" & trendtype=="M" & seasontype=="A") |
            (additive.only & (errortype=="M" | trendtype=="M" | seasontype=="M")))
        stop("Forbidden model combination")
    }
    
    data.positive <- (min(y) > 0)

    if(!data.positive & errortype=="M")
        stop("Inappropriate model for data with negative or zero values")

    if(!is.null(damped))
    {
        if(damped & trendtype=="N")
            stop("Forbidden model combination")
    }

    # Fit model (assuming only one nonseasonal model)
    if(errortype=="Z")
        errortype <- c("A","M")
    if(trendtype=="Z")
        trendtype <- c("N","A","M")
    if(seasontype=="Z")
        seasontype <- c("N","A","M")
    if(is.null(damped))
        damped <- c(TRUE,FALSE)
    best.ic <- Inf
    for(i in 1:length(errortype))
    {
        for(j in 1:length(trendtype))
        {
            for(k in 1:length(seasontype))
            {
                for(l in 1:length(damped))
                {
                    if(trendtype[j]=="N" & damped[l])
                        next
                    if(restrict)
                    {
                        if(errortype[i]=="A" & (trendtype[j]=="M" | seasontype[k]=="M"))
                            next
                        if(errortype[i]=="M" & trendtype[j]=="M" & seasontype[k]=="A")
                            next
                        if(additive.only & (errortype[i]=="M" | trendtype[j]=="M" | seasontype[k]=="M"))
                            next
                    }
                    if(!data.positive & errortype[i]=="M")
                        next
                    fit <- etsmodel(y,errortype[i],trendtype[j],seasontype[k],damped[l],alpha,beta,gamma,phi,
                        lower=lower,upper=upper,opt.crit=opt.crit,nmse=nmse,bounds=bounds)
                    fit.ic <- switch(ic,aic=fit$aic,bic=fit$bic,aicc=fit$aicc)
                    if(!is.na(fit.ic))
                    {
                        if(fit.ic < best.ic)
                        {
                            model <- fit
                            best.ic <- fit.ic
                            best.e <- errortype[i]
                            best.t <- trendtype[j]
                            best.s <- seasontype[k]
                            best.d <- damped[l]
                        }
                    }
                }
            }
        }
    }
    if(best.ic == Inf)
        stop("No model able to be fitted")

    model$m <- m
    model$method <- paste("ETS(",best.e,",",best.t,ifelse(best.d,"d",""),",",best.s,")",sep="")
    model$components <- c(best.e,best.t,best.s,best.d)
    model$call <- match.call()
    model$initstate <- model$states[1,]
    model$sigma2 <- mean(model$residuals^2,na.rm=TRUE)
    model$x <- orig.y
  model$lambda <- lambda
  if(!is.null(lambda))
  {
    model$fitted <- InvBoxCox(model$fitted,lambda)
  }
  
    #model$call$data <- dataname

    return(structure(model,class="ets"))
}

etsmodel <- function(y, errortype, trendtype, seasontype, damped,
  alpha=NULL, beta=NULL, gamma=NULL, phi=NULL,
  lower, upper, opt.crit, nmse, bounds)
{
    tsp.y <- tsp(y)
    if(is.null(tsp.y))
        tsp.y <- c(1,length(y),1)
    if(seasontype != "N")
      m <- tsp.y[3]
    else
      m <- 1

    # Initialize smoothing parameters
    par <- initparam(alpha,beta,gamma,phi,trendtype,seasontype,damped,lower,upper,m)
    names(alpha) <- names(beta) <- names(gamma) <- names(phi) <- NULL
    par.noopt <- c(alpha=alpha,beta=beta,gamma=gamma,phi=phi)
    if(!is.null(par.noopt))
        par.noopt <- c(na.omit(par.noopt))
    if(!is.na(par["alpha"]))
        alpha <- par["alpha"]
    if(!is.na(par["beta"]))
        beta <- par["beta"]
    if(!is.na(par["gamma"]))
        gamma <- par["gamma"]
    if(!is.na(par["phi"]))
        phi <- par["phi"]

#    if(errortype=="M" | trendtype=="M" | seasontype=="M")
#        bounds="usual"
    if(!check.param(alpha,beta,gamma,phi,lower,upper,bounds,m))
    {
        print(paste("Model: ETS(",errortype,",",trendtype,ifelse(damped,"d",""),",",seasontype,")",sep=""))
        stop("Parameters out of range")
    }

    # Initialize state
    init.state <- initstate(y,trendtype,seasontype)
    nstate <- length(init.state)
    par <- c(par,init.state)
    lower <- c(lower,rep(-Inf,nstate))
    upper <- c(upper,rep(Inf,nstate))

    np <- length(par)
    if(np >= length(y)-1) # Not enough data to continue
        return(list(aic=Inf,bic=Inf,aicc=Inf,mse=Inf,amse=Inf,fit=NULL,par=par,states=init.state))

    # Optimize parameters and state
    if(length(par)==1)
        tmp <- options(warn=-1)$warn # Turn off warning on 1-d Nelder-Mead.

    fred <- optim(par,lik,method="Nelder-Mead",y=y,nstate=nstate, errortype=errortype, trendtype=trendtype,
#    fred <- optim(par,lik,method="BFGS",y=y,nstate=nstate, errortype=errortype, trendtype=trendtype,
            seasontype=seasontype, damped=damped, par.noopt=par.noopt, lowerb=lower, upperb=upper,
        opt.crit=opt.crit, nmse=nmse, bounds=bounds, m=m,pnames=names(par),pnames2=names(par.noopt),
        control=list(maxit=2000))
    fit.par <- fred$par

    names(fit.par) <- names(par)

    if(length(par)==1)
        options(warn=tmp)

    init.state <- fit.par[(np-nstate+1):np]
    # Add extra state
    if(seasontype!="N")
         init.state <- c(init.state, m*(seasontype=="M") - sum(init.state[(2+(trendtype!="N")):nstate]))

    if(!is.na(fit.par["alpha"]))
        alpha <- fit.par["alpha"]
    if(!is.na(fit.par["beta"]))
        beta <- fit.par["beta"]
    if(!is.na(fit.par["gamma"]))
        gamma <- fit.par["gamma"]
    if(!is.na(fit.par["phi"]))
        phi <- fit.par["phi"]
    e <- pegelsresid.C(y,m,init.state,errortype,trendtype,seasontype,damped,alpha,beta,gamma,phi)

    n <- length(y)
    aic <- e$lik + 2*np
    bic <- e$lik + log(n)*np
    aicc <- e$lik + 2*n*np/(n-np-1)

    mse <- e$amse[1]
    amse <- mean(e$amse[1:nmse])

    states=ts(e$states,f=tsp.y[3],s=tsp.y[1]-1/tsp.y[3])
    colnames(states)[1] <- "l"
    if(trendtype!="N")
        colnames(states)[2] <- "b"
    if(seasontype!="N")
        colnames(states)[(2+(trendtype!="N")):ncol(states)] <- paste("s",1:m,sep="")

    tmp <- c("alpha",rep("beta",trendtype!="N"),rep("gamma",seasontype!="N"),rep("phi",damped))
    fit.par <- c(fit.par,par.noopt)
#    fit.par <- fit.par[order(names(fit.par))]
    if(errortype=="A")
        fits <- y-e$e
    else
        fits <- y/(1+e$e)

    return(list(loglik=-0.5*e$lik,aic=aic,bic=bic,aicc=aicc,mse=mse,amse=amse,fit=fred,residuals=ts(e$e,f=tsp.y[3],s=tsp.y[1]),fitted=ts(fits,f=tsp.y[3],s=tsp.y[1]),
        states=states,par=fit.par))
}

initparam <- function(alpha,beta,gamma,phi,trendtype,seasontype,damped,lower,upper,m)
{
    # Set up initial parameters
    par <- numeric(0)
    if(is.null(alpha))
    {
        if(m > 12)
            alpha <- 0.0002
        if(is.null(beta) & is.null(gamma))
            alpha <- lower[1] + .5*(upper[1]-lower[1])
        else if(is.null(gamma))
            alpha <- beta+0.001
        else if(is.null(beta))
            alpha <- 0.999-gamma
        else 
            alpha <- 0.5*(beta - gamma + 1)
        if(alpha < lower[1] | alpha > upper[1])
            stop("Inconsistent parameter limits")
        par <- alpha
        names(par) <- "alpha"
    }
    if(is.null(beta))
    {
        if(trendtype !="N")
        {
            if(m > 12)
                beta <- 0.00015
            else
                beta <- lower[2] + .1*(upper[2]-lower[2])
            if(beta > alpha)
                beta <- min(alpha - 0.0001,0.0001)
            if(beta < lower[2] | beta > upper[2])
                stop("Can't find consistent starting parameters")
            par <- c(par,beta)
            names(par)[length(par)] <- "beta"
        }
    }
    if(is.null(gamma))
    {
        if(seasontype !="N")
        {
            if(m > 12)
                gamma <- 0.0002
            else
                gamma <- lower[3] + .01*(upper[3]-lower[3])
            if(gamma > 1-alpha)
                gamma <- min(0.999-alpha,0.001)
            if(gamma < lower[3] | gamma > upper[3])
                stop("Can't find consistent starting parameters")
            par <- c(par,gamma)
            names(par)[length(par)] <- "gamma"
        }
    }
    if(is.null(phi))
    {
        if(damped)
        {
            phi <- lower[4] + .99*(upper[4]-lower[4])
            par <- c(par,phi)
            names(par)[length(par)] <- "phi"
        }
    }

    return(par)
}


check.param <- function(alpha,beta,gamma,phi,lower,upper,bounds,m)
{
    if(bounds != "admissible")
    {
        if(!is.null(alpha))
        {
            if(alpha < lower[1] | alpha > upper[1])
                return(0)
        }
        if(!is.null(beta))
        {
            if(beta < lower[2] | beta > alpha | beta > upper[2])
                return(0)
        }
        if(!is.null(phi))
        {
            if(phi < lower[4] | phi > upper[4])
                return(0)
        }
        if(!is.null(gamma))
        {
            if(gamma < lower[3] | gamma > 1-alpha | gamma > upper[3])
                return(0)
        }
    }
    if(bounds != "usual")
    {
        if(!admissible(alpha,beta,gamma,phi,m))
            return(0)
    }
    return(1)
}

initstate <- function(y,trendtype,seasontype)
{
    if(seasontype!="N")
    {
        # Do decomposition
        m <- frequency(y)
        if(length(y)>=3*m)
            y.d <- decompose(y,type=switch(seasontype,A="additive",M="multiplicative"))
        else
            y.d <- list(seasonal= switch(seasontype,A=y-mean(y),M=y/mean(y)))
        init.seas <- rev(y.d$seasonal[2:m])
        names(init.seas) <- paste("s",0:(m-2),sep="")
        if(seasontype=="A")
            y.sa <- y-y.d$seasonal
        else
            y.sa <- y/y.d$seasonal
    }
    else
    {
        m <- 1 
        init.seas <- NULL
        y.sa <- y
    }
    maxn <- min(max(10,2*m),length(y.sa))
    fit <- lsfit(1:maxn,y.sa[1:maxn])
    l0 <- fit$coef[1]
    if(trendtype=="A")
        b0 <- fit$coef[2]
    else if(trendtype=="M")
    {
        b0 <- 1+fit$coef[2]/fit$coef[1]
        if(abs(b0) > 1e10) # Avoid infinite slopes
            b0 <- sign(b0)*1e10
    }
    else
        b0 <- NULL
    names(l0) <- "l"
    if(!is.null(b0))
        names(b0) <- "b"
    return(c(l0,b0,init.seas))
}


lik <- function(par,y,nstate,errortype,trendtype,seasontype,damped,par.noopt,lowerb,upperb,
    opt.crit,nmse,bounds,m,pnames,pnames2)
{
    names(par) <- pnames
    names(par.noopt) <- pnames2
    alpha <- c(par["alpha"],par.noopt["alpha"])["alpha"]
    if(is.na(alpha))
        stop("alpha problem!")
    if(trendtype!="N")
    {
        beta <- c(par["beta"],par.noopt["beta"])["beta"]
        if(is.na(beta))
            stop("beta Problem!")
    }
    else
        beta <- NULL
    if(seasontype!="N")
    {
        gamma <- c(par["gamma"],par.noopt["gamma"])["gamma"]
        if(is.na(gamma))
            stop("gamma Problem!")
    }
    else
    {
      m <- 1
      gamma <- NULL
    }
    if(damped)
    {
        phi <- c(par["phi"],par.noopt["phi"])["phi"]
        if(is.na(phi))
            stop("phi Problem!")
    }
    else
        phi <- NULL

    if(!check.param(alpha,beta,gamma,phi,lowerb,upperb,bounds,m))
        return(1e12)

    np <- length(par)

    init.state <- par[(np-nstate+1):np]
    # Add extra state
    if(seasontype!="N")
        init.state <- c(init.state, m*(seasontype=="M") - sum(init.state[(2+(trendtype!="N")):nstate]))
    # Check states
    if(seasontype=="M")
    {
        seas.states <- init.state[-(1:(1+(trendtype!="N")))]
        if(min(seas.states) < 0)
            return(1e8)
    }

    e <- pegelsresid.C(y,m,init.state,errortype,trendtype,seasontype,damped,alpha,beta,gamma,phi)

    if(is.na(e$lik))
        return(1e8)
    if(e$lik < -1e10) # Avoid perfect fits
        return(-1e10)

#    points(alpha,e$lik,col=2)

    if(opt.crit=="lik")
        return(e$lik)
    else if(opt.crit=="mse")
        return(e$amse[1])
    else if(opt.crit=="amse")
        return(mean(e$amse[1:nmse]))
    else if(opt.crit=="sigma")
        return(mean(e$e^2))
    else if(opt.crit=="mae")
      return(mean(abs(e$e)))
}

print.ets <- function(x,...)
{
    cat(paste(x$method, "\n\n"))
    cat(paste("Call:\n", deparse(x$call), "\n\n"))
    ncoef <- length(x$initstate)
    if(!is.null(x$lambda))
      cat("  Box-Cox transformation: lambda=",round(x$lambda,4), "\n\n")

    cat("  Smoothing parameters:\n")
    cat(paste("    alpha =", round(x$par["alpha"], 4), "\n"))
    if(x$components[2]!="N")
        cat(paste("    beta  =", round(x$par["beta"], 4), "\n"))
    if(x$components[3]!="N")
        cat(paste("    gamma =", round(x$par["gamma"], 4), "\n"))
    if(x$components[4]!="FALSE")
        cat(paste("    phi   =", round(x$par["phi"], 4), "\n"))

    cat("\n  Initial states:\n")
    cat(paste("    l =", round(x$initstate[1], 4), "\n"))
    if (x$components[2]!="N")
        cat(paste("    b =", round(x$initstate[2], 4), "\n"))
    else
    {
        x$initstate <- c(x$initstate[1], NA, x$initstate[2:ncoef])
        ncoef <- ncoef+1
    }
    if (x$components[3]!="N")
    {
        cat("    s=")
        if (ncoef <= 8)
            cat(round(x$initstate[3:ncoef], 4))
        else
        {
            cat(round(x$initstate[3:8], 4))
            cat("\n           ")
            cat(round(x$initstate[9:ncoef], 4))
        }
        cat("\n")
    }

    cat("\n  sigma:  ")
    cat(round(sqrt(x$sigma2),4))
    stats <- c(x$aic,x$aicc,x$bic)
    names(stats) <- c("AIC","AICc","BIC")
    cat("\n\n")
    print(stats)
#    cat("\n  AIC:    ")
#    cat(round(x$aic,4))
#    cat("\n  AICc:   ")
#    cat(round(x$aicc,4))
#    cat("\n  BIC:    ")
#    cat(round(x$bic,4))
}


pegelsresid.C <- function(y,m,init.state,errortype,trendtype,seasontype,damped,alpha,beta,gamma,phi)
{
    n <- length(y)
    p <- length(init.state)
    x <- numeric(p*(n+1))
    x[1:p] <- init.state
    e <- numeric(n)
    lik <- 0;
    if(!damped)
        phi <- 1;
    if(trendtype == "N")
        beta <- 0;
    if(seasontype == "N")
        gamma <- 0;

    Cout <- .C("etscalc",
        as.double(y),
        as.integer(n),
        as.double(x),
        as.integer(m),
        as.integer(switch(errortype,"A"=1,"M"=2)),
        as.integer(switch(trendtype,"N"=0,"A"=1,"M"=2)),
        as.integer(switch(seasontype,"N"=0,"A"=1,"M"=2)),
        as.double(alpha),
        as.double(beta),
        as.double(gamma),
        as.double(phi),
        as.double(e),
        as.double(lik),
        as.double(numeric(10)),
        PACKAGE="forecast")
    if(!is.na(Cout[[13]]))
    {
        if(abs(Cout[[13]]+99999) < 1e-7)
            Cout[[13]] <- NA
    }
    tsp.y <- tsp(y)
    e <- ts(Cout[[12]])
    tsp(e) <- tsp.y 

    return(list(lik=Cout[[13]], amse=Cout[[14]], e=e, states=matrix(Cout[[3]], nrow=n+1, ncol=p, byrow=TRUE)))
}

admissible <- function(alpha,beta,gamma,phi,m)
{
    if(is.null(phi))
        phi <- 1
    if(phi < 0 | phi > 1+1e-8)
        return(0)
    if(is.null(gamma))
    {
        if(alpha < 1-1/phi | alpha > 1+1/phi)
            return(0)
        if(!is.null(beta))
        {
            if(beta < alpha * (phi-1) | beta > (1+phi)*(2-alpha))
                return(0)
        }
    }
    else # Seasonal model
    {
        if(is.null(beta))
            beta <- 0
        if(gamma < max(1-1/phi-alpha,0) | gamma > 1+1/phi-alpha)
            return(0)
        if(alpha < 1-1/phi-gamma*(1-m+phi+phi*m)/(2*phi*m))
            return(0)
        if(beta < -(1-phi)*(gamma/m+alpha))
            return(0)

        # End of easy tests. Now use characteristic equation
        P <- c(phi*(1-alpha-gamma),alpha+beta-alpha*phi+gamma-1,rep(alpha+beta-alpha*phi,m-2),(alpha+beta-phi),1)
        roots <- polyroot(P)
        if(max(abs(roots)) > 1+1e-10)
            return(0)
    }
    # Passed all tests
    return(1)
}

### PLOT COMPONENTS
plot.ets <- function(x,...)
{
  if(!is.null(x$lambda))
    y <- BoxCox(x$x,x$lambda)
  else
    y <- x$x
  if(x$components[3]=="N" & x$components[2]=="N")
  {
      plot(cbind(observed=y, level=x$states[,1]),
          main=paste("Decomposition by",x$method,"method"),...)
  }
  else if(x$components[3]=="N")
  {
      plot(cbind(observed=y, level=x$states[,1], slope=x$states[,"b"]),
          main=paste("Decomposition by",x$method,"method"),...)
  }
  else if(x$components[2]=="N")
  {
      plot(cbind(observed=y, level=x$states[,1], season=x$states[,"s1"]),
          main=paste("Decomposition by",x$method,"method"),...)
  }
  else
  {
      plot(cbind(observed=y, level=x$states[,1], slope=x$states[,"b"],
          season=x$states[,"s1"]),
          main=paste("Decomposition by",x$method,"method"),...)
  }
}

summary.ets <- function(object,...)
{
   print(object)
   cat("\nIn-sample error measures:\n")
   print(accuracy(object))
}

coef.ets <- function(object,...)
{
   object$par
}

logLik.ets <- function(object,...)
{
    structure(object$loglik,df=length(object$par),class="logLik")
}

\name{forecast.lm}
\alias{forecast.lm}
\title{Forecast a linear model with possible time series components}
\usage{\method{forecast}{lm}(object, newdata, h=10, level=c(80,95), fan=FALSE, 
    lambda=object$lambda, ...)
}

\arguments{
\item{object}{Object of class "lm", usually the result of a call to \code{\link[stats]{lm}} or \code{\link{tslm}}.}
\item{newdata}{An optional data frame in which to look for variables with which to predict. If omitted, it is assumed that the only variables are trend and season, and \code{h} forecasts are produced.}
\item{level}{Confidence level for prediction intervals.}
\item{fan}{If TRUE, level is set to seq(50,99,by=1). This is suitable for fan plots.}
\item{h}{Number of periods for forecasting. Ignored if \code{newdata} present.}
\item{lambda}{Box-Cox transformation parameter. Ignored if NULL. Otherwise, forecasts back-transformed via an inverse Box-Cox transformation.}
\item{...}{Other arguments passed to \code{\link[stats]{predict.lm}()}.}
}

\description{\code{forecast.lm} is used to predict linear models, especially those involving trend and seasonality components.}

\details{\code{forecast.lm} is largely a wrapper for \code{\link[stats]{predict.lm}()} except that it allows variables "trend" and "season" which are created on the fly from the time series characteristics of the data. Also, the output is reformatted into a \code{forecast} object.}


\value{An object of class "\code{forecast}".

The function \code{summary} is used to obtain and print a summary of the
results, while the function \code{plot} produces a plot of the forecasts and prediction intervals.

The generic accessor functions \code{fitted.values} and \code{residuals} extract useful features of
the value returned by \code{forecast.lm}.

An object of class \code{"forecast"} is a list containing at least the following elements:
\item{model}{A list containing information about the fitted model}
\item{method}{The name of the forecasting method as a character string}
\item{mean}{Point forecasts as a time series}
\item{lower}{Lower limits for prediction intervals}
\item{upper}{Upper limits for prediction intervals}
\item{level}{The confidence values associated with the prediction intervals}
\item{x}{The historical data for the response variable.}
\item{residuals}{Residuals from the fitted model. That is x minus fitted values.}
\item{fitted}{Fitted values}
}

\seealso{\code{\link{tslm}}, \code{\link[stats]{lm}}.}

\author{Rob J Hyndman}

\examples{
y <- ts(rnorm(120,0,3) + 1:120 + 20*sin(2*pi*(1:120)/12), frequency=12)
fit <- tslm(y ~ trend + season)
plot(forecast(fit, h=20))
}
\keyword{stats}

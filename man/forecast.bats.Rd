\name{forecast.bats}
\alias{forecast.bats}
\alias{forecast.tbats}
\title{Forecasting using BATS and TBATS models}

\usage{forecast(object, h=10, level=c(80,95), fan=FALSE, ts.frequency=max(object$seasonal.periods))}

\arguments{
\item{object}{A \code{bats} or \code{tbats} object. The the result of a call to \code{\link{bats}}. }
\item{object}{An object of class "\code{bats}". Usually the result of a call to \code{\link{bats}}.}
\item{h}{Number of periods for forecasting}
\item{level}{Confidence level for prediction intervals.}
\item{fan}{If TRUE, level is set to seq(50,99,by=1). This is suitable for fan plots.}
}

\description{Forecasts \code{h} steps ahead with a BATS model. Prediction intervals are also produced.}

\value{An object of class "\code{forecast}".

The function \code{summary} is used to obtain and print a summary of the
results, while the function \code{plot} produces a plot of the forecasts and prediction intervals.

The generic accessor functions \code{fitted.values} and \code{residuals} extract useful features of
the value returned by \code{forecast.bats}.

An object of class \code{"forecast"} is a list containing at least the following elements:
\item{model}{A copy of the \code{bats} object}
\item{method}{The name of the forecasting method as a character string}
\item{mean}{Point forecasts as a time series}
\item{lower}{Lower limits for prediction intervals}
\item{upper}{Upper limits for prediction intervals}
\item{level}{The confidence values associated with the prediction intervals}
\item{x}{The original time series (either \code{object} itself or the time series used to create the model stored as \code{object}).}
\item{residuals}{Residuals from the fitted model. That is x minus fitted values.}
\item{fitted}{Fitted values (one-step forecasts)}
}

\references{De Livera, A.M., Hyndman, R.J., & Snyder, R. D. (2010), Forecasting time series with complex seasonal patterns using exponential smoothing, Journal of the American Statistical Association, to appear.
}

\seealso{\code{\link{bats}}, \code{forecast.ets}.}

\author{Rob J Hyndman and Slava Razbash}

\examples{
fit<-bats(air.passengers)
plot(forecast(fit))
\dontrun{
taylor.fit<-bats(taylor)
plot(forecast(taylor.fit))
}
}


\keyword{ts}

\name{msts}
\alias{msts}

\title{Multi-Seasonal Time Series}
\usage{msts(data, seasonal.periods, ts.frequency=floor(max(seasonal.periods)), ... )}

\description{msts is an S3 class for multi seasonal time series objects, intended to be used for models that support multiple seasonal periods. The msts class inherits from the ts class and has an additional "msts" attribute which contains the vector of seasonal periods. All methods that work on a ts class, should also work on a msts class.}


\arguments{
\item{data}{A numeric vector, ts object, matrix or data frame. It is intended that the time series data is univariate, otherwise treated the same as ts().}
\item{seasonal.periods}{A vector of the seasonal periods of the msts.}
\item{ts.frequency}{The seasonal periods that should be used as frequency of the underlying ts object. The default value is \code{max(seasonal.periods)}.}
\item{...}{Arguements to be passed to the underlying call to \code{ts()}. For example \code{start=c(1987,5)}.}
}

\value{An object of class \code{c("msts", "ts")}.}

\author{Rob J Hyndman and Slava Razbash}

\examples{air.passengers<-msts(AirPassengers, seasonal.periods=c(4,12,18), ts.frequency=12, start=start(AirPassengers) ) }

\keyword{ts}

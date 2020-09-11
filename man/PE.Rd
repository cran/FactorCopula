\name{PE}
\alias{PE}
\docType{data}
\title{
Political-economic risk of 62 countries for the year 1987
}
\description{
Quinn (2004) used 5 mixed variables, namely the continuous variable black-market premium in each country (used as a proxy for illegal economic activity), the continuous variable productivity as measured by real gross domestic product per worker in 1985 international prices, the binary variable independence of the national judiciary (1 if the judiciary is judged to be independent and 0 otherwise), and the ordinal variables measuring the lack of expropriation risk and lack of corruption.
}
\usage{data(PE)}
\format{
  A data frame with 62 observations (countries) on the following 5 variables:
  \describe{
    \item{\code{BM}}{Black-market premium.}
    \item{\code{GDP}}{Gross domestic product.}
    \item{\code{IJ}}{Independent judiciary.}
    \item{\code{XPR}}{Lack of expropriation risk.}
    \item{\code{CPR}}{Lack of corruption.}
  } 
}
\source{

Quinn, K. M. (2004). 
Bayesian factor analysis for mixed ordinal and continuous responses. 
\emph{Political Analysis}, \bold{12}, 338--353.
}

\keyword{datasets}

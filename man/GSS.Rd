\name{GSS}
\alias{GSS}
\docType{data}
\title{
The 1994 General Social Survey
}
\description{
Hoff (2007) analysed seven demographic variables of 464 male respondents to the 1994 General Social Survey. 
Of these seven, two were continuous (income and age of the respondents), three were ordinal with 5 categories (highest degree of the survey respondent, income and highest degree of  respondent's parents),  and two were count variables (number of children of the survey respondent and respondent's parents). 
}
\usage{data(GSS)}
\format{
  A data frame with 464 observations on the following 7 variables:
  \describe{
    \item{\code{INCOME}}{Income of the respondent in 1000s of dollars, binned into 21 ordered categories.}
    \item{\code{DEGREE}}{Highest degree ever obtained (0:None, 1:HS, 2:Associates, 3:Bachelors, 4:Graduate).}
    \item{\code{CHILDREN}}{Number of children of the survey respondent.}
    \item{\code{PINCOME}}{Financial status of respondent's parents when respondent was 16 (on a 5-point scale).}
    \item{\code{PDEGREE}}{Highest degree of the survey respondent's parents (0:None, 1:HS, 2:Associates, 3:Bachelors, 4:Graduate).}
    \item{\code{PCHILDREN}}{Number of children of the survey respondent's parents - 1.}
    \item{\code{AGE}}{Age of the respondents in years.}
  } 
}
\source{

Hoff, P. D. (2007).
Extending the rank likelihood for semiparametric copula estimation.
\emph{The Annals of Applied Statistics}, \bold{1}, 265--283.
}

\keyword{datasets}

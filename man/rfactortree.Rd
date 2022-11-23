\name{rFactorTree}
\alias{rFactorTree}
\alias{r1factortree}
\alias{r2factortree}
\title{Simulation of 1- and 2-factor tree copula models for item response data}
\description{Simulating item response data from the 1- and 2-factor tree copula models.}
\usage{
r1factortree(n, d, A, copname1, copnametree, theta1, delta,K)
r2factortree(n, d, A, copname1, copname2, copnametree,theta1, theta2, delta,K)
}

\arguments{
\item{n}{Sample size.}

\item{d}{Number of observed variables/items.}

\item{A}{\eqn{d \times d} vine array with \eqn{1,...,d} on diagonal, note only the first row and diagnoal values are used for the 1-truncated vine model} 

\item{theta1}{copula parameter vector of size \eqn{d} for items with the first factor.}

\item{theta2}{copula parameter vector of size \eqn{d} for items with the second factor.}

\item{delta}{copula parameter vector of size \eqn{d-1} for the 1-truncated vine tree (conditional dependence).}

\item{copname1}{A name of a bivariate copula that link each of the oberved variabels with the first factor (note only a single copula family for all items with the factor). Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel.}

\item{copname2}{A name of a bivariate copula that link each of the oberved variabels with the second factor (note only a single copula family for all items with the factor). Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel.}

\item{copnametree}{A name of a bivariate copula that link each of the oberved variabels with one another given the factors in the 1-truncated vine (note only a single copula family for all tree). Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel.}
 
\item{K}{Number of categories for the observed variables/items.}
}


\value{
Data matrix of dimension \eqn{n \times d}, where \eqn{n} is the sample size, and \eqn{d} is the total number of observed variables/items.
}
\references{
Joe, H. (2014).
\emph{Dependence Modelling with Copulas}. Chapman & Hall, London.

Kadhem, S.H. and Nikoloulopoulos, A.K. (2022b)
Factor tree copula models for item response data.
\emph{Arxiv e-prints}, <arXiv: 2201.00339>. \url{https://arxiv.org/abs/2201.00339}.

}

\author{
Sayed H. Kadhem \cr
Aristidis K. Nikoloulopoulos \email{a.nikoloulopoulos@uea.ac.uk}
}

\examples{
# ---------------------------------------------------
# ---------------------------------------------------
#Sample size
n = 500

#Ordinal Variables  ---------------------------------
d = 5

#Categories for ordinal  ----------------------------
K = 5
# ---------------------------------------------------
#              1-2-factor tree copula model
# ---------------------------------------------------
#Copula parameters
theta1 = rep(3, d)
theta2 = rep(2, d)
delta = rep(1.5, d-1)

#Copula names
copulaname_1f = "gum"
copulaname_2f = "gum"
copulaname_vine = "gum"

#vine array
#Dvine
d=5
A=matrix(0,d,d)
A[1,]=c(1,c(1:(d-1)))
diag(A)=1:d



#----------------- Simulating data ------------------
#1-factor tree copula
data_1ft = r1factortree(n, d, A, copulaname_1f, copulaname_vine, 
theta1, delta,K)
#2-factor tree copula
data_2ft = r2factortree(n, d, A, copulaname_1f, copulaname_2f, 
copulaname_vine, theta1,theta2, delta,K)
}

\keyword{datagen}


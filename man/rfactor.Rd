\name{rfactor}
\alias{rfactor}
\alias{r1factor}
\alias{r2factor}
\title{Simulation of factor copula models for mixed continuous and discrete data}
\description{Simulating dependent standard uniform and ordinal response data from factor copula models.}
\usage{
r1factor(n, d1, d2, categ, theta, copF1)
r2factor(n, d1, d2, categ, theta, delta, copF1, copF2)
}
\arguments{

\item{n}{Sample size.}

\item{d1}{Number of standard uniform variables.}

\item{d2}{Number of ordinal variables.}

\item{categ}{A vector of categories for the ordinal variables.}

\item{theta}{Copula parameters for the 1st factor.}
\item{delta}{Copula parameters for the 2nd factor.}

\item{copF1}{\eqn{(d_1+d_2)}-vector with the names of bivariate copulas that link the each of the oberved variabels with the 1st factor. Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel, \dQuote{joe} for Joe, \dQuote{rjoe} for reflected Joe, \dQuote{1rjoe} for 1-reflected Joe, \dQuote{2rjoe} for 2-reflected Joe, \dQuote{BB1} for BB1, \dQuote{rBB1} for reflected BB1, \dQuote{BB7} for BB7, \dQuote{rBB7} for reflected BB7, \dQuote{BB8} for BB8, \dQuote{rBB8} for reflected BB8, \dQuote{BB10} for BB10, \dQuote{rBB10} for reflected BB10.}

\item{copF2}{\eqn{(d_1+d_2)}-vector with the names of bivariate copulas that link the each of the oberved variabels with the 2nd factor. Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel, \dQuote{joe} for Joe, \dQuote{rjoe} for reflected Joe, \dQuote{1rjoe} for 1-reflected Joe, \dQuote{2rjoe} for 2-reflected Joe, \dQuote{BB1} for BB1, \dQuote{rBB1} for reflected BB1, \dQuote{BB7} for BB7, \dQuote{rBB7} for reflected BB7, \dQuote{BB8} for BB8, \dQuote{rBB8} for reflected BB8, \dQuote{BB10} for BB10, \dQuote{rBB10} for reflected BB10.}
}


\value{
Data matrix of dimension \eqn{n\times d}, where \eqn{n} is the sample size, and \eqn{d=d_1+d_2} is the total number of variables.
}
\references{
Kadhem, S.H. and Nikoloulopoulos, A.K. (2019)
Factor copula models for mixed data.
\emph{Arxiv e-prints}, <arXiv:1907.07395>. \url{https://arxiv.org/abs/1907.07395}.

}
\author{
Sayed H. Kadhem \email{s.kadhem@uea.ac.uk}\cr
Aristidis K. Nikoloulopoulos \email{a.nikoloulopoulos@uea.ac.uk}
}

\examples{


# ---------------------------------------------------
# ---------------------------------------------------
#              One-factor copula model
# ---------------------------------------------------
# ---------------------------------------------------
#Sample size ----------------------------------------
n = 100

#Continuous Variables  ------------------------------
d1 = 5

#Ordinal Variables  ---------------------------------
d2 = 3

#Categories for ordinal  ----------------------------
categ = c(3,4,5)

#Copula parameters  ---------------------------------
theta = rep(2, d1+d2)

#Copula names  --------------------------------------
copnamesF1 = rep("gum", d1+d2)

#----------------- Simulating data ------------------
datF1 = r1factor(n, d1=d1, d2=d2, categ, theta, copnamesF1)

#------------  Plotting continuous data -------------
pairs(qnorm(datF1[, 1:d1]))

# ---------------------------------------------------
# ---------------------------------------------------
#              Two-factor copula model
# ---------------------------------------------------
# ---------------------------------------------------
#Sample size ----------------------------------------
n = 100

#Continuous Variables  ------------------------------
d1 = 5

#Ordinal Variables  ---------------------------------
d2 = 3

#Categories for ordinal  ----------------------------
categ = c(3,4,5)

#Copula parameters  ---------------------------------
theta = rep(2.5, d1+d2)
delta = rep(1.5, d1+d2)

#Copula names  --------------------------------------
copnamesF1 = rep("gum", d1+d2)
copnamesF2 = rep("gum", d1+d2)

#----------------- Simulating data ------------------
datF2 = r2factor(n, d1=d1, d2=d2, categ, theta, delta, 
                copnamesF1, copnamesF2)

#-----------------  Plotting  data ------------------
pairs(qnorm(datF2[,1:d1]))


}

\keyword{datagen}


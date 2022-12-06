\name{rStructuredFactor}
\alias{rStructuredFactor}
\alias{rBifactor}
\alias{rSecond_order}
\title{Simulation of bi-factor and second-order copula models for item response data}
\description{Simulating item response data from the bi-factor and second-order copula models for item response data.}
\usage{
rBifactor(n, d, grpsize, categ, copnames1,copnames2, theta1, theta2)
rSecond_order(n, d, grpsize, categ, copnames1, copnames2, theta1, theta2)
}

\arguments{
\item{n}{Sample size.}

\item{d}{Number of observed variables/items.}

\item{grpsize}{vector indicating the size for each group, e.g., c(4,4,4) indicating four items in all three groups.}

\item{categ}{A vector of categories for the observed variables/items.}

\item{theta1}{\bold{For the bi-factor model}: copula parameter vector of size \eqn{d} for items with the common factor. \bold{For the second-order copulas}: copula parameter vector of size \eqn{G} for the common factor and group-specific factors.}
\item{theta2}{\bold{For the bi-factor model}: copula parameter vector of size \eqn{d} for items with the group-specific factor. \bold{For the second-order copulas}: copula parameter vector of size \eqn{d} for items with the group-specific factor. }

\item{copnames1}{\bold{For the bi-factor copula}: \eqn{d}-vector with the names of bivariate copulas that link the each of the oberved variabels with the common factor. \bold{For the second-order factor copula}: \eqn{G}-vector with the names of bivariate copulas that link the each of the group-specific factors with the common factor, where \eqn{G} is the number of groups of items. Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel.}
 
\item{copnames2}{\bold{For the bi-factor copula}: \eqn{d}-vector with the names of bivariate copulas that link the each of the oberved variabels with the group-specific factor. \bold{For the second-order factor copula}: \eqn{d}-vector with the names of bivariate copulas that link the each of the oberved variabels with the group-specific factor. Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel.}
}


\value{
Data matrix of dimension \eqn{n\times d}, where \eqn{n} is the sample size, and \eqn{d} is the total number of observed variables/items.
}
\references{
Kadhem, S.H. and Nikoloulopoulos, A.K. (2022)
Bi-factor and second-order copula models for item response data.
\emph{Psychometrika}. \doi{doi:10.1007/s11336-022-09894-2}.
}
\author{
Sayed H. Kadhem \email{s.kadhem@uea.ac.uk}\cr
Aristidis K. Nikoloulopoulos \email{a.nikoloulopoulos@uea.ac.uk}
}

\examples{
# ---------------------------------------------------
# ---------------------------------------------------
#Sample size
n = 500

#Ordinal Variables  ---------------------------------
d = 9
grpsize=c(3,3,3)
ngrp=length(grpsize)

#Categories for ordinal  ----------------------------
categ = rep(3,d)

# ---------------------------------------------------
# ---------------------------------------------------
#              Bi-factor copula model
# ---------------------------------------------------
# ---------------------------------------------------
#Copula parameters
theta = rep(2.5, d)
delta = rep(1.5, d)

#Copula names
copulanames1 = rep("gum", d)
copulanames2 = rep("gum", d)

#----------------- Simulating data ------------------
data_Bifactor = rBifactor(n, d, grpsize, categ, copulanames1,
copulanames2, theta, delta)

# ---------------------------------------------------
# ---------------------------------------------------
#              Second-order copula model
# ---------------------------------------------------
# ---------------------------------------------------
#Copula parameters
theta= rep(1.5, ngrp)
delta = rep(2.5, d)

#Copula names
copulanames1 = rep("gum", ngrp)
copulanames2 = rep("gum", d)

data_Second_order = rSecond_order(n, d, grpsize, categ,
copulanames1, copulanames2, theta, delta)
}

\keyword{datagen}


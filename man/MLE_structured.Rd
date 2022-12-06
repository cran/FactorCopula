\name{mle.StructuredFactor}
\alias{mle.StructuredFactor}
\alias{mleBifactor}
\alias{mleSecond_order}
\title{Maximum likelihood estimation of the bi-factor and second-order copula models for item response data}
\description{We approach the estimation  of the bi-factor and second-order copula models for item response data with the IFM method of Joe (2005).}
\usage{
mleBifactor(y, copnames1, copnames2, gl, ngrp, grpsize,
hessian, print.level)
mleSecond_order(y, copnames1, copnames2, gl, ngrp, grpsize,
hessian, print.level)
}
\arguments{
\item{y}{\eqn{n \times d} matrix with the item reponse data, where \eqn{n} and \eqn{d} is the number of observations and variables, respectively.}

\item{copnames1}{\bold{For the bi-factor copula}: \eqn{d}-vector with the names of bivariate copulas that link the each of the oberved variabels with the common factor. \bold{For the second-order factor copula}: \eqn{G}-vector with the names of bivariate copulas that link the each of the group-specific factors with the common factor, where \eqn{G} is the number of groups of items. Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel.}
 
\item{copnames2}{\bold{For the bi-factor copula}: \eqn{d}-vector with the names of bivariate copulas that link the each of the oberved variabels with the group-specific factor. \bold{For the second-order factor copula}: \eqn{d}-vector with the names of bivariate copulas that link the each of the oberved variabels with the group-specific factor. Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel.}

\item{gl}{Gauss legendre quardrature nodes and weights.}

\item{ngrp}{number of non-overlapping groups.}

\item{grpsize}{vector indicating the size for each group, e.g., c(4,4,4) indicating four items in all three groups.}

\item{hessian}{If TRUE, the hessian of the negative log-likelihood is calculated during the minimization process.}

\item{print.level}{Determines the level of printing which is done during the minimization process; same as in \code{nlm}.}
}


\value{A list containing the following components:
\item{cutpoints}{The estimated univariate cutpoints.}
\item{taus}{The estimated copula parameters in Kendall's tau scale.}
\item{SEs}{The SEs of the Kendall's tau estimates.}
\item{loglik}{The maximized joint log-likelihood.}
}


\details{Estimation is achieved by maximizing  the joint log-likelihood
over the copula parameters  with the univariate cutpoints fixed as estimated at the first step of the proposed two-step estimation approach.
}


\references{
Joe, H. (2005)
Asymptotic efficiency of the two-stage estimation method for copula-based models.
\emph{Journal of Multivariate Analysis}, \bold{94}, 401--419.
\doi{10.1016/j.jmva.2004.06.003}.

Kadhem, S.H. and Nikoloulopoulos, A.K. (2022)
Bi-factor and second-order copula models for item response data.
\emph{Psychometrika}. \doi{doi:10.1007/s11336-022-09894-2}.
}
\author{
Sayed H. Kadhem \email{s.kadhem@uea.ac.uk}\cr
Aristidis K. Nikoloulopoulos \email{a.nikoloulopoulos@uea.ac.uk}
}


\examples{
\donttest{
#------------------------------------------------
# Setting quadreture points
nq <- 25
gl <- gauss.quad.prob(nq)
#------------------------------------------------
#                     TAS Data
#------------------             -----------------
data(TAS)
#using a subset of the data
#group1: 1,3,6,7,9,13,14
grp1=c(1,3,6)
#group2: 2,4,11,12,17
grp2=c(2,4,11)
#group3: 5,8,10,15,16,18,19,20
grp3=c(5,8,10)
#Rearrange items within testlets
set.seed(123)
i=sample(1:nrow(TAS),500)
ydat=TAS[i,c(grp1,grp2,grp3)]

d=ncol(ydat);d
n=nrow(ydat);n

#size of each group
g1=length(grp1)
g2=length(grp2)
g3=length(grp3)

grpsize=c(g1,g2,g3)#group size
#number of groups
ngrp=length(grpsize)

#BI-FACTOR
copX0 = rep("bvt2", d)
copXg = c(rep("rgum", g1), rep("bvt3", g2+g3))
mle_Bifactor =  mleBifactor(y = ydat, copX0, copXg, gl, ngrp, grpsize, hessian=FALSE, print.level=2)

}
}
\keyword{models}
\keyword{multivariate}
\keyword{maximum likelihood}



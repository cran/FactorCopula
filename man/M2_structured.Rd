\name{M2.StructuredFactor}
\alias{M2.StructuredFactor}
\alias{M2Bifactor}
\alias{M2Second_order}
\title{Goodness-of-fit of bi-factor and second-order copula models for item response data}
\description{ The limited information \eqn{M_2} statistic (Maydeu-Olivares and Joe, 2006)  of bi-factor and second-order copula models for item response data.}
\usage{
M2Bifactor(y,cpar, copnames1, copnames2, gl, ngrp, grpsize)
M2Second_order(y,cpar, copnames1, copnames2, gl, ngrp, grpsize)
}
\arguments{
\item{y}{\eqn{n \times d} matrix with the ordinal reponse data, where \eqn{n} and \eqn{d} is the number of observations and variables, respectively.}

\item{cpar}{A list of estimated copula parameters.}

\item{copnames1}{\bold{For the bi-factor copula}: \eqn{d}-vector with the names of bivariate copulas that link each of the oberved variabels with the common factor. \bold{For the second-order factor copula}: \eqn{G}-vector with the names of bivariate copulas that link the each of the group-specific factors with the common factor, where \eqn{G} is the number of groups of items. Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{2, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel.}
 
\item{copnames2}{\bold{For the bi-factor copula}: \eqn{d}-vector with the names of bivariate copulas that link the each of the oberved variabels with the group-specific factor. \bold{For the second-order factor copula}: \eqn{d}-vector with the names of bivariate copulas that link the each of the oberved variabels with the group-specific factor. Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{2, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel.}


\item{gl}{Gauss legendre quardrature nodes and weights.}

\item{ngrp}{number of non-overlapping groups.}

\item{grpsize}{vector indicating the size for each group, e.g., c(4,4,4) indicating four items in all three groups.}

}

\value{A list containing the following components:

\item{M2}{The \eqn{M_2} statistic which has a null asymptotic distribution that is \eqn{\chi^2} with \eqn{s-q} degrees of freedom, where \eqn{s} is the number of univariate and bivariate margins that do not include the category 0 and \eqn{q} is the number of model parameters.}
\item{df}{\eqn{s-q}.}
\item{p-value}{The resultant \eqn{p}-value.}
}



\details{The \eqn{M_2} statistic has been developed for goodness-of-fit testing in multidimensional contingency tables by Maydeu-Olivares and Joe (2006). We use the \eqn{M_2} to assess the overall fit for the bi-factor and second-order copula models for item resposne data (Kadhem & Nikoloulopoulos, 2021).
}

\references{
Kadhem, S.H. and Nikoloulopoulos, A.K. (2023)
Bi-factor and second-order copula models for item response data.
\emph{Psychometrika}, \bold{88}, 132--157. \doi{doi:10.1007/s11336-022-09894-2}.

Maydeu-Olivares, A. and Joe, H. (2006).
Limited information goodness-of-fit testing in multidimensional contingency tables.
\emph{Psychometrika}, \bold{71}, 713--732. \doi{10.1007/s11336-005-1295-9}.

}
\author{
Sayed H. Kadhem \email{s.kadhem@uea.ac.uk}\cr
Aristidis K. Nikoloulopoulos \email{a.nikoloulopoulos@uea.ac.uk}
}


\examples{
\donttest{
#------------------------------------------------
# Setting quadreture points
nq <- 15
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

#------------------------------------------------
#                       M2
#------------------------------------------------
#BI-FACTOR
tauX0 = c(0.49,0.16,0.29,#0.09,0.47,0.49,0.30,
          0.46,0.41,0.33,#0.29,0.24,
          0.10,0.16,0.14)#,0.12,0.03,0.03,0.10,0.10)
tauXg = c(0.09,0.37,0.23,#0.53,0.24,0.32,0.27,
          0.53,0.58,0.20,#0.23,0.25,0.34,0.33,
          0.30,0.19,0.24)#,0.29,0.43,0.26)
copX0 = rep("bvt2", d)
copXg = c(rep("rgum", g1), rep("bvt3", g2+g3))
#converting taus to cpars
cparX0=mapply(function(x,y) tau2par(x,y),x=copX0,y=tauX0)
cparXg=mapply(function(x,y) tau2par(x,y),x=copXg,y=tauXg)
cpar=c(cparX0,cparXg)

m2_Bifactor = M2Bifactor(y=ydat, cpar, copX0, copXg, gl, ngrp, grpsize)

}
}
\keyword{multivariate}



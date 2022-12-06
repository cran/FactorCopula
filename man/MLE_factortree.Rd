\name{mle.FactorTree}
\alias{mle.FactorTree}
\alias{mle1FactorTree}
\alias{mle2FactorTree}
\title{Maximum likelhood estimation of factor tree copula models}
\description{We use a two-stage estimation approach  toward the estimation  of  factor tree copula models for item response data.}
\usage{
mle1FactorTree(y, A, cop, gl, hessian, print.level) 
mle2FactorTree(y, A, cop, gl, hessian, print.level) 
}
\arguments{
\item{y}{\eqn{n \times d} matrix with the ordinal reponse data, where \eqn{n} and \eqn{d} is the number of observations and ordinal variables, respectively.}

\item{A}{\eqn{d \times d} vine array with \eqn{1,...,d} on diagonal, note only the first row and diagnoal values are used for the 1-truncated vine model} 

\item{cop}{\eqn{(2d-1)}-vector with the names of bivariate copulas that link each of the oberved variabels with the 1st factor (1-factor part of the model), and conditional dependence of variables given the latent factor (1-truncated vine tree part of the model) . Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel.}

\item{gl}{Gauss legendre quardrature nodes and weights.}


\item{hessian}{If TRUE, the hessian of the negative log-likelihood is calculated during the minimization process.}

\item{print.level}{Determines the level of printing which is done during the minimization process; same as in \code{nlm}.}
}


\value{A list containing the following components:
\item{cutpoints}{The estimated univariate cutpoints.}
\item{loglik}{The maximized joint log-likelihood.}
\item{taus}{The estimated copula parameters in Kendall's tau scale.}
\item{SEs}{The SEs of the Kendall's tau estimates.}
}


\details{Estimation is achieved by maximizing  the joint log-likelihood 
over the copula parameters  with the univariate cutpoints fixed as estimated at the first step of the proposed two-step estimation approach. 
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
\donttest{
#------------------------------------------------
# Setting quadreture points
nq <- 5  
gl <- gauss.quad.prob(nq) 
#------------------------------------------------
#                    PTSD Data
#------------------             -----------------
data(PTSD)
ydat=PTSD
n=nrow(ydat)
d=ncol(ydat)
#------------------------------------------------
#                   Estimation
#------------------             -----------------
#selecting vine tree based on polychoric
rmat=polychoric0(ydat)$p
A.polychoric=selectFactorTrVine(y=ydat,rmat,alg=3)

#---------------- 1-factor tree  ----------------
# 1-factor tree copula model
copf1 <- rep("frk",d)
coptree <- rep("frk",d-1)
cop <- c(copf1,coptree)
est1factortree <- mle1FactorTree(y=ydat, A=A.polychoric$VineTreeA, cop, 
gl, hessian=FALSE, print.level=2) 

#---------------- 2-factor tree  ----------------
# 2-factor tree copula model
copf1 <- rep("frk",d)
copf2 <- rep("frk",d)
coptree <- rep("frk",d-1)
cop <- c(copf1,copf2,coptree)

est2factortree <- mle2FactorTree(y=ydat, A=A.polychoric$VineTreeA, 
cop, gl, hessian=FALSE, print.level=2)
     
}
}
\keyword{models}
\keyword{multivariate}
\keyword{maximum likelihood}



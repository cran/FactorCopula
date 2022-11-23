\name{Vuong.FactorTree}
\alias{Vuong.FactorTree}
\alias{vuong_FactorTree}
\title{
Vuong's test for the comparison of 1- and 2-factor tree copula models for item response data}
\description{
The Vuong's test (Vuong,1989) is the sample version of the difference in Kullback-Leibler divergence between two models and can be used to differentiate two  parametric models which could be non-nested. 
For the Vuong's test we provide the 95\% confidence interval of the Vuong's test statistic (Joe, 2014, page 258). If the interval does not contain  0, then the best fitted model according to the AICs is better if the interval is completely above 0.
}
\usage{
vuong_FactorTree(models, y, A.m1, tau.m1, copnames.m1, 
tau.m2, copnames.m2,A.m2,nq)
}
\arguments{
\item{models}{choose a number from (1,2,3,4,5,6,7). 
\bold{1}: Model 1 is 1-factor tree copula, Model 2 is 1-factor tree copula.
\bold{2}: Model 1 is 2-factor tree copula, Model 2 is 2-factor tree copula.
\bold{3}: Model 1 is 1-factor tree copula, Model 2 is 2-factor tree copula.
\bold{4}: Model 1 is 1-factor copula, Model 2 is 1-factor tree copula.
\bold{5}: Model 1 is 1-factor copula, Model 2 is 2-factor tree copula.
\bold{6}: Model 1 is 2-factor copula, Model 2 is 1-factor tree copula.
\bold{7}: Model 1 is 2-factor copula, Model 2 is 2-factor tree copula.
}

\item{y}{\eqn{n \times d} matrix with the item reponse data, where \eqn{n} and \eqn{d} is the number of observations and variables, respectively.}

\item{A.m1}{\eqn{d \times d} vine array for Model 1, if it is 1-factor tree, or 2-factor tree otherwise keep as NULL.} 

\item{A.m2}{\eqn{d \times d} vine array for Model 2, note only the first row and diagnoal values are used for the 1-truncated vine model.} 

\item{tau.m1}{vector of copula paramters in Kendall's \eqn{\tau} for model 1.}

\item{tau.m2}{vector of copula paramters in Kendall's \eqn{\tau} for model 2.}

\item{copnames.m1}{vector of names of copula families for model 1.}

\item{copnames.m2}{vector of names of copula families for model 2.}

\item{nq}{Number of Gauss legendre quardrature points.}
}

\value{
A vector containing the following components:
\item{z}{the test statistic.}
\item{p.value}{ the \eqn{p}-value.}
\item{CI.left}{ lower/left endpoint of 95\% confidence interval.}
\item{CI.right}{ upper/right endpoint of 95\% confidence interval.}
}


\references{

Joe, H. (2014).
\emph{Dependence Modelling with Copulas}. Chapman and Hall/CRC.

Kadhem, S.H. and Nikoloulopoulos, A.K. (2022b)
Factor tree copula models for item response data.
\emph{Arxiv e-prints}, <arXiv: 2201.00339>. \url{https://arxiv.org/abs/2201.00339}.

Vuong, Q.H. (1989).
Likelihood ratio tests for model selection and non-nested hypotheses.
\emph{Econometrica}, \bold{57}, 307--333. 
}

\author{
Sayed H. Kadhem \cr
Aristidis K. Nikoloulopoulos \email{a.nikoloulopoulos@uea.ac.uk}
}

\examples{
\donttest{
#------------------------------------------------
# Setting quadreture points
nq <- 25
#------------------------------------------------
#                    PTSD Data
#------------------             -----------------
data(PTSD)
ydat=PTSD

d=ncol(ydat);d
n=nrow(ydat);n

#vine tree structure
#selecting vine tree based on polychoric corr
rmat=polychoric0(ydat)$p
A.polychoric=selectFactorTrVine(y=ydat,rmat,alg=3)
A.polychoric=A.polychoric$VineTreeA
#------------------------------------------------
# M1 1-factor tree - M2 1-factor tree
tau.m1 = rep(0.4,d*2-1)
copnames.m1 = rep("bvn",d*2-1)
tau.m2 = rep(0.4,d*2-1)
copnames.m2 = rep("rgum",d*2-1)
vuong.1factortree = vuong_FactorTree(models=1, ydat, 
A.m1=A.polychoric, tau.m1, copnames.m1, tau.m2, 
copnames.m2,A.m2=A.polychoric,nq)
#------------------------------------------------
# M1 2-factor tree - M2 2-factor tree
tau.m1 = rep(0.4,d*3-1)
copnames.m1 = rep("bvn",d*3-1)
tau.m2 = rep(0.4,d*3-1)
copnames.m2 = rep("rgum",d*3-1)

vuong.2factortree = vuong_FactorTree(models=2, ydat, 
A.m1=A.polychoric, tau.m1, copnames.m1, tau.m2, 
copnames.m2,A.m2=A.polychoric,nq)

#------------------------------------------------
# M1 1-factor tree - M2 2-factor tree
tau.m1 = rep(0.4,d*2-1)
copnames.m1 = rep("bvn",d*2-1)

tau.m2 = rep(0.4,d*3-1)
copnames.m2 = rep("rgum",d*3-1)

vuong.12factortree = vuong_FactorTree(models=3, ydat, 
A.m1=A.polychoric, tau.m1, copnames.m1, tau.m2, 
copnames.m2,A.m2=A.polychoric,nq)
                
}
}

\keyword{maximum likelihood}

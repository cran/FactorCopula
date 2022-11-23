\name{Vuong.StructuredFactor}
\alias{Vuong.StructuredFactor}
\alias{vuong_structured}
\title{
Vuong's test for the comparison of bi-factor and second-order copula models}
\description{
The Vuong's test (Vuong,1989) is the sample version of the difference in Kullback-Leibler divergence between two models and can be used to differentiate two  parametric models which could be non-nested. 
For the Vuong's test we provide the 95\% confidence interval of the Vuong's test statistic (Joe, 2014, page 258). If the interval does not contain  0, then the best fitted model according to the AICs is better if the interval is completely above 0.
}
\usage{
vuong_structured(models, cpar.m1, copnames.m1, cpar.m2, 
copnames.m2, y, grpsize)
}
\arguments{
\item{models}{choose a number from (1,2,3,4). \bold{1}: Model1 is bifactor, Model2 is bifactor. \bold{2}: Model1 is second-order, Model2 is second-order. \bold{3}: Model1 is second-order, Model2 is bifactor. \bold{4}: Model1 is bifactor, Model2 is nested.}

\item{cpar.m1}{vector of copula paramters for model 1, starting with copula parameters that link the items with common factor (bifactor), or group factors with common factor (second-order).}

\item{cpar.m2}{vector of copula paramters for model 2, starting with copula parameters that link the items with common factor (bifactor), or group factors with common factor (second-order).}

\item{copnames.m1}{vector of names of copula families for model 1, starting with copulas that link the items with common factor (bifactor), or group factors with common factor (second-order).}

\item{copnames.m2}{vector of names of copula families for model 2, starting with copulas that link the items with common factor (bifactor), or group factors with common factor (second-order).}

\item{y}{matrix of ordinal data.}

\item{grpsize}{vector indicating the size for each group, e.g., c(4,4,4) indicating four items in all three groups.}
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

Kadhem, S.H. and Nikoloulopoulos, A.K. (2021)
Bi-factor and second-order copula models for item response data.
\emph{Arxiv e-prints}, <arXiv:2102.10660>. \url{https://arxiv.org/abs/2102.10660}.

Vuong, Q.H. (1989).
Likelihood ratio tests for model selection and non-nested hypotheses.
\emph{Econometrica}, \bold{57}, 307--333. 
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
grp1=c(1,3,6,7,9,13,14)
grp2=c(2,4,11,12,17)
grp3=c(5,8,10,15,16,18,19,20)
ydat=TAS[,c(grp1,grp2,grp3)]

d=ncol(ydat);d
n=nrow(ydat);n

#Rearrange items within testlets
g1=length(grp1)
g2=length(grp2)
g3=length(grp3)

grpsize=c(g1,g2,g3)#group size
#number of groups
ngrp=length(grpsize)
#------------------------------------------------
# M1 bifactor - M2 bifactor
cpar.m1 = rep(0.6,d*2)
copnames.m1 = rep("bvn",d*2)
cpar.m2 = rep(1.6,d*2)
copnames.m2 = rep("rgum",d*2)
vuong.bifactor = vuong_structured(models=1, cpar.m1, copnames.m1, 
                 cpar.m2, copnames.m2, 
                 y=ydat, grpsize)

# M1 seconod-order - M2 seconod-order
cpar.m1 = rep(0.6,d+ngrp)
copnames.m1 = rep("bvn",d+ngrp)
cpar.m2 = rep(1.6,d+ngrp)
copnames.m2 = rep("rgum",d+ngrp)
vuong.second_order = vuong_structured(models=2, cpar.m1, 
    copnames.m1, cpar.m2, copnames.m2, y=ydat, grpsize)
 
# M1 seconod-order - M2 bifactor
cpar.m1 = rep(0.6,d+ngrp)
copnames.m1 = rep("bvn",d+ngrp)
cpar.m2 = rep(1.6,d*2)
copnames.m2 = rep("rgum",d*2)
vuong.2ndO_bif = vuong_structured(models=3, cpar.m1, copnames.m1, 
                 cpar.m2, copnames.m2, 
                 y=ydat, grpsize)     
                 
# M1 bifactor - M2 seconod-order
cpar.m1 = rep(0.6,d*2)
copnames.m1 = rep("bvn",d*2)
cpar.m2 = rep(1.6,d+ngrp)
copnames.m2 = rep("rgum",d+ngrp)
vuong.bif_2ndO = vuong_structured(models=4, cpar.m1, copnames.m1, 
                 cpar.m2, copnames.m2, 
                 y=ydat, grpsize)
}
}

\keyword{maximum likelihood}

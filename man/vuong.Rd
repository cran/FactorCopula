\name{vuong}
\alias{vuong}
\alias{vuong.1f}
\alias{vuong.2f}
\title{
Vuong's test for the comparison of factor  copula models}
\description{
Vuong (1989)'s test for the comparison of non-nested factor copula models for mixed data.  
We compute the Vuong's test between the factor copula model with BVN copulas (that is the standard factor  model) and a competing factor copula model to reveal if the latter provides better fit than the standard factor model.  
}
\usage{
vuong.1f(cpar.bvn, cpar, copF1, continuous, ordinal, count, gl, param)
vuong.2f(cpar.bvn, cpar, copF1, copF2, continuous, ordinal, count, gl, param)
}
\arguments{
\item{cpar.bvn}{copula parameters of the factor copula model with BVN copulas.}
\item{cpar}{copula parameters of the competing factor copula model.}
\item{copF1}{copula names for the first factor of the competing factor copula model.}
\item{copF2}{copula names for the second factor of the competing factor copula model.}
\item{continuous}{ matrix of continuous data.}
\item{ordinal}{ matrix of ordinal data.}
\item{count}{ matrix of count data.}

\item{gl}{ gauss-legendre quardature points.}
\item{param}{parameterization of estimated copula parameters. If FALSE, then cpar are the actual copula parameters without any transformation/reparamterization.}

}

\value{
A vector containing the following components:
\item{z}{the test statistic.}
\item{p.value}{ the \eqn{p}-value.}
\item{CI.left}{ lower/left endpoint of 95\% confidence interval.}
\item{CI.right}{ upper/right endpoint of 95\% confidence interval.}
}


\references{
Kadhem, S.H. and Nikoloulopoulos, A.K. (2019)
Factor copula models for mixed data.
\emph{Arxiv e-prints}, <arXiv:1907.07395>. \url{https://arxiv.org/abs/1907.07395}.

Vuong, Q.~H. (1989).
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
#                     PE Data
#------------------             -----------------
data(PE)
continuous.PE1 = -PE[,1]
continuous.PE2 = PE[,2]
continuous.PE <- cbind(continuous.PE1, continuous.PE2)
categorical.PE <- PE[, 3:5]
d <- ncol(PE)
#------------------------------------------------
#                   Estimation
#------------------             -----------------
# factor copula model with BVN copulas
cop1f.PE.bvn <- rep("bvn", d)
PE.bvn1f <- mle1factor(continuous.PE, categorical.PE, 
count=NULL, copF1=cop1f.PE.bvn, gl, hessian = T)

# Selected factor copula model
cop1f.PE <- c("joe", "joe", "rjoe", "joe", "gum")
PE.selected1f <- mle1factor(continuous.PE, categorical.PE, 
count=NULL, copF1=cop1f.PE, gl, hessian = T)
#------------------------------------------------
#                   Vuong's test
#------------------             -----------------
v1f.PE.selected <- vuong.1f(PE.bvn1f$cpar$f1,
PE.selected1f$cpar$f1,cop1f.PE, continuous.PE, 
categorical.PE, count=NULL, gl, param=F)

#------------------------------------------------
#------------------------------------------------
#                    GSS Data
#------------------             -----------------
data(GSS)
attach(GSS)
continuous.GSS <- cbind(INCOME,AGE)
ordinal.GSS <- cbind(DEGREE,PINCOME,PDEGREE) 
count.GSS <- cbind(CHILDREN,PCHILDREN)
d <- ncol(GSS)

#------------------------------------------------
#                   Estimation
#------------------             -----------------
# factor copula model with BVN copulas
# one-factor copula model
cop1f.GSS.bvn <- rep("bvn", d)
GSS.bvn1f <- mle1factor(continuous.GSS, ordinal.GSS, 
count.GSS, copF1 = cop1f.GSS.bvn, gl, hessian = T)

# two-factor copula model
cop1f.GSS.bvn = cop2f.GSS.bvn = rep("bvn", d)
GSS.bvn2f <- mle2factor.bvn(continuous.GSS, ordinal.GSS,
count.GSS, copF1 = cop1f.GSS.bvn, copF2 = cop2f.GSS.bvn, gl, SpC =7)

# Selected factor copula  model
# one-factor copula model
cop1f.GSS <- c("joe","2rjoe","bvt3","bvt3",
               "rgum","2rjoe","2rgum")
GSS.selected1f <- mle1factor(continuous.GSS, ordinal.GSS, 
count.GSS, copF1 = cop1f.GSS, gl, hessian = T)

# two-factor copula model
cop2f1.GSS <- c("rgum","rjoe","bvn","1rjoe","1rjoe","rjoe","gum")
cop2f2.GSS <- c("gum","2rjoe","rjoe","gum","bvt5","bvn","2rgum")
GSS.selected2f <- mle2factor(continuous.GSS, ordinal.GSS,
count.GSS, copF1 = cop2f1.GSS, copF2 = cop2f2.GSS, gl, hessian = T)

#------------------------------------------------
#                   Vuong's test
#------------------             -----------------
#1-factor
v1f.GSS.selected <- vuong.1f(GSS.bvn1f$cpar$f1, 
GSS.selected1f$cpar$f1, cop1f.GSS, continuous.GSS, 
ordinal.GSS, count=count.GSS, gl, param=F)

#2-factor
v2f.GSS.selected <- vuong.2f(GSS.bvn2f$cpar, 
GSS.selected2f$cpar, cop2f1.GSS, cop2f2.GSS, 
continuous.GSS, ordinal.GSS, count=count.GSS, gl, param=F)


}
}

\keyword{maximum likelihood}

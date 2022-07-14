\name{mle.Factor}
\alias{mle1factor}
\alias{mle2factor}
\alias{mle2factor.bvn}
\title{Maximum likelhood estimation of factor copula models for mixed data}
\description{We use a two-stage etimation approach  toward the estimation  of  factor copula models for mixed continuous and discrete data.}
\usage{
mle1factor(continuous, ordinal, count, copF1, gl, hessian, print.level)
mle2factor(continuous, ordinal, count, copF1, copF2, gl, hessian, print.level)
mle2factor.bvn(continuous, ordinal, count, copF1, copF2, gl, SpC, print.level)
}
\arguments{

\item{continuous}{\eqn{n \times d_1} matrix with the continuous reponse data, where \eqn{n} and \eqn{d_1} is the number of observations and continous  variables, respectively.}

\item{ordinal}{\eqn{n \times d_2} matrix with the ordinal reponse data, where \eqn{n} and \eqn{d_2} is the number of observations and ordinal variables, respectively.}

\item{count}{\eqn{n \times d_3} matrix with the count reponse data, where \eqn{n} and \eqn{d_3} is the number of observations and count variables, respectively.}

\item{copF1}{\eqn{(d_1+d_2+d_3)}-vector with the names of bivariate copulas that link the each of the oberved variabels with the 1st factor. Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel, \dQuote{joe} for Joe, \dQuote{rjoe} for reflected Joe, \dQuote{1rjoe} for 1-reflected Joe, \dQuote{2rjoe} for 2-reflected Joe, \dQuote{BB1} for BB1, \dQuote{rBB1} for reflected BB1, \dQuote{BB7} for BB7, \dQuote{rBB7} for reflected BB7, \dQuote{BB8} for BB8, \dQuote{rBB8} for reflected BB8, \dQuote{BB10} for BB10, \dQuote{rBB10} for reflected BB10.}

\item{copF2}{\eqn{(d_1+d_2+d_3)}-vector with the names of bivariate copulas that link the each of the oberved variabels with the 2nd factor. Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel, \dQuote{joe} for Joe, \dQuote{rjoe} for reflected Joe, \dQuote{1rjoe} for 1-reflected Joe, \dQuote{2rjoe} for 2-reflected Joe, \dQuote{BB1} for BB1, \dQuote{rBB1} for reflected BB1, \dQuote{BB7} for BB7, \dQuote{rBB7} for reflected BB7, \dQuote{BB8} for BB8, \dQuote{rBB8} for reflected BB8, \dQuote{BB10} for BB10, \dQuote{rBB10} for reflected BB10.}


\item{gl}{Gauss legendre quardrature nodes and weights.}

\item{SpC}{Special case for the 2-factor copula model with BVN copulas. Select a bivariate copula at the 2nd factor to be fixed to independence. e.g. "SpC = 1" to set the first copula at the 2nd factor to independence.  }
\item{hessian}{If TRUE, the hessian of the negative log-likelihood is calculated during the minimization process.}

\item{print.level}{Determines the level of printing which is done during the minimization process; same as in \code{nlm}.}



}


\value{A list containing the following components:
\item{cutpoints}{The estimated univariate cutpoints (fitting the univariate probit model).}
\item{negbinest}{The estimated univariate parametes for the count responses (fitting the negative binomial distribution).}
\item{loglik}{The maximized joint log-likelihood.}
\item{cpar}{Estimated copula parameters in a list form.}
\item{taus}{The estimated copula parameters in Kendall's tau scale.}
\item{SEs}{The SEs of the Kendall's tau estimates.}
}


\details{Estimation is achieved by maximizing  the joint log-likelihood 
over the copula parameters  with the univariate parameters/distributions fixed as estimated at the first step of the proposed two-step estimation approach. 
}


\references{
Kadhem, S.H. and Nikoloulopoulos, A.K. (2021)
Factor copula models for mixed data.
\emph{British Journal of Mathematical and Statistical Psychology}, \bold{74}, 365--403. \doi{10.1111/bmsp.12231}.

Krupskii, P. and Joe, H. (2013)
Factor copula models for multivariate data.
\emph{Journal of Multivariate Analysis}, \bold{120}, 85--101. \doi{10.1016/j.jmva.2013.05.001}.


Nikoloulopoulos, A.K. and Joe, H. (2015)
Factor copula models with item response data.
\emph{Psychometrika}, \bold{80}, 126--150. \doi{10.1007/s11336-013-9387-4}.

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
#------------------------------------------------
#                   Estimation
#------------------             -----------------
#------------------ One-factor  -----------------
# one-factor copula model
cop1f.PE <- c("joe", "joe", "rjoe", "joe", "gum")
est1factor.PE <- mle1factor(continuous.PE, categorical.PE, 
                            count=NULL, copF1=cop1f.PE, gl, hessian = TRUE)
est1factor.PE                    
#------------------------------------------------
#------------------------------------------------
#                     GSS Data
#------------------             -----------------
data(GSS)
attach(GSS)
continuous.GSS <- cbind(INCOME, AGE)
ordinal.GSS <- cbind(DEGREE, PINCOME, PDEGREE) 
count.GSS <- cbind(CHILDREN, PCHILDREN)

#------------------------------------------------
#                   Estimation
#------------------             -----------------
#------------------ One-factor  -----------------
# one-factor copula model
cop1f.GSS <- c("joe","2rjoe","bvt3","bvt3",
          "rgum","2rjoe","2rgum")
est1factor.GSS <- mle1factor(continuous.GSS, ordinal.GSS, 
                        count.GSS, copF1 = cop1f.GSS, gl, hessian = TRUE)

     
}
}
\keyword{models}
\keyword{multivariate}



\name{mapping}
\alias{par2tau}
\alias{tau2par}
\title{
Mapping of Kendall's tau and copula parameter
}
\description{
Bivariate copulas: mapping of Kendall's tau and copula parameter.
}
\usage{
par2tau(copulaname, cpar)
tau2par(copulaname, tau)
}
\arguments{
\item{copulaname}{Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel, \dQuote{joe} for Joe, \dQuote{rjoe} for reflected Joe, \dQuote{1rjoe} for 1-reflected Joe, \dQuote{2rjoe} for 2-reflected Joe, \dQuote{BB1} for BB1, \dQuote{rBB1} for reflected BB1, \dQuote{BB7} for BB7, \dQuote{rBB7} for reflected BB7, \dQuote{BB8} for BB8, \dQuote{rBB8} for reflected BB8, \dQuote{BB10} for BB10, \dQuote{rBB10} for reflected BB10.}
\item{cpar}{Copula parameter(s).}
\item{tau}{Kendall's tau.}
}

\value{
Kendall's tau or copula parameter. 
}

\references{
Joe H (1997) \emph{Multivariate Models and Dependence Concepts}. Chapman & Hall

Joe H (2014) \emph{Dependence Modeling with Copulas}. Chapman & Hall/CRC.

Joe H (2014) \emph{CopulaModel:  Dependence Modeling with Copulas}. Software for book: \emph{Dependence Modeling with Copulas}, Chapman & Hall/CRC, 2014.
}

\examples{


# 1-param copulas
#BVN copula
cpar.bvn = tau2par("bvn", 0.5)
tau.bvn = par2tau("bvn", cpar.bvn)

#Frank copula
cpar.frk = tau2par("frk", 0.5)
tau.frk = par2tau("frk", cpar.frk)

#Gumbel copula
cpar.gum = tau2par("gum", 0.5)
tau.gum = par2tau("gum", cpar.gum)

#Joe copula
cpar.joe = tau2par("joe", 0.5)
tau.joe = par2tau("joe", cpar.joe)

# 2-param copulas
#BB1 copula
tau.bb1 = par2tau("bb1", c(0.5,1.5))

#BB7 copula
tau.bb7 = par2tau("bb7", c(1.5,1))

#BB8 copula
tau.bb8 = par2tau("bb8", c(3,0.8))

#BB10 copula
tau.bb10 = par2tau("bb10", c(3,0.8))


}


\keyword{parameters}

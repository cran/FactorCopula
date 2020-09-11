\name{select}
\alias{select}
\alias{select1F}
\alias{select2F}
\title{Model selection of the factor copula models for mixed data}
\description{An heuristic algorithm that automatically selects the bivariate parametric copula families that link the observed to the latent variables.}
\usage{
select1F(continuous, ordinal, count, copnamesF1, gl)
select2F(continuous, ordinal, count, copnamesF1, copnamesF2, gl)
}
\arguments{

\item{continuous}{\eqn{n \times d_1} matrix with the continuous reponse data, where \eqn{n} and \eqn{d_1} is the number of observations and continous  variables, respectively.}

\item{ordinal}{\eqn{n \times d_2} matrix with the ordinal reponse data, where \eqn{n} and \eqn{d_2} is the number of observations and ordinal variables, respectively.}

\item{count}{\eqn{n \times d_3} matrix with the count reponse data, where \eqn{n} and \eqn{d_3} is the number of observations and count variables, respectively.}

\item{copnamesF1}{A vector with the names of possible candidates of bivariate copulas that link the each of the oberved variabels with the 1st factor. Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel, \dQuote{joe} for Joe, \dQuote{rjoe} for reflected Joe, \dQuote{1rjoe} for 1-reflected Joe, \dQuote{2rjoe} for 2-reflected Joe, \dQuote{BB1} for BB1, \dQuote{rBB1} for reflected BB1, \dQuote{BB7} for BB7, \dQuote{rBB7} for reflected BB7, \dQuote{BB8} for BB8, \dQuote{rBB8} for reflected BB8, \dQuote{BB10} for BB10, \dQuote{rBB10} for reflected BB10.}

\item{copnamesF2}{A list with the names of possible candidates of bivariate copulas that link the each of the oberved variabels with the 1st and 2nd factors. Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel, \dQuote{joe} for Joe, \dQuote{rjoe} for reflected Joe, \dQuote{1rjoe} for 1-reflected Joe, \dQuote{2rjoe} for 2-reflected Joe, \dQuote{BB1} for BB1, \dQuote{rBB1} for reflected BB1, \dQuote{BB7} for BB7, \dQuote{rBB7} for reflected BB7, \dQuote{BB8} for BB8, \dQuote{rBB8} for reflected BB8, \dQuote{BB10} for BB10, \dQuote{rBB10} for reflected BB10.}

\item{gl}{Gauss legendre quardrature nodes and weights.}
}


\value{A list containing the following components:
\item{``1st factor''}{The selected bivariate linking copulas for the 1st factor.}
\item{``2nd factor''}{The selected bivariate linking copulas for the 2nd factor.}
\item{AIC}{Akaike information criterion.}
\item{taus}{The estimated copula parameters in Kendall's tau scale.}
}


\details{The linking copulas at each factor are selected with a sequential algorithm under the initial assumption that linking copulas are Frank, and then sequentially copulas with non-tail quadrant independence are assigned to any of pairs where necessary to account for tail asymmetry (discrete data) or tail dependence (continuous data).
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
\donttest{
#------------------------------------------------
#                   Estimation
#------------------             -----------------
# Setting quadreture points
nq<-25  
gl<-gauss.quad.prob(nq) 
#------------------------------------------------
#                     PE Data
#------------------             -----------------
data(PE)
continuous.PE1 = -PE[,1]
continuous.PE <- cbind(continuous.PE1, PE[,2])
categorical.PE <- PE[, 3:5]

#------------------ One-factor  -----------------
# listing the possible copula candidates:
d <- ncol(PE)
copulasF1 <- rep(list(c("bvn", "bvt3", "bvt5", "frk", "gum", 
"rgum", "rjoe","joe", "1rjoe","2rjoe", "1rgum","2rgum")), d)
out1F.PE <- select1F(continuous.PE, categorical.PE, 
count=NULL, copulasF1, gl)

#------------------------------------------------
#------------------------------------------------
#                     GSS Data
#------------------             -----------------
data(GSS)
attach(GSS)
continuous.GSS <- cbind(INCOME, AGE)
ordinal.GSS <- cbind(DEGREE, PINCOME, PDEGREE) 
count.GSS <- cbind(CHILDREN, PCHILDREN)

#------------------ One-factor  -----------------
# listing the possible copula candidates:
d <- ncol(GSS)
copulasF1 <- rep(list(c("bvn", "bvt3", "bvt5", "frk", "gum",
"rgum", "rjoe","joe", "1rjoe","2rjoe", "1rgum","2rgum")), d)
out1F.GSS <- select1F(continuous.GSS, ordinal.GSS, count.GSS, copulasF1, gl)

#------------------ two-factor  -----------------
# listing the possible copula candidates:
copulasF1 = copulasF2 = rep(list(c("bvn", "bvt3", "bvt5", "frk",
"gum", "rgum", "rjoe","joe", "1rjoe","2rjoe", "1rgum","2rgum")), d)
out2F.GSS <- select2F(continuous.GSS, ordinal.GSS, 
count.GSS, copulasF1, copulasF2, gl)


}
}
\keyword{multivariate}



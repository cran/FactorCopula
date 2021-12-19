\name{M2.Factor}
\alias{M2.1F}
\alias{M2.2F}
\title{Goodness-of-fit of factor copula models for mixed data}
\description{ The limited information \eqn{M_2} statistic (Maydeu-Olivares and Joe, 2006)  of factor copula models for mixed continuous and discrete data.}
\usage{
M2.1F(tcontinuous, ordinal, tcount, cpar, copF1, gl)
M2.2F(tcontinuous, ordinal, tcount, cpar, copF1, copF2, gl, SpC)
}
\arguments{

\item{tcontinuous}{\eqn{n \times d_1} matrix with the transformed continuous to ordinal reponse data, where \eqn{n} and \eqn{d_1} is the number of observations and transformed continous  variables, respectively.}

\item{ordinal}{\eqn{n \times d_2} matrix with the ordinal reponse data, where \eqn{n} and \eqn{d_2} is the number of observations and ordinal variables, respectively.}

\item{tcount}{\eqn{n \times d_3} matrix with the transformed count to ordinal reponse data, where \eqn{n} and \eqn{d_3} is the number of observations and transformed count variables, respectively.}

\item{cpar}{A list of estimated copula parameters.}

\item{copF1}{\eqn{(d_1+d_2+d_3)}-vector with the names of bivariate copulas that link the each of the oberved variabels with the 1st factor. Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel, \dQuote{joe} for Joe, \dQuote{rjoe} for reflected Joe, \dQuote{1rjoe} for 1-reflected Joe, \dQuote{2rjoe} for 2-reflected Joe, \dQuote{BB1} for BB1, \dQuote{rBB1} for reflected BB1, \dQuote{BB7} for BB7, \dQuote{rBB7} for reflected BB7,\dQuote{BB8} for BB8, \dQuote{rBB8} for reflected BB8, \dQuote{BB10} for BB10, \dQuote{rBB10} for reflected BB10.}

\item{copF2}{\eqn{(d_1+d_2+d_3)}-vector with the names of bivariate copulas that link the each of the oberved variabels with the 2nd factor. Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel, \dQuote{joe} for Joe, \dQuote{rjoe} for reflected Joe, \dQuote{1rjoe} for 1-reflected Joe, \dQuote{2rjoe} for 2-reflected Joe, \dQuote{BB1} for BB1, \dQuote{rBB1} for reflected BB1, \dQuote{BB7} for BB7, \dQuote{rBB7} for reflected BB7,\dQuote{BB8} for BB8, \dQuote{rBB8} for reflected BB8, \dQuote{BB10} for BB10, \dQuote{rBB10} for reflected BB10.}

\item{gl}{Gauss legendre quardrature nodes and weights.}

\item{SpC}{Special case for the 2-factor copula model with BVN copulas. Select a bivariate copula at the 2nd factor to be fixed to independence. e.g. "SpC = 1" to set the first copula at the 2nd factor to independence.  }
}

\value{A list containing the following components:

\item{M2}{The \eqn{M_2} statistic which has a null asymptotic distribution that is \eqn{\chi^2} with \eqn{s-q} degrees of freedom, where \eqn{s} is the number of univariate and bivariate margins that do not include the category 0 and \eqn{q} is the number of model parameters.}

\item{df}{\eqn{s-q}.}

\item{p-value}{The resultant \eqn{p}-value.}

}



\details{The \eqn{M_2} statistic has been developed for goodness-of-fit testing in multidimensional contingency tables by Maydeu-Olivares and Joe (2006). 
Nikoloulopoulos and Joe (2015) have used the \eqn{M_2} statistic to assess the goodness-of-fit of factor copula models for ordinal data. We build on the aforementioned papers and propose a methodology to assess the overall goodness-of-fit of factor copula models for mixed continuous and discrete responses. Since the \eqn{M_2} statistic has been developed for multivariate ordinal data, we propose to first transform the continuous and count variables to ordinal and then calculate the \eqn{M_2} statistic at the maximum likelihood estimate before transformation. 
}

\references{
Kadhem, S.H. and Nikoloulopoulos, A.K. (2021)
Factor copula models for mixed data.
\emph{British Journal of Mathematical and Statistical Psychology}, \bold{74}, 365--403. \doi{10.1111/bmsp.12231}.

Maydeu-Olivares, A. and Joe, H. (2006).
Limited information goodness-of-fit testing in multidimensional contingency tables.
\emph{Psychometrika}, \bold{71}, 713--732. \doi{10.1007/s11336-005-1295-9}.

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
                            count=NULL, copF1=cop1f.PE, gl, hessian = T)
#------------------------------------------------
#                       M2  
#------------------------------------------------
#Transforming the continuous to ordinal data:
ncontinuous.PE = continuous2ordinal(continuous.PE, 5)
# M2 statistic for the one-factor copula model:

m2.1f.PE <- M2.1F(ncontinuous.PE, categorical.PE, tcount=NULL, 
       cpar=est1factor.PE$cpar, copF1=cop1f.PE, gl)
       
#------------------------------------------------
#                    GSS Data
#------------------             -----------------
data(GSS)
attach(GSS)
continuous.GSS <- cbind(INCOME,AGE)
ordinal.GSS <- cbind(DEGREE,PINCOME,PDEGREE) 
count.GSS <- cbind(CHILDREN,PCHILDREN)

#------------------------------------------------
#                   Estimation
#------------------             -----------------
# one-factor copula model
cop1f.GSS <- c("joe","2rjoe","bvt3","bvt3",
          "rgum","2rjoe","2rgum")
est1factor.GSS <- mle1factor(continuous.GSS, ordinal.GSS,
                        count.GSS, copF1=cop1f.GSS, gl, hessian = T)

#two-factor copula model
cop1.2f <- c("rgum","rjoe","bvn","1rjoe",
           "1rjoe","rjoe","gum")
cop2.2f <- c("gum","2rjoe","rjoe","gum",
           "bvt5","bvn","2rgum")
est2factor.GSS <- mle2factor(continuous.GSS, ordinal.GSS, 
        count.GSS, copF1=cop1.2f, copF2=cop2.2f, gl, hessian = T)

#------------------------------------------------
#                 Transformation
#------------------             -----------------
# Transforming the continuous to ordinal data:

# continuous1: Income
continuous1 = as.integer(cut(continuous.GSS[,1], 
c(0,10,19,29,40,100), include.lowest = T))
continuous1 = continuous1 - 1

# continuous2: AGE
continuous2 = as.integer(cut(continuous.GSS[,2] , 
c(0, 24, 44, 64, 100), include.lowest = T))
continuous2 = continuous2 - 1

# Combining the transformed continuous variables. 
ncontinuous.GSS <- cbind(continuous1, continuous2)

#---------------- COUNT VARIABLE -----------------
# count1 : CHILDREN
count1  = count.GSS[,1]
count1[count1 > 3] = 3

# count2: PCHILDREN
count2  = count.GSS[,2]
count2[count2 > 7] = 7

# Combining both transformed count variables
ncount.GSS = cbind(count1, count2)

#------------------------------------------------
#                       M2  
#------------------------------------------------
# M2 statistic for the one-factor copula model:

m2.1f.GSS <- M2.1F(ncontinuous.GSS, ordinal.GSS,
                      ncount.GSS, cpar = est1factor.GSS$cpar, 
                      copF1 = cop1f.GSS, gl)

#------------------------------------------------
# M2 statistic for the two-factor copula model:

m2.2f.GSS <- M2.2F(ncontinuous.GSS, ordinal.GSS, ncount.GSS, 
                      cpar = est2factor.GSS$cpar, copF1 = cop1.2f,     
                      copF2 = cop2.2f, gl)
                      
}
}
\keyword{multivariate}



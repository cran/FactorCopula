\name{semicorr}
\alias{semicorr}
\title{
Diagnostics to detect tail dependence or tail asymmetry.
}
\description{
The sample versions of the correlation \eqn{\rho_N}, upper semi-correlation \eqn{\rho_N^+} (correlation in the joint upper  quadrant) and lower semi-correlation \eqn{\rho_N^-} (correlation  in the joint lower quadrant). These  are sample linear  (when both variables are continuous),  polychoric (when both variables are ordinal),  and polyserial (when one variable is continuous and the
other is ordinal) correlations.  
}
\usage{
semicorr(dat,type)}
\arguments{
\item{dat}{Data frame of mixed continuous and ordinal data. }
\item{type}{A vector with 1's for the location of continuous data and 2's for the location of ordinal data.}
}

\value{
A matrix containing the following components for semicorr():

\item{rho}{\eqn{\rho_N}.}
\item{lrho}{\eqn{\rho_N^-}.}
\item{urho}{\eqn{\rho_N^+}.}
}


\references{


Joe, H. (2014).
\emph{Dependence Modelling with Copulas}. Chapman and Hall/CRC.

Kadhem, S.H. and Nikoloulopoulos, A.K. (2019)
Factor copula models for mixed data.
\emph{Arxiv e-prints}, <arXiv:1907.07395>. \url{https://arxiv.org/abs/1907.07395}.
}


\author{
Sayed H. Kadhem \email{s.kadhem@uea.ac.uk}\cr
Aristidis K. Nikoloulopoulos \email{a.nikoloulopoulos@uea.ac.uk}
}

\examples{

#------------------------------------------------
#                     PE Data
#------------------             -----------------
data(PE)
#correlation
continuous.PE1 <- -PE[,1] 
continuous.PE <- cbind(continuous.PE1, PE[,2])
categorical.PE <- data.frame(apply(PE[, 3:5], 2, factor))
nPE <- cbind(continuous.PE, categorical.PE)

#-------------------------------------------------
# Semi-correlations-------------------------------
#-------------------------------------------------
# Exclude the dichotomous variable
sem.PE = nPE[,-3]
semicorr.PE = semicorr(dat = sem.PE, type = c(1,1,2,2)) 
#------------------------------------------------
#------------------------------------------------
#                    GSS Data
#------------------             -----------------
data(GSS)
attach(GSS)
continuous.GSS <- cbind(INCOME,AGE)
ordinal.GSS <- cbind(DEGREE,PINCOME,PDEGREE) 
count.GSS <- cbind(CHILDREN,PCHILDREN)

# Transforming the COUNT variables to ordinal
# count1 : CHILDREN
count1  = count.GSS[,1]
count1[count1 > 3] = 3

# count2: PCHILDREN
count2  = count.GSS[,2]
count2[count2 > 7] = 7

# Combining both transformed count variables
ncount.GSS = cbind(count1, count2)

# Combining ordinal and transformed count variables
categorical.GSS <- cbind(ordinal.GSS, ncount.GSS)
categorical.GSS <- data.frame(apply(categorical.GSS, 2, factor))

# combining continuous and categorical variables
nGSS = cbind(continuous.GSS, categorical.GSS)
#-------------------------------------------------
# Semi-correlations-------------------------------
#-------------------------------------------------
semicorr.GSS = semicorr(dat = nGSS, type = c(1, 1, rep(2,5)))




}

\keyword{correlation}

\name{discrepancy}
\alias{discrepancy}
\title{
Diagnostics to detect a factor dependence structure 
}
\description{
The diagnostic method in Joe (2014, pages 245-246)  to show that each dataset has  a factor  structure based on linear factor analysis. The  correlation matrix \eqn{\R_{\mathrm{observed}}}  has been obtained based on the sample correlations  from the bivariate pairs of the observed variables. 
These are the linear  (when both variables are continuous),  polychoric (when both variables are ordinal),  and polyserial (when one variable is continuous and the other is ordinal) sample correlations among the observed variables. 
The resulting \eqn{\R_{\mathrm{observed}}}  is generally positive definite  if the sample size is not small enough; if not one has to convert it to positive definite. We calculate various measures of discrepancy between \eqn{\R_{\mathrm{observed}}}  and \eqn{\R_{\mathrm{model}}} (the resulting correlation matrix of linear factor analysis), such as  the maximum absolute correlation
difference \eqn{D_1=\max|\R_{\mathrm{model}} - \R_{\mathrm{observed}}|}, the average absolute correlation
difference \eqn{D_2=\mathrm{avg}| \R_{\mathrm{model}} - \R_{\mathrm{observed}}|}, and  the correlation matrix discrepancy measure \eqn{D_3=\log\bigl( \det(\R_{\mathrm{model}}) \bigr) - \log\bigl( \det(\R_{\mathrm{observed}})\bigr) + \mathrm{tr}( \R^{-1}_{\mathrm{model}} \R_{\mathrm{observed}} ) - d}.
}
\usage{
discrepancy(cormat, n, f3)}
\arguments{
\item{cormat}{\eqn{\R_{\mathrm{observed}}}.}
\item{n}{Sample size.}
\item{f3}{If TRUE, then the linear 3-factor analysis is fitted.}
}

\value{
A matrix with the calculated discrepancy measures for different number of factors. 
}


\references{


Joe, H. (2014).
\emph{Dependence Modelling with Copulas}. Chapman & Hall, London.

Kadhem, S.H. and Nikoloulopoulos, A.K. (2021)
Factor copula models for mixed data.
\emph{British Journal of Mathematical and Statistical Psychology}, \bold{74}, 365--403. \doi{10.1111/bmsp.12231}.
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
u.PE <- apply(continuous.PE, 2, rank)/(nrow(PE)+1)
z.PE <- qnorm(u.PE)
categorical.PE <- data.frame(apply(PE[, 3:5], 2, factor))
nPE <- cbind(z.PE, categorical.PE)

#-------------------------------------------------
# Discrepancy measures----------------------------
#-------------------------------------------------
#correlation matrix for mixed data
cormat.PE <- as.matrix(polycor::hetcor(nPE, std.err=FALSE))
#discrepancy measures
out.PE = discrepancy(cormat.PE, n = nrow(nPE), f3 = FALSE)

#------------------------------------------------
#------------------------------------------------
#                    GSS Data
#------------------             -----------------
data(GSS)
attach(GSS)
continuous.GSS <- cbind(INCOME,AGE)
continuous.GSS <- apply(continuous.GSS, 2, rank)/(nrow(GSS)+1)
z.GSS <- qnorm(continuous.GSS)
ordinal.GSS <- cbind(DEGREE,PINCOME,PDEGREE) 
count.GSS <- cbind(CHILDREN,PCHILDREN)

# Transforming the count variables to ordinal
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
nGSS = cbind(z.GSS, categorical.GSS)

#-------------------------------------------------
# Discrepancy measures----------------------------
#-------------------------------------------------
#correlation matrix for mixed data
cormat.GSS <- as.matrix(polycor::hetcor(nGSS, std.err=FALSE))
#discrepancy measures
out.GSS = discrepancy(cormat.GSS, n = nrow(nGSS), f3 = TRUE)

}


\keyword{correlation}

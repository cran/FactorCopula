\name{polychoric0}
\alias{polychoric0}
\Rdversion{1.1}
\title{
Polychoric correlation
}
\description{
Polychoric correlation
}
\usage{
polychoric0(odat,iprint=FALSE,prlevel=0)    
}
\arguments{
\item{odat}{nxd matrix of ordinal responses in 0,...,(ncat-1) or 1,...,ncat}
\item{iprint}{flag for printing of intermediate results, including
bivariate tables for observed versus expected assuming discretized
bivariate Gaussian}
\item{prlevel}{print.level for nlm for numerical optimization}
}
\details{
Polychoric correlation for ordinal random variables. The number of categories can vary. 
}
\value{
 $polych is
a polychoric correlation matrix  based on two-stage estimate;
$iposdef is an indicator if the 2-stage correlation matrix estimate
is positive definite. 

}
\examples{
data(PTSD)
ydat=PTSD
rmat=polychoric0(ydat)$p
}

\keyword{dependence}
\keyword{maximum likelihood}


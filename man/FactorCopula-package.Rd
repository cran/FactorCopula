\name{FactorCopula-package}
\alias{FactorCopula-package}
\docType{package}
\title{Factor copula models for mixed continuous and disrete data}
\description{
Estimation, model selection and goodness-of-fit of factor copula models for mixed continuous and discrete data in Kadhem and Nikoloulopoulos (2019) <arXiv:1907.07395>.
}
\details{
This package contains \R functions for: 
\itemize{
\item diagnostics based on semi-correlations (Kadhem and Nikoloulopoulos, 2019; Joe, 2014) to detect tail dependence or tail asymmetry; 
\item diagnostics  to show that a  dataset  has a factor structure based on linear factor analysis (Kadhem and Nikoloulopoulos, 2019; Joe, 2014);
\item estimation of the factor copula models in Krupskii and Joe (2013), Nikoloulopoulos and Joe (2015), and Kadhem and Nikoloulopoulos (2019);

\item model selection of the factor copula models in Krupskii and Joe (2013), Nikoloulopoulos and Joe (2015), and Kadhem and Nikoloulopoulos (2019) using the heuristic algorithm in Kadhem and Nikoloulopoulos (2019) that automatically selects the bivariate parametric copula families that link the observed to the latent variables;

\item goodness-of-fit of the factor copula models in Krupskii and Joe (2013), Nikoloulopoulos and Joe (2015), and Kadhem and Nikoloulopoulos (2019) using the \eqn{M_2} statistic (Maydeu-Olivares and Joe, 2006). Note that the continuous  and count data have to be transformed to ordinal. 


}



}
\author{
Sayed H. Kadhem \email{s.kadhem@uea.ac.uk}\cr
Aristidis K. Nikoloulopoulos \email{a.nikoloulopoulos@uea.ac.uk}
}
\references{


Joe, H. (2014).
\emph{Dependence Modelling with Copulas}. Chapman and Hall/CRC.

Maydeu-Olivares, A. and Joe, H. (2006).
Limited information goodness-of-fit testing in multidimensional contingency tables.
\emph{Psychometrika}, \bold{71}, 713--732. \doi{10.1007/s11336-005-1295-9}.

Kadhem, S.H. and Nikoloulopoulos, A.K. (2019)
Factor copula models for mixed data.
\emph{Arxiv e-prints}, <arXiv:1907.07395>. \url{https://arxiv.org/abs/1907.07395}.

Krupskii, P. and Joe, H. (2013)
Factor copula models for multivariate data.
\emph{Journal of Multivariate Analysis}, \bold{120}, 85--101. \doi{10.1016/j.jmva.2013.05.001}.


Nikoloulopoulos, A.K. and Joe, H. (2015)
Factor copula models with item response data.
\emph{Psychometrika}, \bold{80}, 126--150. \doi{10.1007/s11336-013-9387-4}.


}
\keyword{package}

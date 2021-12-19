\name{FactorCopula-package}
\alias{FactorCopula-package}
\docType{package}
\title{Factor, bi-factor and second-order copula models}
\description{
Estimation, model selection and goodness-of-fit of (1) factor copula models for mixed continuous and discrete data in Kadhem and Nikoloulopoulos (2021a); (2) bi-factor and second-order copula models for item response data in Kadhem and Nikoloulopoulos (2021b).
}
\details{
This package contains \R functions for: 
\itemize{
\item diagnostics based on semi-correlations (Kadhem and Nikoloulopoulos, 2021a,b; Joe, 2014) to detect tail dependence or tail asymmetry; 
\item diagnostics  to show that a  dataset  has a factor structure based on linear factor analysis (Kadhem and Nikoloulopoulos, 2021a,b ; Joe, 2014);
\item estimation of the factor copula models in Krupskii and Joe (2013), Nikoloulopoulos and Joe (2015), and Kadhem and Nikoloulopoulos (2021a, 2021b);

\item model selection of the factor copula models in Krupskii and Joe (2013), Nikoloulopoulos and Joe (2015) and Kadhem and Nikoloulopoulos (2021a, 2021b) using the heuristic algorithms in Kadhem and Nikoloulopoulos (2021a, 2021b) that automatically selects the bivariate parametric copula families that link the observed to the latent variables;

\item goodness-of-fit of the factor copula models in Krupskii and Joe (2013), Nikoloulopoulos and Joe (2015) and Kadhem and Nikoloulopoulos (2021a, 2021b) using the \eqn{M_2} statistic (Maydeu-Olivares and Joe, 2006). Note that the continuous  and count data have to be transformed to ordinal. 
}



}
\author{
Sayed H. Kadhem \email{s.kadhem@uea.ac.uk}\cr
Aristidis K. Nikoloulopoulos \email{a.nikoloulopoulos@uea.ac.uk}
}
\references{


Joe, H. (2014).
\emph{Dependence Modelling with Copulas}. Chapman & Hall, London.

Maydeu-Olivares, A. and Joe, H. (2006).
Limited information goodness-of-fit testing in multidimensional contingency tables.
\emph{Psychometrika}, \bold{71}, 713--732. \doi{10.1007/s11336-005-1295-9}.


Kadhem, S.H. and Nikoloulopoulos, A.K. (2021a)
Factor copula models for mixed data.
\emph{British Journal of Mathematical and Statistical Psychology}, \bold{74}, 365--403. \doi{10.1111/bmsp.12231}.


Kadhem, S.H. and Nikoloulopoulos, A.K. (2021b)
Bi-factor and second-order copula models for item response data.
\emph{Arxiv e-prints}, <arXiv:2102.10660>. \url{https://arxiv.org/abs/2102.10660}.

Krupskii, P. and Joe, H. (2013)
Factor copula models for multivariate data.
\emph{Journal of Multivariate Analysis}, \bold{120}, 85--101. \doi{10.1016/j.jmva.2013.05.001}.


Nikoloulopoulos, A.K. and Joe, H. (2015)
Factor copula models with item response data.
\emph{Psychometrika}, \bold{80}, 126--150. \doi{10.1007/s11336-013-9387-4}.


}
\keyword{package}

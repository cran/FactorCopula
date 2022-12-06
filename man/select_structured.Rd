\name{Select.StructuredFactor}
\alias{Select.StructuredFactor}
\alias{selectBifactor}
\alias{selectSecond_order}
\title{Model selection of the bi-factor and second-order copula models for item response data}
\description{A heuristic algorithm that automatically selects the bivariate parametric copula families for the bi-factor and second-order copula models for item response data.}
\usage{
selectBifactor(y, grpsize, copnames, gl)
selectSecond_order(y, grpsize, copnames, gl)
}
\arguments{
\item{y}{\eqn{n \times d} matrix with the item reponse data, where \eqn{n} and \eqn{d} is the number of observations and variables, respectively.}

\item{grpsize}{vector indicating the size for each group, e.g., c(4,4,4) indicating four items in all three groups.}

\item{copnames}{A vector with the names of possible candidates of bivariate copulas to be selected for the bi-factor and second-order copula models for item response data. Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel.}

\item{gl}{Gauss legendre quardrature nodes and weights.}
}


\value{A list containing the following components:
\item{``common factor''}{The selected bivariate linking copulas for the common factor (Bi-factor: copulas link items with the common factor. Second-order: copulas link group-specific factors with the common factor).}
\item{``group-specific factor''}{The selected bivariate linking copulas for the items with group-speicifc factors.}
\item{log-likelihood}{The maximized joint log-likelihood.}
\item{taus}{The estimated copula parameters in Kendall's tau scale.}
}


\details{The linking copulas at each factor are selected with a sequential algorithm under the initial assumption that linking copulas are BVN, and then sequentially copulas with tail dependence are assigned to any of pairs where necessary to account for tail asymmetry.
}


\references{
Kadhem, S.H. and Nikoloulopoulos, A.K. (2022)
Bi-factor and second-order copula models for item response data.
\emph{Psychometrika}. \doi{doi:10.1007/s11336-022-09894-2}.
}
\author{
Sayed H. Kadhem \email{s.kadhem@uea.ac.uk}\cr
Aristidis K. Nikoloulopoulos \email{a.nikoloulopoulos@uea.ac.uk}
}


\examples{
\donttest{
#------------------------------------------------
# Setting quadreture points
nq <- 15
gl <- gauss.quad.prob(nq)
#------------------------------------------------
#                     TAS Data
#------------------             -----------------
data(TAS)
#using a subset of the data
#group1: 1,3,6,7,9,13,14
grp1=c(1,3)
#group2: 2,4,11,12,17
grp2=c(2,4)
#group3: 5,8,10,15,16,18,19,20
grp3=c(5,8)
#Rearrange items within testlets
set.seed(123)
i=sample(1:nrow(TAS),500)
ydat=TAS[i,c(grp1,grp2,grp3)]

#size of each group
g1=length(grp1)
g2=length(grp2)
g3=length(grp3)
grpsize=c(g1,g2,g3)

# listing possible copula candidates:
copnames=c("bvt2","bvt3")

Bifactor_model = selectBifactor(ydat, grpsize, copnames, gl)
}
}
\keyword{multivariate}



\name{Select.FactorTree}
\alias{Select.FactorTree}
\alias{selectFactorTrVine}
\alias{copselect1FactorTree}
\alias{copselect2FactorTree}
\title{Model selection of the 1- and 2-factor tree copula models for item response data}
\description{Heuristic algorithms that automatically select the bivariate parametric copula families and 1-truncated vine tree structure for the 1- and 2-factor tree copula models for item response data.}
\usage{
selectFactorTrVine(y,rmat,alg)
copselect1FactorTree(y, A, copnames, gl)
copselect2FactorTree(y, A, copnames, gl)
}
\arguments{
\item{y}{\eqn{n \times d} matrix with the item reponse data, where \eqn{n} and \eqn{d} is the number of observations and variables, respectively.}

\item{rmat}{Polychoric correlation matrix of y.}

\item{alg}{1-truncated vine selection using partial and polychoric correlations. Choose a number from (1,2,3). 
\bold{1}: Factor correlation selection algorithm for 1-factor tree copula.
\bold{2}: Factor correlation selection algorithm for 2-factor tree copula.
\bold{3}: Polychoric correlation selection algorithm for 1-factor/2-factor tree copulas, here any other correlation matrix (rmat) can be used.}

\item{A}{\eqn{d \times d} vine array with \eqn{1,...,d} on diagonal, note only the first row and diagnoal values are used for the 1-truncated vine model} 

\item{copnames}{A vector with the names of possible candidates of bivariate copulas to be selected for the 1-factor/2-factor copula models for item response data. Choices are \dQuote{bvn} for BVN, \dQuote{bvt\eqn{\nu}} with \eqn{\nu = \{1, \ldots, 9\}} degrees of freedom for t-copula, \dQuote{frk} for Frank, \dQuote{gum} for Gumbel, \dQuote{rgum} for reflected Gumbel, \dQuote{1rgum} for 1-reflected Gumbel, \dQuote{2rgum} for 2-reflected Gumbel.}

\item{gl}{Gauss legendre quardrature nodes and weights.}
}


\value{A list containing the following components:
\item{``1Factor copulas''}{The selected bivariate linking copulas for the first factor.}
\item{``2Factor copulas''}{The selected bivariate linking copulas for the second factor.}
\item{``Vine tree copulas''}{The selected bivariate linking copulas for the 1-truncated vine tree.}
\item{``log-likelihood''}{The maximized joint log-likelihood.}
\item{``estimated taus''}{The estimated copula parameters in Kendall's tau scale.}
\item{F1treeA}{Vine array for the 1-factor tree copula.}
\item{F2treeA}{Vine array for the 2-factor tree copula.}
\item{pmat_f1}{Factor correlation matrix for the 1-factor tree copula; \eqn{(Y_{j1}Y_{j2}|V_1)}}
\item{pmat_f2}{Factor correlation matrix for the 2-factor tree copula;
\eqn{(Y_{j1}Y_{j2}|V_1V_2)}}
}

\details{
The model selection involves two separate steps. At the first step we select the 1-truncated vine or Markov tree structure assuming our models are constructed with BVN copulas. We find the maximum spanning tree which is a tree on all nodes that maximizes the pairwise dependencies using the well-known algorithm of Prim (1957).   That is we find the tree with \eqn{d-1} edges \eqn{\mathcal{E}} that minimizes 
\eqn{\sum_{\mathcal{E}}\log(1-r_{jk}^2).}
The minimum spanning tree algorithm of Prim (1957)
guarantees to find the optimal solution when edge weights between nodes  are given by \eqn{\log(1-r_{jk}^2)}. 
We use two  different measures of pairwise dependence \eqn{r_{jk}}. The first measure is the estimated  polychoric correlation (Olsson, 1979). The second measures of pairwise dependence    is  the correlation between each of the items and the first factor along with the partial correlation  between each of the items and the second  factor  given the first factor (Nikoloulopoulos and Joe, 2015). We call polychoric and factor correlation selection algorithm when the pairwise dependence is the polychoric and factor correlation, respectively.


At the the second step, we sequentially select suitable bivariate copulas to account for any tail dependence/asymmetry. We start from the 1st tree that includes copulas connecting the observed variables to the first factor and then we iterate over a set of copula candidates and select suitable bivariate copulas in the subsequent trees. We select only one copula family for all the edges in each tree to ease interpretation.
}


\references{



Fontenla, M. (2014).
\emph{optrees: Optimal Trees in Weighted Graphs}.
R package version 1.0.


Kadhem, S.H. and Nikoloulopoulos, A.K. (2022b)
Factor tree copula models for item response data.
\emph{Arxiv e-prints}, <arXiv: 2201.00339>. \url{https://arxiv.org/abs/2201.00339}.

Nikoloulopoulos, A.K. and Joe, H. (2015)
Factor copula models with item response data.
\emph{Psychometrika}, \bold{80}, 126--150. \doi{10.1007/s11336-013-9387-4}.

Olsson, F. (1979) Maximum likelihood estimation of the polychoric correlation coefficient. \emph{Psychometrika}, \bold{44}, 443--460.

Prim, R. C. (1957) Shortest connection networks and some generalizations. \emph{The Bell System Technical Journal}, \bold{36}, 1389--1401.

}
\author{
Sayed H. Kadhem \cr
Aristidis K. Nikoloulopoulos \email{a.nikoloulopoulos@uea.ac.uk}
}


\examples{
\donttest{
#------------------------------------------------
# Setting quadreture points
nq <- 25
gl <- gauss.quad.prob(nq)
#------------------------------------------------
#                    PTSD Data
#------------------             -----------------
data(PTSD)
ydat=PTSD

#------------------------------------------------
#vine tree structure
#selecting vine tree based on polychoric corr
rmat=polychoric0(ydat)$p
A.polychoric=selectFactorTrVine(y=ydat,rmat,alg=3)

#selecting bivariate copulas for the 1-2-factor tree
# listing the possible copula candidates:
copnames=c("bvn", "bvt2", "bvt3",
"gum", "rgum")

f1tree_model = copselect1FactorTree(ydat, A.polychoric$VineTreeA, 
copnames, gl)
f2tree_model = copselect2FactorTree(ydat, A.polychoric$VineTreeA, 
copnames, gl)

}
}
\keyword{multivariate}



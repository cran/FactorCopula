\name{transformation}
\alias{continuous2ordinal}
\alias{count2ordinal}
\title{
Continuous/count to ordinal responses
}
\description{
Transforming a continuous/count to ordinal variable with \eqn{K} categories. 
}
\usage{
continuous2ordinal(continuous, categ)
count2ordinal(count, categ)}
\arguments{
\item{continuous}{Matrix of continuous data.}
\item{count}{Matrix of count data.}
\item{categ}{The number of categories \eqn{K}.}
}

\examples{

#------------------------------------------------
#                     PE Data
#------------------             -----------------
data(PE)
continuous.PE <- PE[, 1:2]

#Transforming the continuous to ordinal data :
tcontinuous = continuous2ordinal(continuous.PE, categ=5)
table(tcontinuous)

#Transforming the count to ordinal data:
set.seed(12345)
count.data = rpois(1000, 3)
tcount = count2ordinal(count.data, 5)
table(tcount)


}

\keyword{univar}

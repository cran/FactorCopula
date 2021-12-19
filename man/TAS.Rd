\name{TAS}
\alias{TAS}
\docType{data}
\title{
Toronto Alexithymia Scale (TAS)
}
\description{
The Toronto Alexithymia Scale is the most utilized measure of alexithymia in empirical research and is composed of \eqn{d = 20} items that can be subdivided into \eqn{G = 3} non-overlapping groups: \eqn{d_1 = 7} items to assess difficulty identifying feelings (DIF), \eqn{d_2 = 5} items to assess difficulty describing feelings (DDF) and \eqn{d_3 = 8} items to assess externally oriented thinking (EOT).
Students were 17 to 25 years old and 58\% of them were female and 42\% were male. They were asked to respond to each item using one of K = 5 categories: ``1 = completely disagree'', ``2 = disagree'', ``3 = neutral'', ``4 = agree'', ``5 = completely agree''.
}
\usage{data(TAS)}
\format{
  A data frame with 1925 observations on the following 20 items:
  \describe{
    \item{\code{DIF}}{items: 1,3,6,7,9,13,14.}
    \item{\code{DDF}}{items: 2,4,11,12,17.}
    \item{\code{EOT}}{items: 5,8,10,15,16,18,19,20.}
  }
}

\source{
Briganti, G. and Linkowski, P. (2020).
Network approach to items and domains from the toronto alexithymia scale.
\emph{Psychological Reports}, \bold{123}, 2038--2052.

Williams, D. and Mulder, J. (2020).
\emph{BGGM: Bayesian Gaussian Graphical Models}.
R package version 1.0.0.
}

\keyword{datasets}

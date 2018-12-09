\name{ars}
\alias{ars}
\title{
Adaptive-Rejection Sampler
}
\description{
An implementation of the adaptive-rejection sampling algorithm
in Gilks, et. al. (1992). It allows one to perform adaptive-rejection sampling
on a log-concave density.
}
\usage{
ars(FUN, n = 1, D = c(-Inf, Inf))
}
\arguments{
  \item{FUN}{log-concave density function. Need not be normalized.}
  \item{n}{optional; number of points to sample}
  \item{D}{optional; domain of the density. If unspecified, domain will be
  -Inf to Inf}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Vector of length n containing samples from specified input density.
}
\references{
Gilks, W. R., and P. Wild. “Adaptive Rejection Sampling for Gibbs Sampling.” Journal of the Royal Statistical Society. Series C (Applied Statistics), vol. 41, no. 2, 1992, pp. 337–348. JSTOR, JSTOR, www.jstor.org/stable/2347565.
}
\author{
Vaibhav Ramamoorthy, Jonathan Lee, Colin Kou, Brandon Mannion
}
\note{
%%  ~~further notes~~
For a density with extreme domain values, it is suggested first to transform the
density to get a smaller mean before feeding it into the algorithm and then to scale it back after the sample is generated.
}
\examples{
sample <- ars(dnorm, n = 10000, D = c(-Inf, Inf))
sample <- ars(dnorm, n = 10000, D = c(0, Inf)) # truncated normal
sample <- ars(function(t) dbeta(t, 1, 1), n = 10000, D = c(0, 1))
}

\name{prob.B0}
\alias{prob.B0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Limiting Distribution Induced by the Bartlett Mutation Model
}
\description{
This discrete distribution, as a limiting distribution induced by the Bartlett model (Bartlett 1955; Zheng 2008),
was first discovered by Zheng (2010). This distribution is often denoted by B0(A,k), where A and k are both 
positive real parameters. It was later found that this distribution can also be regarded as a
continuous mixture of the Luria-Delbruck distribution, where the mixing distribution of the parameter m is a
gamma distribution (Zheng 2011). Therefore, the relation between the Luria-Delbruck distribution and the B0
distribution is similar to that between the Poisson distribution and the negative binomial distribution.
}
\usage{
prob.B0(A, k, n = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
A positive real number.
}
  \item{k}{
A positive real number.
}
  \item{n}{
An non-negative integer specifying the number of probabilities to be computed.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector containing the first n+1 probabilities.
}
\references{

M.S. Bartlett, An Introduction to Stochastic Processes, Cambridge
University Press, London, 1955 (pp.132-135).

Q. Zheng, On Bartlett's formulation of the Luria-Delbruck
mutation model, Mathematical Biosciences 215 (2008) 48-54.


Q. Zheng, A new discrete distribution induced by the
Luria-Delbruck mutation model, Statistics 44 (2010) 529-540.


Q. Zheng, A Bayesian two-level model for fluctuation assay,
Genetica 139 (2011) 1409-1416.

}
\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
prob.LD, prob.MK
}
\examples{
prob.B0(5,0.8,12) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

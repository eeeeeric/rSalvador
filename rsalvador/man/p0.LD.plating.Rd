\name{p0.LD.plating}
\alias{p0.LD.plating}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
P0 method modified for partial plating
}
\description{
The classic P0 method proposed by Luria and Delbruck (1943) does not allow for partial plating. A modified P0
method that accounts for partial plating was known to Stewart (1991) and to Jones (1993). Algorithms for
computing maximum likelihood estimates of m was later proposed by Zheng (2008).
}
\usage{
p0.LD.plating(data, e)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data vector of mutant counts
}
  \item{e}{
A known plating efficiency
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An estimate of m, the expected number of mutations per culture
}
\references{
M.E. Jones, Accounting for plating efficiency when estimating spontaneous mutation rates, Mutation
Research 292 (1993) 187-189.

S.E. Luria, M. Delbruck, Mutations of bacteria
from virus sensitivity to virus resistance, Genetics 28 (1943) 491-511.

F.M. Stewart, Fluctuation analysis: the effect of plating efficiency, Genetica 54 (1991) 51-55.

Q. Zheng, A note on plating efficiency in fluctuation experiments,
Mathematical Biosciences 216 (2008) 150-153.

}
\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
newton.LD.plating
}
\examples{
p0.LD.plating(luria.16.data,e=0.4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

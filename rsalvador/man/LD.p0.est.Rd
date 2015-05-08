\name{LD.p0.est}
\alias{LD.p0.est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Classic P0 Method of Luria and Delbruck
}
\description{
This is one of the two methods proposed in the classic paper of Luria and Delbruck (1943)
to estimate the expected number of mutations per test tube, now commonly known as m. An 
important property of this method is often overlooked: the P0 method does not rely on the
assumption that mutant cells and wild type cells grow at the same rate.
}
\usage{
LD.p0.est(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a vector of non-negative integers representing mutant cell counts.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An estimate of m, the expected number of mutations.
}
\references{
S.E. Luria, M. Delbruck, Mutations of bacteria
from virus sensitivity to virus resistance, Genetics 28 (1943) 491-511.
}
\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
newton.LD
}
\examples{
LD.p0.est(cairns.foster.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

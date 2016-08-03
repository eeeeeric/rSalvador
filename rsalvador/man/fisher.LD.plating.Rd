\name{fisher.LD.plating}
\alias{fisher.LD.plating}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computing Expected Fisher Information for m Under the Standard Luria-Delbruck Distribution with Partial Plating
}
\description{
This function computes approximate expected Fisher information for m under the standard Luria-Delbruck
distribution with partial plating. The algorithm is adapted from that given in Zheng (2002). The expected Fisher
information can be utilized in determining sample size (the number of cultures) in a fluctuation experiment.
}
\usage{
fisher.LD.plating(m, e, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
The excepted number of mutations.
}
  \item{e}{
The parameter of plating efficiency (0<e<1).
}
  \item{n}{
The number of terms kept in summing an infinite series that defines the expected Fisher information.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A positive real number.
}
\references{
Q. Zheng, Statistical and algorithmic methods for fluctuation analysis with SALVADOR as an
implementation, Mathematical Biosciences 176 (2002) 237-252.
}
\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
fisher.MK
}
\examples{fisher.LD.plating(m=50.0,e=0.2,n=2000)}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

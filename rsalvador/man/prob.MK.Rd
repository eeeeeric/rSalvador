\name{prob.MK}
\alias{prob.MK}

\title{
Calculating the Mandelbrot-Koch Probability Distribution
}
\description{
This mutant distribution was first proposed by Mandelbrot (1974). Kock (1982) and
Stewart et al. (1990) made independent contributions. Explicit algorithm for computing
the probabilities were given by Zheng (2005). 
}
\usage{
prob.MK(m, w = 1, n)
}

\arguments{
  \item{m}{Expected number of mutation per culture.

}
  \item{w}{Relative fitness w=(birth rate of mutants)/(birth rate of nonmutants).

}
  \item{n}{The first n+1 probabilities will be computed.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A vector of n+1 probabilities.
}
\references{

A. L. Koch, Mutation and growth rates from Luria-Delbruck fluctuation 
tests, Mutation Research 95 (1982) 129-143.

B. Mandelbrot, A population birth-and-mutation 
process, I: Explicit distributions for the number of mutants in an old culture
of bacteria, J Appl Prob 11 (1974) 437-444.

F. M. Stewart, D. M. Gordon, B. R.  Levin, Fluctuation
analysis: the probability distribution of the number of mutants under
different conditions, Genetics 124 (1990) 175-185.

Q. Zheng, New algorithms for Luria-Delbruck fluctuation analysis,
Mathematical Biosciences 196 (2005) 198-214.

}
\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
prob.LD, prob.Haldane
}
\examples{prob.MK(2.3,1.5,20)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

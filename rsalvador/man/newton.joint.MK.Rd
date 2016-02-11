\name{newton.joint.MK}
\alias{newton.joint.MK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Joint Estimation of m and w Under the Mandelbrot-Koch Model
}
\description{
This function computes maximum likelihood estimates of m and w simultaneously under the
Mandelbrot-Koch model (Mandelbrot 1974; Koch 1982). Here m is the expected number of mutations
per culture, and w is relative fitness of mutants. The Mandelbrot-Koch model allows for
differential growth between mutants and nonmutants. The Newton-Raphson algorithm adopted
is described in detail in Zheng (2005).
}
\usage{
newton.joint.MK(data, tol = 1e-08, init.m = -1, init.w = -1,
   max.iter = 30, show.iter = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A vector containing mutant counts.
}
  \item{tol}{
Tolerance parameter to control the number of numerical iterations.
}
  \item{init.m}{
A starting value for the m parameter, the expected number of mutations per culture.
}
  \item{init.w}{
A starting value for the w parameter, relative fitness of the mutants.
}
  \item{max.iter}{
Maximum number of Newton-Raphson iterations.
}
  \item{show.iter}{
A logical variable enabling the display of the numerical iteration process.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A pair of non-negative numbers representing the maximum likelihood estimates of m and w.
}
\references{
A. L. Koch, Mutation and growth rates from Luria-Delbruck fluctuation 
tests, Mutation Research 95 (1982) 129-143.

B. Mandelbrot, A population birth-and-mutation 
process, I: Explicit distributions for the number of mutants in an old culture
of bacteria, J Appl Prob 11 (1974) 437-444.

Q. Zheng, New algorithms for Luria-Delbruck fluctuation analysis,
Mathematical Biosciences 196 (2005) 198-214.
}
\author{
Qi Zheng <qzheng@sph.tamhse.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
newton.MK
}
\examples{
newton.joint.MK(demerec.data)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

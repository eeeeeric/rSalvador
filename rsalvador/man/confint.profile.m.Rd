\name{confint.profile.m}
\alias{confint.profile.m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Profile Likelihood Confidence Intervals for m Under the Mandelbrot-Koch Model
}
\description{
This function computes profile likelihood confidence intervals for m under the
Mandelbrot-Koch model (Mandelbrot 1974; Koch 1982). Here m is the expected number of mutations
per culture. The Mandelbrot-Koch model allows for differential growth between mutants and
nonmutants, that is, the relative fitness of mutants, w,  can be different from 1.0. The
adopted iterative algorithm is described in detail in Zheng (2005).
}
\usage{
confint.profile.m(data, alpha = 0.05, init.low.m = -1, init.low.w = -1,
    init.up.m = -1, init.up.w = -1, init.m = -1, init.w = -1, 
       max.iter = 30, tol = 1e-09, show.iter = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A vector containing mutant counts.
}
  \item{alpha}{
The confidence level is 1.0-alpha.
}
  \item{init.low.m}{
A starting value for the lower confidence limit of m.
}
  \item{init.low.w}{
A starting value for the lower confidence limit of w.
}
  \item{init.up.m}{
A starting value for the upper confidence limit of m.
}
  \item{init.up.w}{
A starting value for the upper confidence limit of w.
}
  \item{init.m}{
A starting value for the MLE of m, needed for computing the profile likelihood CI.
}
  \item{init.w}{
A starting value for the MLE of w, needed for computing the profile likelihood CI.
}
  \item{max.iter}{
Maximum number of numerical iterations.
}
  \item{tol}{
A tolerance parameter to control the number of iterations.
}
  \item{show.iter}{
A logical variable controlling the display of numerical iterations.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A pair of real numbers representing the profile confidence interval for m.
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
Qi Zheng <qzheng.@sph.tamhsc.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
confint.MK, confint.profile.w
}
\examples{
confint.profile.m(demerec.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

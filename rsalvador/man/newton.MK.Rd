\name{newton.MK}
\alias{newton.MK}

\title{
Computing MLE of m Under the Mandelbrot-Koch Model
}
\description{
This function computes maximum likelihood estimate of m, the expected number of mutations
per culture, under the Mandelbrot-Koch model that allows for differential growth between mutants
and nonmutants (Mandelbrot 1974; Koch 1982). The relative fitness (w) is usually obtained from
a fitness assay (aka competition assay). The algorithm is a simple specialization of a more
general algorithm given in Zheng (2005).
}
\usage{
newton.MK(data, w = 1, tol = 1e-08, init.m = 0, max.iter = 30,
          show.iter = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data vector containing the numbers of mutant cells from a fluctuation experiment.
}
  \item{w}{
Relative fitness define as (growth rate of mutants)/(growth rate of nonmutants).
}
  \item{tol}{
Tolerance parameter to set a desired convergence criterion.
}
  \item{init.m}{
An initial guess of m, it is rarely needed.
}
  \item{max.iter}{
A parameter to set the maximum number of iterations for the Newton-Raphson procedure.
}
  \item{show.iter}{
A logical variable to show the numerical iteration process.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A positive number as an MLE of m.
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
newton.MK(demerec.data, w = 1.25)
}

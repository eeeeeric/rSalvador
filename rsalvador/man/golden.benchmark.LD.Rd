\name{golden.benchmark.LD}
\alias{golden.benchmark.LD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum Likelihood Estimation of the Mutation Rate Accounting for Variation in Nt
}
\description{
Nt is the (non-mutant) cell population size immediately prior to plating. In a fluctuation experiment 
comprising n tubes, the final cell population sizes are Nt(1), Nt(2), ..., Nt(n). Measuring all Nt(i)
is a daunting experimental challenge. However, theoretically, when all Nt(i) are known, a maximum
likelihood estimate of m can be computed. Thus, this scenario serves as a theoretical benchmark
model. This function uses the golden section search method to calculate an MLE of mu (the mutation
rate), incorporating the measurements of Nt(i). The likelihood function is formed using the
Luria-Delbruck distribution as described by Lea and Coulson (1949).
}
\usage{
golden.benchmark.LD(data, Nt, mu.low = 1e-12, mu.up = 5e-07,
               tol = 1e-09, max.iter = 100, show.iter = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A vector containing the mutant counts.
}
  \item{Nt}{
A vector containing the final cell population sizes.
}
  \item{mu.low}{
A lower bound for the mutation rate.
}
  \item{mu.up}{
An upper bound for the mutation rate. Because exp(-m) will cause underflow when m exceeds 744,
it is imperative to ensure that (max Nt(i))*mu.up < 700.
}
  \item{tol}{
A parameter to control numerical accuracy.
}
  \item{max.iter}{
Maximum number of iterations allowable in the golden section search process.
}
  \item{show.iter}{
To exhibit the step-by-step golden section search process.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A real number as an estimate of the mutation rate, not the expected number of mutations (m).
}
\references{
E.A.  Lea, C.A. Coulson, The distribution of the numbers of mutants in
bacterial populations. J Genetics 49 (1949) 264-285.
}
\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
golden.benchmark.LD(demerec.data,rep(1.9e8,30))
}

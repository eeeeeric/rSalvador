\name{mcmc.MK.test}
\alias{mcmc.MK.test}

\title{
Comparison of mutation rates using an MCMC approach
}
\description{
The classic approach to the comparison of mutation rates is the likelihood ratio test (LRT),
which can be done with the LRT.MK functions, for example. The function allows for a Bayesian
approach to the comparison of mutation rates. In the first experiment, the expected number of
mutations per culture is N1*exp(b0), while in the second experiment that quantity is 
N2*exp(b0+b1). Therefore, the two mutation rates are exp(b0) and exp(b0+b1), respectively.
Equality of the two mutation rates is equivalent to b1=0. Note that the ratio of the two
mutation rates is mu2:mu1=exp(b1). For a general introduction to MCMC computing, see
Albert (2007); for an application of MCMC technique in the analysis of fluctuation assay data,
see Zheng (2011).
}


\usage{
mcmc.MK.test(x, y, N1, N2, w1 = 1, w2 = 1, Iter = 100, init.mu0 = -20,
    init.mu1 = 0.05, mu0 = -20.7, mu1 = 2.7, v0 = 30, v1 = 30, s0 = 0.5,
    s1 = 0.5, show.simu = FALSE)
}

\arguments{
  \item{x}{
A data vector containing mutant counts in the first experiment.
}
  \item{y}{
A data vector containing mutant counts in the second experiment.
}
  \item{N1}{
Average final cell population size of the first experiment.
}
  \item{N2}{
Average final cell population size of the second experiment.
}
  \item{w1}{
Relative fitness in the first experiment.
}
  \item{w2}{
Relative fitness in the second experiment.
}
  \item{Iter}{
Number of MCMC simulation cycles.
}
  \item{init.mu0}{
Initial guess of the logarithm of the mutation rate in the first experiment.
}
  \item{init.mu1}{
Initial guess of the difference between log(second mutation rate) and log(first
mutation rate).
}
  \item{mu0}{
Prior mean of log(first mutation rate)
}
  \item{mu1}{
Prior mean of log(second mutation rate) - log(first mutation rate).
}
  \item{v0}{
Prior variance of log(first mutation rate).
}
  \item{v1}{
Prior variance of mu1.
}
  \item{s0}{
Tuning parameter to control the acceptance rate of b0, the variance of a normal proposal
distribution.
}
  \item{s1}{
Tuning parameter to control the acceptance rate of b1, the variance of a normal proposal
distribution.
}
  \item{show.simu}{
A logic variable to control the display of the evolution of the MCMC simulation process.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
list(c(a1,a2), chain), where a1 and a2 are the two acceptance rates, and chain is a Iter by 2
matrix containing the simulated posterior distributions of mu0 and mu1.
}
\references{
J. Albert, Bayesian Computation Using R, Wiley, 2007.

Q. Zheng, A Bayesian approach for correcting for partial plating
in fluctuation experiments, Genetics Research 93 (2011) 351-356.
}

\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}
\note{

}

\seealso{
mcmc.MK
}
\examples{
mcmc.MK.test(demerec.data,ford.data,1.9e8,1.08e9,1,1,100)
}


\name{mcmc.MK}
\alias{mcmc.MK}

\title{
Estimating mutation rates under the MK model using an MCMC approach
}
\description{
An MCMC method is used to simulate the posterior distribution of the logarithm
of the mutation rate. A Gaussian prior distribution is assumed. For a general
introduction to MCMC computing, see Albert (2007); for an application of MCMC
technique in the analysis of fluctuation assay data, see Zheng (2011).
}
\usage{
mcmc.MK(y, Nt, w = 1, Iter = 100, init.mu = 2e-08, b0 = -15,
    v0 = 30, s0 = 0.8, show.simu = FALSE)
}

\arguments{
  \item{y}{
A data vector containing the number of mutants in a fluctuation experiment.
}
  \item{Nt}{
Average final cell population size.
}
  \item{w}{
Relative fitness defined as (division rate of mutants)/(division rate of wide-type cells).
}
  \item{Iter}{
Length of the simulated Markov chain.
}
  \item{init.mu}{
Initial guess of the mutation rate (not the expected number of mutations).
}
  \item{b0}{
Prior mean of log of the mutation rate.
}
  \item{v0}{
Variance in the Gaussian prior distribution.
}
  \item{s0}{
A tuning parameter to control acceptance rate; it is the variance of a Gaussian proposal
distribution.
}
  \item{show.simu}{
A logical variable to control the display of the evolution of the MCMC process.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
list(acc.rate,chain) where acc.rate is the acceptance rate and where chain is the simulated
Markov chain. Note that chain usually needs further processing, e.g., discarding a burn-in
period and thinning, etc.
}
\references{
J. Albert, Bayesian Computation Using R, Wiley, 2007.

Q. Zheng, A Bayesian approach for correcting for partial plating
in fluctuation experiments, Genetics Research 93 (2011) 351-356.
}
\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}

\seealso{
mcmc.MK.test
}
\examples{
mcmc.MK(demerec.data,1.9e8,w=1,Iter=100)
}


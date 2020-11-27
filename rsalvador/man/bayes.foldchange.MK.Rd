\name{bayes.foldchange.MK}
\alias{bayes.foldchange.MK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Bayesian Approach to Mutation Rate Fold Change
}
\description{
The function considers two fluctuation experiments. The mutation rate in one experiment is regarded as a baseline rate,
and the second mutation rate is expressed as mutation rate fold change relative to the baseline rate. Normal prior 
distributions are assumed for both parameters: the baseline rate and the fold change. Posterior distributions of the
logarithms of the two parameters are then simulated using an MCMC method. 
}
\usage{
bayes.foldchange.MK(y1, y2, N1, N2, w1 = 1, w2 = 1,
      s0 = 0.4, s1 = 0.5, init.base = 1e-08, init.fold = 1.6,
      v0 = 100, v1 = 100, iter = 1100, burn = 100, thin = 1,
      alpha = 0.05, short.out = TRUE, show.simu = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y1}{
Mutant count data for the first experiment.
}
  \item{y2}{
Mutant count data for the second experiment.
}
  \item{N1}{
Average final cell population size of the first experiment.
}
  \item{N2}{
Average final cell population size of the second experiment.
}

  \item{w1}{
Mutant relative fitness in the first experiment.
}
  \item{w2}{
Mutant relative fitness in the second experiment.
}
  \item{s0}{
Tuning parameter controlling the acceptance rate for the logarithm of the baseline rate, which is
the variance of a normal proposal distribution.
}
  \item{s1}{
Tuning parameter controlling the acceptance rate for the logarithm of the fold change, which is 
the variance of a normal proposal distribution.
}
  \item{init.base}{
An initial guess of the baseline mutation rate.
}
  \item{init.fold}{
An initial guess of the mutation rate fold change.
}
  \item{v0}{
Prior variance for the logarithm of the baseline mutation rate.
}
  \item{v1}{
Prior variance for the logarithm of the mutation rate fold change.
}
  \item{iter}{
Total number of MCMC iterations.
}
  \item{burn}{
The burn-in period used in preprocessing the MCMC chains.
}
  \item{thin}{
Thinning parameter used in preprocessing the MCMC chains.
}

  \item{alpha}{
The Bayesian confidence interval will have an approximate confidence coefficient of 1-alpha.
}

  \item{short.out}{
If short.out=TRUE, the function gives only a summary of the fold change posterior distribution.
Otherwise, it gives both Markov chains as well as summary statistics.
}
  \item{show.simu}{
A logic variable controlling the display of the MCMC simulation process.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
If short.out=TRUE, the function gives the median and a credible confidence interval for the fold change.
Otherwise, the output is of the form
list(accept.rate, list(stats1, stats2), list(chain1, chain2)), where accept.rate contains the two
acceptance rate, where stats1 and stats2 are the two summaries of the two chains, and where
chain1 and chain2 are the two preprocessed Markov chains
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Qi Zheng <qzheng@tamu.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
y1 = cairns.foster.data; y2 = demerec.data; N1 = 1.5e8; N2 = 1.9e8;
bayes.foldchange.MK(y1, y2, N1, N2, init.base = 2.5e-8,
              init.fold = 2.5, burn = 1, iter = 5)
} %% end of example



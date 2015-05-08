\name{plot.likelihood.LD}
\alias{plot.likelihood.LD}

\title{Plotting the Log-Likelihood Function under the Luria-Delbruck Model}

\description{This function plots the log-likelihood function of data under the Luria-Delbruck
model, using the stochastic formulation first proposed by Lea and Coulson (1948).
A history of the Lea-Coulson formulation is given in Zheng (1999).}

\usage{
plot.likelihood.LD(data, m.low = -1, m.up = -1, plot.pts = 30,
  lik.col = "black", mle.col = "red", title = "", x.lab = "Value of m",
  y.lab = "Log-likelihood", show.secant = TRUE)
}

\arguments{
  \item{data}{A vector of experimental data (non-negative integers).}

  \item{m.low}{plot.likelihood.LD plots the log-likelihood function for m ranging from
between m.low to m.up.}
  \item{m.up}{see above.}
  \item{plot.pts}{number of points used to plot the log-likelihood functions. A larger number
 of points yields a smoother graph, but takes more computing time.}

  \item{lik.col}{color for the log-likelihood function.}

  \item{mle.col}{color for the vertical bar that marks the maximum likelihood estimate of m.}

  \item{title}{title to be given to the graph.}

  \item{x.lab}{x label.}

  \item{y.lab}{y label.}

  \item{show.secant}{A secant line will be drawn at the maximum of the log-likelihood function.}
}

\references{
Q. Zheng,   Progress of a half century in the study of
the Luria-Delbruck distribution, Mathematical Biosciences 162 (1999) 1-32.

E.A. Lea and  C.A. Coulson, The distribution of the
numbers of mutants in bacterial populations.  J. Genetics 49 (1949) 264-285.
}


\author{Qi Zheng <qzheng@sph.tamhsc.edu>}

\seealso{plot.likelihood.LD.plating}

\examples{plot.likelihood.LD(demerec.data)}


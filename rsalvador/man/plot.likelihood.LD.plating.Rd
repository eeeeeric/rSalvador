\name{plot.likelihood.LD.plating}
\alias{plot.likelihood.LD.plating}

\title{Plotting the Log-Likelihood Function that Adjusts for Plating Efficiency}

\description{This functions plots the log-likelihood function that adjusts for plating
efficiency. The plating efficiency e is assumed to be known, the parameter m is regarded
as the argument of the log-likelihood function.}

\usage{
plot.likelihood.LD.plating(data, e, init.m=0, m.low = -1, m.up = -1,
   plot.pts = 30, lik.col = "black", mle.col = "red", title = "",
  x.lab = "Value of m", y.lab = "Log-likelihood", show.secant = TRUE)
}

\arguments{
  \item{data}{A vector of experimental data (non-negative integers).}
  \item{e}{Plating efficiency.}
  \item{init.m}{An initial guess of the parameter m.}
  \item{m.low}{It plots the log-likelihood function for m ranging from m.low to m.up.}
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

\author{Qi Zheng <qzheng@sph.tamhsc.edu>}

\seealso{plot.likelihood.LD}

\examples{plot.likelihood.LD.plating(luria.16.data,e=0.4)}


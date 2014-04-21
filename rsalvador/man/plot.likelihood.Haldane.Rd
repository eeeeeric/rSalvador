\name{plot.likelihood.Haldane}
\alias{plot.likelihood.Haldane}

\title{Plotting the Log-Likelihood Function Under the Haldane Model}

\description{This function plots the log-likelihood function of data under the Haldane
model. The argument of the log-likelihood function is the mutation rate mu. It is based
on the algorithm described in Zheng (2007).}

\usage{
plot.likelihood.Haldane(data,g, init.mu=0,  mu.low = -1, mu.up = -1,
 plot.pts = 30, lik.col = "black", mle.col = "red", title = "",
 x.lab = "Value of mu", y.lab = "Log-likelihood", show.secant = TRUE)
}

\arguments{
  \item{data}{A vector of experimental data (non-negative integers).}
  \item{g}{number of generations.}
  \item{init.mu}{a user specified initial guess of the mutation rate mu.}
  \item{mu.low}{plot.likelihood.Haldane plots the log-likelihood function 
  for mu between mu.low and mu.up.}
  \item{mu.up}{see above.}
  \item{plot.pts}{number of points used to plot the log-likelihood functions. A larger number
 of points yields a smoother graph, but takes more computing time.}

  \item{lik.col}{color for the log-likelihood function.}

  \item{mle.col}{color for the vertical bar that marks the maximum likelihood estimate of m.}

  \item{title}{title to be given to the graph.}

  \item{x.lab}{x label.}

  \item{y.lab}{y label.}

  \item{show.secant}{A secant line will be drawn at the maximum of the log-likelihood function.}
}


\references{ Q. Zheng, On Haldane's formulation of Luria and Delbruck mutation model,
Mathematical Biosciences 209 (2007) 500-513.}

\author{Qi Zheng <qzheng@sph.tamhsc.edu>}

\seealso{plot.likelihood.LD, plot.likelihood.LD.plating}

\examples{plot.likelihood.Haldane(niccum.data,g=25)}


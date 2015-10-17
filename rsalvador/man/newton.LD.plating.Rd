\name{newton.LD.plating}
\alias{newton.LD.plating}

\title{Computing Maximum Likelihood Estimate of m that Adjusts for Plating Efficiency} 

\description{This function computes the m.l.e. of m (the expected number of mutations)
when plating efficiency is imperfect. It employs a Newton-Raphson type algorithm
described in Zheng (2008).}

\usage{
newton.LD.plating(data, e, tol = 1e-08, init.m = 0, max.iter = 30,
 show.iter = FALSE)
}

\arguments{
  \item{data}{a vector of numbers of mutant colonies in a fluctuation experiment.}

  \item{e}{plating efficiency assumed to be known.}

  \item{tol}{tolerance parameter to control the number of iterations.}

  \item{init.m}{an initial guess about m can be specified, but is rarely needed.}

  \item{max.iter}{maximum number of iterations allowed.}

  \item{show.iter}{to view the iteration process.}
}

\value{A positive number.}

\references{
Q. Zheng, A note on plating efficiency in fluctuation experiments,
Mathematical Biosciences 216 (2008) 150-153.}

\author{Qi Zheng <qzheng@sph.tamhsc.edu>}

\note{
%%  ~~further notes~~
}

\seealso{newton.LD}

\examples{     
newton.LD.plating(wh.data[[1]],e=0.4)
}


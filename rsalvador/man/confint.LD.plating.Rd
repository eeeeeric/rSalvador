\name{confint.LD.plating}
\alias{confint.LD.plating}

\title{Confidence Intervals for m that Accounts for  Plating Efficiency}

\description{This function computes likelihood ratio-based confident intervals for m that
accounts for plating efficiency. The algorithm was described in Zheng (2008).}

\usage{
confint.LD.plating(data, e, alpha = 0.05, tol = 1e-9, init.m = 0,
 init.lower = 0, init.upper = 0, max.iter = 30, show.iter = FALSE)}

\arguments{
  \item{data}{A vector of experimental data, all integers.}

  \item{e}{Plating efficiency e (0<e<1)}

  \item{alpha}{A 1-alpha confidence interval is to be constructed.}

  \item{tol}{Tolerance parameter to control the number of iterations.}
  \item{init.m}{An initial guess for m can be supplied, but is rarely needed.}

  \item{init.lower}{An initial guess for the lower limit of the confidence interval, but is
rarely needed.}

  \item{init.upper}{An initial guess for the upper limit of the confidence interval, but is
rarely needed.}

  \item{max.iter}{Maximum number of iterations allowed.}

  \item{show.iter}{To view the iteration process.}

} %% end of arguments

\seealso{confint.LD}

\value{A pair of positive real numbers.}

\author{Qi Zheng <qzheng@sph.tamhsc.edu>}

\references{
Q. Zheng. A note on plating efficiency in fluctuation experiments,
Mathematical Biosciences, 216:150-153 (2008).}




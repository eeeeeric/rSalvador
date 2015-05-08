\name{newton.LD}
\alias{newton.LD}

\title{Computing  MLE of m Under the Luria-Delbruck Model }

\description{
This function computes maximum likelihood estimates of m, the expected number of mutations per
culture, using a Newton-Raphson type algorithm described in Zheng (2005). 
}

\usage{
newton.LD(data, phi = 1, tol = 1e-08, init.m = 0, max.iter = 30,
 show.iter = FALSE)
}

\arguments{
  \item{data}{a data vector, containing numbers of mutant cells in all test tubes.}

  \item{phi}{this parameter is defined by 1-N0/Nt, usually phi is taken as unity.}

  \item{tol}{tolerance parameter to control convergence.}

  \item{init.m}{an initial guess of m, it is rarely needed.}

  \item{max.iter}{maximum number of iterations.}

  \item{show.iter}{exhibition of the intermediate values of m during the iteration process.}
}

\value{A positive number.}

\references{
Q. Zheng, (2005). New algorithms for Luria-Delbruck fluctuation analysis,
Mathematical Biosciences 196:198-214 (2005).
}

\author{Qi Zheng <qzheng@sph.tamhsc.edu>}

\seealso{newton.LD.plating }

\examples{newton.LD(demerec.data,show.iter=TRUE)}


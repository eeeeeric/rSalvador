\name{confint.Haldane}
\alias{confint.Haldane}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Confidence Intervals for Mutation Rates Under the Haldane Model
}
\description{
It uses an algorithm described in Zheng (2007) to compute a confidence
interval for mu, the probability of mutation per cell division under the Haldane model.
}

\usage{
confint.Haldane(data, g, N0 = 1, alpha = 0.05, tol = 1e-06, init.mu = 0,
init.lower = 0, init.upper = 0, max.iter = 30, show.iter = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a vector of experimental data (non-negative integers).
}
  \item{g}{number of generations.
}
  \item{N0}{initial number of nonmutant cells.
}
  \item{alpha}{level of the confidence interval.
}

  \item{tol}{tolerance parameter to control the iteration process. 
}
  \item{init.mu}{an initial guess for mu.
}
  \item{init.lower}{ an initial guess for lower limit of m, rarely needed.
}
  \item{init.upper}{ an initial guess for upper limit of m, rarely needed.
}
  \item{max.iter}{maximum number of iterations.
}
  \item{show.iter}{to show the iteration progress.
}
}

\details{
The algorithm is described in Zheng (2007).
}
\value{a confidence interval for the mutation rate mu.
}
\references{
Q. Zheng, (2007). On Haldane's formulation of Luria and Delbruck's 
mutation model, Mathematical Biosciences 209:500-513.
}
\author{Qi Zheng <qzheng@sph.tamhsc.edu> }
\note{
}

\seealso{
confint.LD
confint.LD.plating
}
\examples{confint.Haldane(demerec.data,g=25,N0=90,show.iter=TRUE)
}



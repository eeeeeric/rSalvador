\name{newton.Haldane}
\alias{newton.Haldane}

\title{Computing MLE of the Mutation Rate Under the Haldane Model}

\description{
This function computes maximum likelihood estimates of mu, the probability of a mutation per
cell division, under the Haldane model. The Haldane model was described by Sarkar (1991), the
maximum likelihood estimate of mu is computed using a Newton-Raphson type iterative algorithm
described in Zheng (2007). 
}

\usage{
newton.Haldane(data, g, N0 = 1, tol = 1e-08, init.mu = 0, max.iter = 30,
 show.iter = FALSE)
}

\arguments{
  \item{data}{a data vector, containing numbers of mutant cells in all tubes.}

  \item{g}{the number of generations.}

  \item{N0}{the number of initial cells.}

  \item{tol}{tolerance parameter to control convergence.}

  \item{init.mu}{an initial guess of mu can be supplied.}

  \item{max.iter}{maximum number of iterations allowed.}

  \item{show.iter}{exhibition of the intermediate values of m during the iteration process.}
}

\value{A positive number.}

\references{

S. Sarkar (1991) Haldane's solution of the Luria-Delbruck distribution, Genetics 127, 257-261.

Q. Zheng (2007)  On Haldane's formulation of Luria and Delbruck's mutation model, Mathematical
Biosciences 209:500-513.
}

\author{Qi Zheng <qzheng@sph.tamhsc.edu>}

\seealso{newton.LD, newton.LD.plating}

\examples{newton.Haldane(niccum.data,g=25,show.iter=TRUE)}


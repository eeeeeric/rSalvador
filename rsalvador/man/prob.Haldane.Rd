\name{prob.Haldane}
\alias{prob.Haldane}

\title{Computing the Haldane Distribution}

\description{The Haldane mutation model is described in Sarkar (1991). Under the Haldane
model, cell division is synchronized. When a nonmutant cell divides, the probability is mu
that one daughter cell will be a mutant cell. This function computes the Haldane distribution
using an algorithm described in Zheng (2007).}

\usage{prob.Haldane(gen, mu, n, N0)}

\arguments{
  \item{gen}{number of generations.}

  \item{mu}{the mutation rate.}

  \item{n}{It calculates the first n+1 probabilities, starting from the zeroth probability.}
  
  \item{N0}{The number of initial nonmutant cells.}

} %% end of arguments



\value{A vector contain the n+1 probabilities.}

\references{
S. Sarkar (1991) Haldane's solution of the Luria-Delbruck distribution,
Genetics 127, 257-261.

Q. Zheng, On Haldane's formulation of Luria and Delbruck mutation model,
Mathematical Biosciences 209 (2007) 500-513.}

\author{Qi Zheng <qzheng@sph.tamhsc.edu>}

\note{
%%  ~~further notes~~
}


\seealso{prob.LD}

\examples{prob.Haldane(12,0.0003,20,1)}


\name{prob.LD.plating}
\alias{prob.LD.plating}

\title{
Calculating the Luria-Delbruck Mutant Cell Distribution that Adjusts for Plating Efficiency
}
\description{
This function computes the mutant cell distribution that adjusts for plating efficiency.
Several authors contributed to the computation of this distribution, among them are
Armitage (1952), Stewart (1991) and Jones (1994). Inspired by these important results,
Zheng (2008) proposed a computationally feasible algorithm, which prob.LD.plating adopts.
}

\usage{
prob.LD.plating(m, e, n)
}

\arguments{
  \item{m}{the expected number of mutations per culture.}

  \item{e}{plating efficiency ranging between 0 and 1.}

  \item{n}{the function computes the first n+1 probabilities, i.e., P0, P1, ..., Pn.}
}

\value{A vector of the first n+1 probabilities.}

\references{

P. Armitage, The statistical theory of bacterial populations subject to
mutation, J. Royal Statistical Society, ser. B, 14:1-44, 1952.

M. E. Jones, Luria-Delbruck fluctuation experiments; Accounting simultaneously
for plating efficiency and differential growth rate, Journal of Theoretical
Biology, 166: 355-363, 1994.

F.M. Stewart, Fluctuation analysis: the effect of plating
efficiency, Genetica 54: 51-55, 1991.

Q. Zheng, A note on plating efficiency in fluctuation experiments,
Mathematical Biosciences, 216:150-153, 2008.
}

\author{Qi Zheng <qzheng@sph.tamhsc.edu>}

\seealso{prob.LD}

\examples{prob.LD.plating(5.8,0.2,20)}

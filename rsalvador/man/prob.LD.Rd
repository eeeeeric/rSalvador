\name{prob.LD}
\alias{prob.LD}

\title{Computing the Luria-Delbruck Mutant Cell Distribution}

\description{This function computes the Luria-Delbruck distribution function. This particular
version of the Luria-Delbruck distribution was first proposed by Lea and Coulson (1949), which
was not an exact form of the intended distribution. Bartlett later (see, e.g. Bartlett, 1955)
modified this distribution by introducing the parameter phi, making the distribution exact.
A historical account of this distribution was given by Zheng (1999). The algorithm for computing
this distribution was given by Ma et al. (1992), although the parameter phi was not considered
at the time.} 

\usage{prob.LD(m, phi, k)}

\arguments{
  \item{m}{expected number of mutations.}

  \item{phi}{This parameter is defined as 1-N0/Nt.}

  \item{k}{It computes the first k+1 probabilities, starting with the zeroth probability.}

} %% end of arguments


\details{Most implementations did not include the clock parameter phi.}

\value{A vector of the k+1 probabilities.}

\references{

M.S. Bartlett, An Introduction to Stochastic Processes, Cambridge University Press,
London, 1955 (pp.132-135).

D.E. Lea and  C.A. Coulson, The distribution of the
numbers of mutants in bacterial populations, J. Genetics 49 (1949) 264-285.

W.T. Ma, G. vH. Sandri, S. Sarkar, Analysis of the Luria-Delbruck distribution
using discrete convolution powers, J. Appl. Prob. 29 (1992) 255-267.

Q. Zheng, Progress of a half century in the study of the Luria-Delbruck
distribution, Mathematical Biosciences 162 (1999) 1-32.
}

\author{Qi Zheng <qzheng@sph.tamhsc.edu>}

\note{
%%  ~~further notes~~
}


\seealso{prob.LD.plating}

\examples{prob.LD(5.8,1,12)}


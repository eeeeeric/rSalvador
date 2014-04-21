\name{simu.cultures}
\alias{simu.cultures}

\title{Simulating Fluctuation Experiments}

\description{This function simulates the numbers of mutant cells in all cultures prior to
plating. The mutation model used here was proposed by Mandelbrot (1966) and Koch (1974).
Under this model, the nonmutant cell growth rate, b1, can differ from that of mutant cells,
b2. The expected number of mutations is m=mu*(Nt-N0)/b1. The simulation algorithm is
described in Zheng (2002).}

\usage{simu.cultures(n, mu, b1, b2, N0, Nt)}

\arguments{

\item{n}{total number of tubes in an experiment.}

\item{mu}{the mutation rate, not the expected number of mutations.}

\item{b1}{wild-type cell growth rate.}

\item{b2}{mutant cell growth rate.}

\item{N0}{initial number of wild-type cells in a tube.}

\item{Nt}{final number of wild-type cells in a tube.}

} %% end of arguments

\value{A vector of n non-negative integers.}

\references{

A.L. Koch, Mutation and growth rates from Luria-Delbruck fluctuation 
tests, Mutation Research 95:129-143 (1982).

B. Mandelbrot, A population birth-and-mutation process, I: Explicit distributions for
the number of mutants in an old culture of bacteria. J. Appl. Prob. 11:437-444 (1974).

Q. Zheng. Statistical and algorithmic methods for fluctuation analysis with SALVADOR
as an implementation. Mathematical Biosciences, 176:237-252 (2002).
}

\author{Qi Zheng <qzheng@sph.tamhsc.edu>}

\note{
}

\examples{simu.cultures(30,10^-8,1,1,20,5*10^8)}


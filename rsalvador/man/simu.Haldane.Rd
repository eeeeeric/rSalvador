\name{simu.Haldane}
\alias{simu.Haldane}

\title{Simulating a Fluctuation Experiment Under the Haldane Model}

\description{This function simulates the numbers of mutant cells in a fluctuation experiment
under the Haldane model as described in Sarkar (1991) and Zheng (2007).}

\usage{simu.Haldane(gen, mu, culture)}

\arguments{

\item{gen}{The number of generations of cell division.}

\item{mu}{The mutation rate}

\item{culture}{The total number of cultures in an experiment.}


} %% end of arguments

\value{A vector of non-negative integers.}

\references{

S. Sarkar (1991) Haldane's solution of the Luria-Delbruck distribution, Genetics 127, 257-261.

Q. Zheng (2007)  On Haldane's formulation of Luria and Delbruck's mutation model, Mathematical
Biosciences 209:500-513.
}

\author{Qi Zheng <qzheng@sph.tamhsc.edu>}

\seealso{newton.LD, newton.LD.plating, simu.Kimmel.}

\note{ }

\examples{simu.Haldane(9,0.03,22)}


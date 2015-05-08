\name{simu.Kimmel}
\alias{simu.Kimmel}

\title{Simulating a Fluctuation Experiment Under a variation of the Haldane Model}

\description{This function simulates the numbers of mutant cells in a fluctuation experiment
under a variation of the Haldane model. Unlike the original Haldane model, this model allows
both daughter cells to be mutant cells when a nonmutant cell divides. (Each daughter cell has
a probability mu to be a mutant cell.) This model was studied by Kimmel and Axelrod (1994).}

\usage{simu.Kimmel(gen, mu, culture)}

\arguments{

\item{gen}{The number of generations.}

\item{mu}{The mutation rate.}

\item{culture}{The total number of cultures in an experiment.}


} %% end of arguments

\value{A vector of non-negative integers.}

\references{
M. Kimmel, D.E. Axelrod, 1994. Fluctuation test for two-stage
mutations: application to gene amplification, Mutation Research 306, 45-60.
}

\author{Qi Zheng <qzheng@sph.tamhsc.edu>}

\seealso{simu.Haldane}

\note{}

\examples{simu.Kimmel(9,0.03,22)}


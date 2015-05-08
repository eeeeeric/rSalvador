\name{log.likelihood.Haldane}
\alias{log.likelihood.Haldane}

\title{Computing the Log-Likelihood Function Under the Haldane Model
}
\description{
This function computes the log-likelihood function for data under the Haldane model. The
algorithm is based on Zheng (2007). 
}
\usage{
log.likelihood.Haldane(data, g, mu, N0 = 1)
}

\arguments{
  \item{data}{
A vector of data.
}
  \item{g}{
Number of generations.
}
  \item{mu}{
The mutation rate.
}
  \item{N0}{
Initial number of nonmutant cells.
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\references{
Q. Zheng (2007)  On Haldane's formulation of Luria and Delbruck's mutation model, Mathematical
Biosciences 209:500-513.
}

\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}

\seealso{
log.likelihood.LD, log.likelihood.LD.plating
}
\examples{
log.likelihood.Haldane(niccum.data,25,1e-6)
}

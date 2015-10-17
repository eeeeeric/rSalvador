\name{LRT.LD}
\alias{LRT.LD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood Ratio Test under the Luria-Delbruck Model
}
\description{
The function computes the likelihood ratio test statistic under the Luria-Delbruck
model. Under this model, relative fitness w=1. The algorithm is a generalization
of a simpler algorithm given in Zheng (2015).
}
\usage{
LRT.LD(X1, X2, R = 1, phi1 = 1, phi2 = 1, init.mc = 0, init.m1 = 0,
     init.m2 = 0, tol = 1e-09, show.iter = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X1}{
Data from the first experiment.
}
  \item{X2}{
%%     ~~Describe \code{X2} here~~
}
  \item{R}{
Data from the second experiment.
}
  \item{phi1}{
"Time" parameter for the first experiment, usually set to unit.
}
  \item{phi2}{
"Time" parameter for the second experiment, usually set to unit.
}
  \item{init.mc}{
Initial guess of m under the null hypothesis -- the combined model.
}
  \item{init.m1}{
Initial guess of m for the first experiment.
}
  \item{init.m2}{
Initial guess of m for the second experiment.
}
  \item{tol}{
Tolerance parameter to control numerical accuracy.
}
  \item{show.iter}{
Logical variable to control the display of the numerical iteration process.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A pair of real numbers: (test statistic, p-value).
}

\references{
Q. Zheng, Methods for comparing mutation rates using fluctuation assay data, Mutation
Research/Fundamental and Molecular Mechanisms of Mutagenesis 777 (2015) 20-22.
}
\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}

\seealso{
LRT.MK
}
\examples{
x=simu.cultures(20,1e-8,1,1,50,2e8)
y=simu.cultures(20,1e-8,1,1,50,3e8)
LRT.LD(x,y,1.5)
}

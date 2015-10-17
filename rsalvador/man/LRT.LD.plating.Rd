\name{LRT.LD.plating}
\alias{LRT.LD.plating}

\title{
Likelihood Ratio Test under the Classic Luria-Delbruck Model with Partial Plating
}
\description{
The function computes the likelihood ratio test statistics. The algorithm is a
generalization of a simpler algorithm given in Zheng (2015).
}
\usage{
LRT.LD.plating(X1, X2, R = 1, e1 = 0.1, e2 = 0.1, init.mc = 0,
      init.m1 = 0, init.m2 = 0, tol = 1e-09, show.iter = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X1}{Data from the first experiment.
}
  \item{X2}{
Data from the second experiment.
}
  \item{R}{
Ratio of final cell population sizes, R=(Nt of second experiment)/(Nt of first experiment).
}
  \item{e1}{
Plating efficiency in the first experiment.
}
  \item{e2}{
Plating efficiency in the second experiment.
}
  \item{init.mc}{
Initial guess of m under the null hypothesis.
}
  \item{init.m1}{
Initial guess for m1, m of the first experiment.
}
  \item{init.m2}{
Initial guess for m2, m of the second experiment.
}
  \item{tol}{
Tolerance parameter to control numerical accuracy.
}
  \item{show.iter}{
Logical variable to allow for display of the numerical iteration process.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A pair of real numbers: (test statistic, p-value).
}
\references{
Q. Zheng, Methods for comparing mutation rates using fluctuation assay data, Mutation
Research/Fundamental and Molecular Mechanisms of Mutagenesis 777 (2015) 20-22.
}
\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}
\note{
%%  ~~further notes~~
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
R=0.5/2.3
LRT.LD.plating(wh.data[[1]],wh.data[[2]],R=R,e1=0.4,e2=0.4)
}

\name{LRT.MK}
\alias{LRT.MK}

\title{
Likelihood Ratio Test under the Mandelbrot-Koch Model
}
\description{
The function computes the likelihood ratio test statistics. The algorithm is a
generalization of a simpler algorithm given in Zheng (2015).
}
\usage{
LRT.MK(X1, X2, R = 1, w1 = 1, w2 = 1, init.mc = 0, init.m1 = 0,
  init.m2 = 0, tol = 1e-09, show.iter = FALSE)
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
  \item{w1}{
Relative fitness in the first experiment.
}
  \item{w2}{
Relative experiment in the second experiment.
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
x=simu.cultures(20,1e-8,1,1.2,50,2e8)
y=simu.cultures(20,1e-8,1,1.4,50,3e8)
LRT.MK(x,y,R=1.5,w1=1.2,w2=1.4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

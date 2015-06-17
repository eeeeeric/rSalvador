\name{likely.average}
\alias{likely.average}

\title{
A family of Estimators of m Based on Delbruck's Concept of the Likely Average
}
\description{
The concept of the likely average (Luria and Delbruck 1943) spawns a family of estimators.
The well-known method of Lea and Coulson (1949) is a member of that family. However, most
of these estimators are given solely for their pedagogical values and historical interest.
}
\usage{
likely.average(data, b = 1.24, classic = FALSE, init.m = 0, tol = 1e-09,
    max.iter = 25, use.median = TRUE, show.iter = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A vector of experimental data.
}
  \item{b}{
Setting b=1.24 (default) gives the Lea-Coulson estimator when the variable classic is set to
False.
}
  \item{classic}{
If classic is true, the estimator depends on the number of cultures in an experiment. That was
Delbruck's original idea.
}
  \item{init.m}{
An initial guess of m.
}
  \item{tol}{
Tolerance limit to control accuracy.
}
  \item{max.iter}{
Maximum number of iterations in the Newton procedure.
}
  \item{use.median}{
Either the sample median (default) or the sample mean is used in the equation.
}
  \item{show.iter}{
Show iteration history.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A real number as an estimate of m.
}
\references{
E. A.  Lea, C. A.  Coulson, The distribution of the numbers of mutants in bacterial
populations.  J Genetics 49 (1949) 264-285.

S. E. Luria, M. Delbruck, Mutations of bacteria from virus sensitivity to virus
resistance. Genetics 28 (1943) 491-511.
}
\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
newton.LD
}
\examples{
likely.average(demerec.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

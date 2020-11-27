\name{compare.LD}
\alias{compare.LD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compare mutation rates with a common Nt
}
\description{
This function compares the expected number of mutations per test tube (m) between two fluctuation
experiments. Because a common Nt is assumed, it compares the mutation rates between the two
experiments.
}
\usage{
compare.LD(x1, x2, init.m0 = 0, init.m1 = 0, init.m2 = 0,
           phi = 1, show.iter=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
a data vector representing mutant counts in the first experiment.
}
  \item{x2}{
a data vector representing mutant counts in the second experiment.
}
  \item{init.m0}{
initial value for estimating a common m shared by the two experiments.
}
  \item{init.m1}{
initial value for estimating m for the first experiment.
}
  \item{init.m2}{
initial value for estimating m for the second experiment.
}
  \item{phi}{
the common clock parameter, often set to unity.
}
  \item{show.iter}{
a logic variable controlling the display of the iterative numerical process.
}



}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a pair of positive real numbers: (chi-squared test statistic, p-value).
}
\references{
Q Zheng, Methods for comparing mutation rates using fluctuation assay data,
Mutat. Res. Fundam. Mol. Mech. Mutagen. 777 (2015) 20-22.
}
\author{
Qi Zheng <qzheng@.tamu.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
x=simu.cultures(15, 1e-7, 1, 1, 5, 1e7)
y=simu.cultures(15, 3e-7, 1, 1, 5, 1e7)
compare.LD(x, y)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{simu.Kimmel2}
\alias{simu.Kimmel2}

\title{Simulating the Numbers of Mutant Cells for the Last Two Generations Under the Kimmel Model}

\description{
This function is the same as simu.Kimmel, except that it returns the numbers of mutant cells
for the last two generations.
}
\usage{
simu.Kimmel2(gen, mu, culture = 1)
}

\arguments{
  \item{gen}{The last generation. }
  \item{mu}{The mutation rate.}
  \item{culture}{The number of cultures to be simulated. }
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A pair of non-negative integer vectors.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
}
\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
simu.Haldane2
}
\examples{
simu.Kimmel2(9,0.08,20)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

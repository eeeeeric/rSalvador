\name{simu.Bartlett}
\alias{simu.Bartlett}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulating a Fluctuation Experiment Under the Bartlett Model
}
\description{
This function simulates the dynamics of (Xt,Yt), where Xt is the number of nonmutant cells
at time t, and where Yt is the number of mutant cells at time t. This model was first proposed
by Bartlett (1955) and further studied by Zheng (2008). The simulation algorithm is similar
to that given on page 281 of Renshaw (1991).
}
\usage{
simu.Bartlett(b1, b2, mu, N0, T, max.events = 1e+10, show.growth = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b1}{
Wild cell growth rate
}
  \item{b2}{
Mutant cell growth rate
}
  \item{mu}{
The mutation parameter "mu", giving a mutation rate of alpha=mu/(b1+mu).
}
  \item{N0}{
Initial number of wild cells.
}
  \item{T}{
Terminating time.
}
  \item{max.events}{
Maximum number of simulation iterations; if exceeded, return value will be NA.
}
  \item{show.growth}{ 
A logical variable; if set to TRUE, cellular dynamics will be shown.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A pair of non-negative integers in the form of (wild,mutant).
}
\references{

M.S. Bartlett, An Introduction to Stochastic Processes, Cambridge University Press, (1955).

E. Renshaw, Modeling Biological Populations in Space and Time, Cambridge University Press (1991).

Q. Zheng, On Bartlett's formulation of the Luria-Delbruck mutation model, Mathematical
Biosciences 215:48-54 (2008).


}
\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
simu.cultures, simu.Haldane
}
\examples{
simu.Bartlett(1,1,5e-3,1,5.8,show.growth=TRUE)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

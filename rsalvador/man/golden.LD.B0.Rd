\name{golden.LD.B0}
\alias{golden.LD.B0}

\title{
Accounting for Varying Nt Under the Luria-Delbruck Model
}
\description{
This function uses the B0 distribution (Zheng 2010, 2011) to account for variation in Nt, the final
cell population size prior to plating. The function employs the golden section search algorithm
to find a value of m that maximizes the log likelihood functions.
}
\usage{
golden.LD.B0(data, cv = 0.1, m.low = 0.1, m.up = 30, tol = 1e-8,
   max.iter = 60, show.iter = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A vector containing the mutant count data.
}
  \item{cv}{
Known (independently estimated) coefficient of variation for Nt, final cell population size.
}
  \item{m.low}{
A lower starting value of m for the golden section search.
}
  \item{m.up}{
An upper starting value of m for the golden section search.
}
  \item{tol}{
A tolerance parameter to control numerical precision.
}
  \item{max.iter}{
Maximum number of iterations.
}
  \item{show.iter}{
A logical variable controlling the exhibition of the numerical iteration process.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A positive real number as an estimate of m, the expected number of mutations per culture.
}
\references{
Q. Zheng, A new discrete distribution induced by the
Luria-Delbruck mutation model, Statistics 44 (2010) 529-540.

Q. Zheng, A Bayesian two-level model for fluctuation assay,
Genetica 139 (2011) 1409-1416.
}
\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
newton.LD, newton.LD.plating
}
\examples{
golden.LD.B0(demerec.data,cv=0.064)
}

\name{newton.foldchange.LD.plating}
\alias{newton.foldchange.LD.plating}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mutation Rate Fold Change Under a Luria-Delbruck Model with Imperfect Plating Efficiency
}
\description{
The algorithm regards the mutation rate in the first experiment as a baseline. The mutation rate fold change is estimated
by the maximum likelihood method. This function then uses the profile likelihood method to compute confidence limits for
the mutation rate fold change. This approach is similar to the one described in Zheng (2005).
}
\usage{
newton.foldchange.LD.plating(x, y, Nx, Ny, e1 = 0.1, e2 = 0.1, tol = 1e-08, init.base = 1.5e-08, 
                             init.fold = 1.6, max.iter = 30, no.log = TRUE, show.iter = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Mutant count data from the first experiment.
}
  \item{y}{
Mutant count data from the second experiment.
}
  \item{Nx}{
Final cell number Nt in the first experiment.
}
  \item{Ny}{
Final cell number Nt in the second experiment.
}
  \item{e1}{
Plating efficiency in the first experiment.
}
  \item{e2}{
Plating efficiency in the second experiment.
}
  \item{tol}{
Tolerance controlling the numerical iteration process.
}
  \item{init.base}{
An initial guess of the baseline mutation rate, that is, the mutation rate in the first experiment.
}
  \item{init.fold}{
An initial guess of the mutation rate fold change.
}
  \item{max.iter}{
Maximum number of iterations allowed in the numerical iteration process.
}
  \item{no.log}{
The algorithm computes the logarithm of a fold change, but it returns the fold change directly unless no.log is FALSE.
}
  \item{show.iter}{
A logic variable controlling the display of the numerical iteration process.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of the form (baseline mutation rate, fold change)
}
\references{
Q. Zheng, New algorithms for Luria-Delbruck fluctuation analysis,
Mathematical Biosciences 196 (2005) 198-214.
}
\author{
Qi Zheng <qzheng@tamu.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
confint.foldchange.LD.plating
}
\examples{
y1 = wh.data[[1]]; y2 = wh.data[[2]]; N1 = 2.3e8; N2 = 0.5e8;
newton.foldchange.LD.plating(y1, y2, N1, N2, e1 = 0.4, e2 = 0.4, init.base = 2e-8, show.iter = TRUE)
}

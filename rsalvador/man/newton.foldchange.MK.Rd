\name{newton.foldchange.MK}
\alias{newton.foldchange.MK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating Mutation Rate Fold Change Under an MK Model
}
\description{
The algorithm regards the mutation rate in the first experiment as a baseline. The mutation rate fold change
is estimated by the maximum likelihood method. This approach is similar to the one described in Zheng (2005).
}
\usage{
newton.foldchange.MK(x, y, Nx, Ny, w1 = 1.0, w2 = 1.0, tol = 1e-09, init.base = 1.5e-08, init.fold = 1.6,
                     max.iter = 30, no.log = TRUE, show.iter = FALSE)
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
  \item{w1}{
Relative fitness of mutants in the first experiment.
}
  \item{w2}{
Relative fitness of mutants in the second experiment.
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
Q. Zheng, New algorithms for Luria-Delbruck fluctuation analysis, Mathematical Biosciences 196 (2005) 198-214.
}
\author{
Qi Zheng <qzheng@tamu.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
confint.foldchange.MK   
}
\examples{
y1 = cairns.foster.data; y2 = demerec.data; N1 = 1.5e8; N2 = 1.9e8;
newton.foldchange.MK(y1, y2, N1, N2, w1 =0.95, w2 = 1.05, init.base = 2.5e-8, init.fold = 2.5, show.iter = TRUE)
}

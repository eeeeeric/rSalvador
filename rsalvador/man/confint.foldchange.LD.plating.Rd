\name{confint.foldchange.LD.plating}
\alias{confint.foldchange.LD.plating}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

CI for Mutation Rate Fold Change Under a Luria-Delbruck Model with Imperfect Plating Efficiency

}
\description{
The algorithm regards the mutation rate in the first experiment as a baseline. The mutation rate fold change is estimated
by the maximum likelihood method. This function then uses the profile likelihood method to compute confidence limits for
the mutation rate fold change. This approach is similar to the one described in Zheng (2005).
}
\usage{
confint.foldchange.LD.plating(x, y, Nx, Ny, e1 = 0.1, e2 = 0.1,
       alpha = 0.05, init.base = 1e-08, init.fold = 1.6, 
       init.low.base = -9, init.up.base = -9, init.low.fold = -9,
       init.up.fold = -9, max.iter = 30, tol = 1e-09, show.iter = FALSE)
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
  \item{alpha}{
The confidence interval for mutation rate foldl change will have an approximate confidence coefficient of 1-alpha.
}
  \item{init.base}{
An initial guess of the baseline mutation rate, that is, the mutation rate in the first experiment.
}
  \item{init.fold}{
An initial guess of the mutation rate fold change.
}
  \item{init.low.base}{
An initial guess of the lower confidence limit for the baseline mutation rate.
}

  \item{init.up.base}{
An initial guess of the upper confidence limit for the baseline mutation rate.
}
  \item{init.low.fold}{
An initial guess of the left confidence limit for the mutation rate fold change.
}
  \item{init.up.fold}{
An initial guess of the right confidence limit for the mutation rate fold change.
}
  \item{max.iter}{
Maximum number of iterations allowed in the numerical iteration process.
}
  \item{tol}{
Tolerance controlling the numerical iteration process.
}
  \item{show.iter}{
A logic variable controlling the display of the numerical iteration process.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of the form (lower confidence limit, fold change, upper confidence limit)
}
\references{
Q. Zheng, New algorithms for Luria-Delbruck fluctuation analysis,
Mathematical Biosciences 196 (2005) 198-214.
}
\author{
Qi Zheng <qzheng@tamu.edu>
}
\note{
newton.foldchange.LD.plating
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
y1 = wh.data[[1]]; y2 = wh.data[[2]]; N1 = 2.3e8; N2 = 0.5e8;
confint.foldchange.LD.plating(y1, y2, N1, N2, 
   e1 = 0.4, e2 = 0.4, init.base = 2e-8, show.iter = TRUE)
}



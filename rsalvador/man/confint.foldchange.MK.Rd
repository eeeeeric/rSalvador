\name{confint.foldchange.MK}
\alias{confint.foldchange.MK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence Interval for Mutation Rate Fold Change Under an MK Model
}
\description{
The algorithm regards the mutation rate in the first experiment as a baseline. The mutation rate fold change is
estimated by the maximum likelihood method. This function then uses the profile likelihood method to compute
confidence limits for the mutation rate fold change. This approach is similar to the one described in Zheng (2005).
}
\usage{
confint.foldchange.MK(x, y, Nx, Ny, w1 = 1.0, w2 = 1.0, alpha = 0.05,
     init.base = 1e-08, init.fold = 1.6, init.low.base = -9,
     init.up.base = -9, init.low.fold = -9, init.up.fold = -9,
     max.iter = 30, tol = 1e-09, show.iter = FALSE)
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
Final cell number Nt for the first experiment.
}
  \item{Ny}{
Final cell number Nt for the second experiment.
}
  \item{w1}{
Relative fitness of mutants in the first experiment.
}
  \item{w2}{
Relative fitness of mutants in the second experiment.
}
  \item{alpha}{
The confidence interval will have an approximate confidence coefficient of 1-alpha.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
y1 = cairns.foster.data; y2 = demerec.data; N1 = 1.5e8; N2 = 1.9e8;
confint.foldchange.MK(y1, y2, N1, N2, w1 =0.95, w2 = 1.05,
     init.base = 2.5e-8, init.fold = 2.5, show.iter = TRUE)
}

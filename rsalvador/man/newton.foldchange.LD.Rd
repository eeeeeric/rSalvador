\name{newton.foldchange.LD}
\alias{newton.foldchange.LD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computing mutation rate fold change
}
\description{
The algorithm regards the mutation rate in the first experiment as a baseline. The mutation rate fold change is estimated
by the maximum likelihood method.
}
\usage{
newton.foldchange.LD(x, y, Nx, Ny, phi.x = 1, phi.y = 1, tol = 1e-08, init.base = 1.5e-08,
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
Final cell number Nt for the first experiment.
}
  \item{Ny}{
Final cell number Nt for the second experiment.
}
  \item{phi.x}{
The phi parameter in the mutant distribution for the first experiment, it is almost always 1.0 in practice.
}
  \item{phi.y}{
The phi parameter in the mutant distribution for the second experiment, it is almost always 1.0 in practice.
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
%% ~put references to the literature/web site here ~
}
\author{
Qi Zheng <qzheng@tamu.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
confint.foldchange.LD
}
\examples{
y1=cairns.foster.data; y2=demerec.data; N1=1.5e8; N2=1.9e8;
newton.foldchange.LD(y1,y2,N1,N2,init.base=2.5e-8,init.fold=2.5,show.iter=TRUE)
}

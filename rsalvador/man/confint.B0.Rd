\name{confint.B0}
\alias{confint.B0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood ratio-based confidence interval for m, under the Luria-Delbruck model, 
accounting for variation in Nt
}
\description{
This function computes a likelihood ratio-based confidence interval for m, taking into consideration
of the variation in Nt, the cell population size immediately prior to plating.
The algorithm is given in Zheng (2016), which is adapted from an algorithm given in Zheng (2011).
}
\usage{
confint.B0(data, cv = 0.1, alpha = 0.05, tol = 1e-08, init.m = 5.1,
   init.lower = 0, init.upper = 0, max.iter = 30, show.iter = FALSE)
}

\arguments{
  \item{data}{
A data vector representing the numbers of mutants in the test tubes.
}
  \item{cv}{
Coefficient of variation indicating the variation in Nt.
}
  \item{alpha}{
Significance level.
}
  \item{tol}{
Tolerance parameter to control numerical accuracy.
}
  \item{init.m}{
Initial guess for the parameter m.
}
  \item{init.lower}{
Initial guess for the lower limit of the confidence interval.
}
  \item{init.upper}{
Initial guess for the upper limit of the confidence interval.
}
  \item{max.iter}{
Maximum number of numerical iterations.
}
  \item{show.iter}{
To show the numerical iteration process.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A pair of real numbers representing a confidence interval for m.
}
\references{
Q. Zheng, A Bayesian two-level model for fluctuation assay, Genetica 139 (2011) 1409-1416.

Q. Zheng, A second look at the final number of cells in a fluctuation experiment, Journal
of Theoretical Biology 401 (2016) 54-61.
}
\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
newton.B0,confint.LD, confint.MK, confint.LD.plating
}
\examples{
confint.B0(demerec.data,cv=0.1)
}

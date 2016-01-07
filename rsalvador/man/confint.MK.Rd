\name{confint.MK}
\alias{confint.MK}

\title{
Likelihood-Based Confidence Interval under the Mandelbrot-Koch Model
}
\description{
This function computes confidence intervals for m under the Mandelbrot-Koch model
that allows for differential growth between mutants and nonmutants. A value of
the relative fitness w is obtained from a fitness assay. The algorithm used is a
specialization of a more general algorithm described in Zheng (2005).
}
\usage{
confint.MK(data, w = 1, alpha = 0.05, tol = 1e-08, init.m = 0,
 init.lower = 0, init.upper = 0, max.iter = 30, show.iter = FALSE)
}

\arguments{
  \item{data}{A vector of experimental data (non-negative integers as mutant counts).

}
  \item{w}{
Relative fitness defined as (mutant growth rate)/(nonmutant growth rate).
}
  \item{alpha}{The confidence coefficient is 1-alpha. 

}
  \item{tol}{
Tolerance parameter to control numerical accuracy.
}
  \item{init.m}{
Initial guess for the parameter m, often the default is sufficient.
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
A pair of real numbers as the confidence limits.
}
\references{
Q. Zheng, New algorithms for Luria-Delbruck fluctuation analysis,
Mathematical Biosciences 196 (2005) 198-214.
}
\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
confint.LD
}
\examples{
confint.MK(demerec.data,w=0.75)
}


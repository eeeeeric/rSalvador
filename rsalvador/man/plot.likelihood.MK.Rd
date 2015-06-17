\name{plot.likelihood.MK}
\alias{plot.likelihood.MK}

\title{
Plotting the Log-Likelihood Function under the Mandelbrot-Koch Model
}
\description{
This function plots the log-likelihood function under the Mandelbrot-Kock mutation model.
A value for the relative fitness parameter w is often obtained via a fitness assay. The
algorithm for computing the maximum likelihood estimate of m is a specialization of an
algorithm given in Zheng (2005).
}
\usage{
plot.likelihood.MK(data, w = 1, m.low = -1, m.up = -1, init.m = 0,
  plot.pts = 30, lik.col = "black", mle.col = "red", title = "",
  x.lab = "Value of m", y.lab = "Log-likelihood", show.secant = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A vector of experimental data.
}
  \item{w}{
Relative fitness.
}
  \item{m.low}{
The range of m for the plot is from m.low to m.up.
}
  \item{m.up}{
See above.
}
  \item{init.m}{
An initial guess of the parameter m.
}
  \item{plot.pts}{
Number of points used to plot the log-likelihood function.
}
  \item{lik.col}{
Color for the log-likelihood function.
}
  \item{mle.col}{
Color of the vertical line segment indicating the maximum likelihood estimate of m.
}
  \item{title}{
Plot title.
}
  \item{x.lab}{
Plot x-label.
}
  \item{y.lab}{
Plot y-label.
}
  \item{show.secant}{
A secent line will be drawn at the maximum of the log-likelihood function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{None. A plot is created.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
plot.likelihood.MK(demerec.data,w=1.2)
}

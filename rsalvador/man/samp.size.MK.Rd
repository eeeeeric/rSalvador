\name{samp.size.MK}
\alias{samp.size.MK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating sample size under the Mandelbrot-Koch model
}
\description{
This function calculates required sample size based on rough prior knowledge of m and w. The user can
estimate the required sample size based on the psi score as defined in Zheng (2017).
}
\usage{
samp.size.MK(m = 1.2, w = 0.9, psi = 0.25, trunc = 3000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
Rough prior knowledge about m, the expected number of mutations per culture
}
  \item{w}{
Experimentalist's knowledge about w, relative fitness of mutants
}
  \item{psi}{
The psi score defined as the half width of a 95\% confidence interval for m divided by an anticipated magnitude of m
}
  \item{trunc}{
Number of terms kept in computing expected Fisher information
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A positive integer representing required sample size, i.e., the number of cultures
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Q. Zheng, Sample size determination for the fluctuation experiment, Mutation
Research/Fundamental and Molecular Mechanisms of Mutagenesis 795 (2017) 10-14.
}
\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
samp.size.LD.plating
}
\examples{
samp.size.MK(m=4,w=0.75,psi=0.25)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

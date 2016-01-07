\name{wierdl.est}
\alias{wierdl.est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wierdl method of accounting for variation in Nt
}
\description{
Nt is the (non-mutant) cell population size immediately prior to plating. In a fluctuation experiment 
comprising n tubes, the final cell population sizes are Nt(1), Nt(2), ..., Nt(n). The method
proposed by Wierdl et al. (1996) appears to be the first to account for variation in Nt. However,
this method requires the experimentalist to measure Nt for all tubes, which is a daunting
experimental challenge. The Wierdl method is a generalization of the method of the median devised by
Lea and Coulson (1949).
}
\usage{
wierdl.est(y, Nt, init.m = 0, tol = 1e-09, max.iter = 25,
   show.iter = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
A vector containing the mutant counts.
}
  \item{Nt}{
A vector containing the final cell population sizes.
}
  \item{init.m}{
An initial guess for m.
}
  \item{tol}{
Tolerance parameter to control numerical accuracy.
}
  \item{max.iter}{
A parameter to control the number of iterations in solving the Lea-Coulson equation.
}
  \item{show.iter}{
To show the iteration process.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A real number as an estimate of the mutation rate, not the usual mean number of mutations (m).
}
\references{
E.A.  Lea, C.A. Coulson, The distribution of the numbers of mutants in
bacterial populations. J Genetics 49 (1949) 264-285.

M. Wierdl, C.N. Green, A. Datta, S. Jinks-Robertson, T.D. Petes, Destabilization
of simple repetitive DNA sequences by transcription in yeast. Genetics 143 (1996) 713-721.
}
\author{
Qi Zheng <qzheng@sph.tamhsc.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
newton.B0
}
\examples{
Nt=rep(1.9e8, 30); wierdl.est(demerec.data,Nt)
}

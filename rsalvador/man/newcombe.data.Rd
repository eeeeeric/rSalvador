\name{newcombe.data}
\alias{newcombe.data}
\docType{data}
\title{
Classic Experimental Data of Howard B. Newcombe
}
\description{
The data frame contains 8 fluctuation experiments, labeled exptA, exptB, ... , exptH. Each
experiment consists of 25 cultures (0.2ml). Initial numbers of cells (N0) are:
(10,10,10,10,1e4,1e4,1e4,1e4). Terminal cell population sizes (Nt) are:
(3.1e8, 4.6e8, 2.5e8, 2.8e8, 4.2e8, 3.7e8, 3.2e8, 3.8e8).
}
% \usage{data(newcombe.data)}
 
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
HB Newcombe, Delayed phenotypic expression of spontaneous
mutations in Escherichia Coli, Genetics 33 (1948) 447-476.
}
\examples{
LD.p0.est(newcombe.data[[1]])
newton.LD(newcombe.data[[1]])
}
\keyword{datasets}

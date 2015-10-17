\name{boshoff.data}
\alias{boshoff.data}
\docType{data}
\title{Experimental Data by Boshoff et al.
}
\description{
Data from 6 fluctuation experiments by Boshoff et al. (2003). The final cell population sizes
Nt are: 7.2e8, 3.3e8, 1.02e9, 7.65e8, 2.1e9 and 1.9e9. According to the authors, some large
mutants counts are not very reliable. Data are from Supplementary Table of Boshoff (2003).
}

\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
HIM Boshoff, MB Reed, CE Barry III, V Mizrahi, DnaE2 polymerase contributes to in vivo 
survival and the emergence of drug resistance in Mycobacterium tuberculosis, Cell 113
(2003) 183-193.
}
\examples{
newton.LD(boshoff.data[[1]])/(7.2e8)
}
\keyword{datasets}

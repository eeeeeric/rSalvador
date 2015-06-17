\name{ford.data}
\alias{ford.data}
\docType{data}
\title{
Data of Mycobacterium tuberculosis Mutation by Fort et al.
}
\description{
One of several fluctuation experiments performed by Ford et al. (2011) to study Mycobacterium
tuberculosis mutation rates. Nt=1.08e9.
}

\details{
}

\references{
CB Ford, PL Lin, MR Chase, RR Shah, O Iartchouk, J Galagan, N. Mohaideen, TR Ioerger, 
JC Sacchettini, M Lipsitch, JL Flynn, SM Fortune, Use of whole genome sequencing to
estimate the mutation rate of Mycobacterium tuberculosis during latent infection,
Nature Genetics 43 (2011) 482-486. See supplementary inforaamtion. 
}
\examples{
newton.LD(ford.data)
confint.LD(ford.data)/(1.08e9)
}
\keyword{datasets}

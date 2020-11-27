\name{wh.data}
\alias{wh.data}
\docType{data}
\title{
Experimental Data of Werngren and Hoffner
}
\description{
There are 13 fluctuation experiments in this data object. From each
condensed culture of 2.5ml a sample of 1ml was plated, and hence the plating
efficiency is e=0.4. Terminal
cell population sizes are c(2.3, 0.5, 1.3, 1.0, 1.6, 1.0, 1.0, 0.8, 1.2, 0.6, 0.8, 0.9, 0.9)
times 1e8 per ml. Thus, for example, Nt=1.15e9 for the first experiment. See Zheng and
Werngren (2018) for further information about the Werngren-Hoffner data.
}
 
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
J Werngren, SE Hoffner, Drug-susceptible Mycobacterium tuberculosis Beijing genotype
does not develop mutation-conferred resistance to rifampin at an elevated rate,
Journal of Clinical Microbiology 41 (2003) 1520-1524.

Q Zheng, J Werngren, An unbiased attitude is vital to exploring the Beijing genotype of Mycobacterium
tuberculosis, Tuberculosis 111 (2018) 193-197.




}
\examples{
newton.LD.plating(wh.data[[1]],e = 0.4)
}
\keyword{datasets}

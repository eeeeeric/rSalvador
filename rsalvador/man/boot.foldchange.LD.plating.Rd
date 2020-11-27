\name{boot.foldchange.LD.plating}
\alias{boot.foldchange.LD.plating}
\title{
A Bootstrap Approach to Mutation Rate Fold Change Under an LD Model with Imperfect Plating Efficiency
}
\description{
A bootstrapping sampling approach is utilized to construct a confidence interval for mutation rate fold change under an
LD mutation model in which the plating efficiency is imperfect. The mutation rate fold change is defined as the mutation r
ate in the second experiment divided by the mutation rate in the second experiment. A similar approach was first taken
by Russell and March (2011)
}
\usage{
boot.foldchange.LD.plating(y1, y2, N1, N2, e1 = 0.1, e2 = 0.1,
         alpha = 0.05, init.m1 = 0, init.m2 = 0, na.ok = FALSE, nboot = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y1}{
Mutant count data from the first experiment.
}
  \item{y2}{
Mutant count data from the second experiment.
}
  \item{N1}{
Average final cell population size of the first experiment.
}
  \item{N2}{
Average final cell population size of the second experiment.
}
  \item{e1}{
Plating efficiency in the first experiment, must be smaller than 1.0.
}
  \item{e2}{
Plating efficiency in the second experiment, must be smaller than 1.0.
}
  \item{alpha}{
The confidence interval will have an approximate confidence coefficient of 1-alpha.
}
  \item{init.m1}{
An initial guess of the expected number of mutations per culture in the first experiment.
}
  \item{init.m2}{
An initial guess of the expected number of mutations per culture in the second experiment.
}

  \item{na.ok}{
Not all bootstrap samples will allow the computational process for mutation rates to converge, often due to
poorly chosen starting values. If na.ok is set to TRUE, such recalcitrant samples will be discarded before
calculating a Bayesian confident interval. Otherwise, it stops further data processing and simply returns
an "NA" as output.
}



  \item{nboot}{
Number of bootstrap samples to be drawn.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of the form (lower confidence limit, median, upper confidence limit)
}
\references{
M.S. Russell, J. C. March, Bootstrap estimation of confidence intervals on
mutation rate ratios. Environmental and Molecular Mutagenesis 52 (2011) 385-396.
}
\author{
Qi Zheng <qzheng@tamu.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
boot.foldchange.MK
}
\examples{
y1 = wh.data[[1]]; y2 = wh.data[[2]]; N1 = 2.3e8; N2 = 0.5e8;
boot.foldchange.LD.plating(y1, y2, N1, N2, e1 = 0.4, e2 = 0.4)
}

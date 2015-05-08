\name{log.likelihood.LD}
\alias{log.likelihood.LD}

\title{Computing the Log-Likelihood Function Under the Luria-Delbruck Model}

\description{This function computes the log-likelihood function for data under the Lea-Coulson
formulation of the Luria-Delbruck mutation model.}

\usage{log.likelihood.LD(data, m, phi = 1)}

\arguments{
  \item{data}{A vector of data.}
  \item{m}{The parameter m is regarded as the argument of the log-likelihood function. 
}
  \item{phi}{This parameter is defined as phi=1-N0/Nt.}
} %% end of arguments

\author{Qi Zheng <qzheng@sph.tamhsc.edu>}

\seealso{log.likelihood.LD.plating}

\examples{log.likelihood.LD(demerec.data,m=8.5)}

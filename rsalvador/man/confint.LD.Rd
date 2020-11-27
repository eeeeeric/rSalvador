\name{confint.LD}
\alias{confint.LD}

\title{Likelihood Ratio-Based Confidence Interval under the Luria-Delbruck Model
}
\description{
This function employs an algorithm described in Zheng (2005) to compute a likelihood
ratio-based confidence interval for m, the expected number of mutations per test tube.
}
\usage{
confint.LD(data, alpha = 0.05, phi = 1, tol = 1e-9, init.m = 0,
 init.lower = 0, init.upper = 0, max.iter = 30, show.iter = FALSE)}

\arguments{
  \item{data}{a vector of numbers of mutant colonies

}
  \item{alpha}{significance level

}
  \item{phi}{this parameter is defined by phi=1-N0/Nt, is usually assumed to be 1.0.
}
  \item{tol}{tolerance parameter 
}
  \item{init.m}{an initial guess for m, rarely needed

}
  \item{init.lower}{an initial guess for lower limit of m, rarely needed
%%     ~~Describe \code{init.lower} here~~
}
  \item{init.upper}{an initial guess for upper limit of m, rarely needed
}
  \item{max.iter}{maximum number of iterations
}
  \item{show.iter}{ to show iteration progress
%%     ~~Describe \code{show.iter} here~~
}
}
\details{
The algorithm is described in Zheng (2005).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Q. Zheng, (2005). New algorithms for Luria-Delbruck fluctuation analysis,
Mathematical Biosciences 196:198-214 (2005).
}
\author{Qi Zheng <qzheng@sph.tamhsc.edu> }
\note{
}

\seealso{
confint.LD.plating
}
\examples{ confint.LD(demerec.data,show.iter=TRUE)
}

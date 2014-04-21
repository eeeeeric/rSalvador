\name{import.text.data}
\alias{import.text.data}

\title{Importing Data from a Text File}

\description{
The text data file to be imported can have several header lines at the beginning, followed by
a single-column of experimental data. The parameter jump specifies the number of header lines.
}
\usage{import.text.data(filename, jump = 1, col = 1)}


\arguments{
\item{filename}{Experimental data should be listed in a single column, but the first few lines
can be the experimenter's note.}

\item{jump}{Number of lines to skip.}

\item{col}{Which column contains data.}
} %% end of argument

\value{A vector of non-negative integers imported form the text file.}

\author{Qi Zheng <qzheng@sph.tamhsc.edu>}


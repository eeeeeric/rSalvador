\name{import.excel.data}
\alias{import.excel.data}

\title{ Importing Data from an Excel File}

\description{
Biologists often keep their experimental results in Excel files. In addition to the numbers of
mutant cells, such files usually contain other experimental details. rSalvador requires
excel files to satisfy two conditions. First, the first row of the file should be reserved as
a header, and second, the first non-blank column should be reserved to keep the numbers of
mutant cells that rSalvador intends to import. The parameter col (which has a default value
of 1) can be used to choose which column of the file is to be imported.}

\usage{
import.excel.data(filename, col = 1)
}

\arguments{
  \item{filename}{ the data file to be imported.
}
  \item{col}{the desired column to be imported.
}
}

\value{A vector of non-negative integers.}

\author{Qi Zheng <qzheng@sph.tamhsc.edu>}

\seealso{import.text.data}


\name{export.text.data}
\alias{export.text.data}

\title{Exporting Data to an External Text File}

\description{It is sometimes desirable to save experimental data already in rSalvador
to an external text file. For example, data imported by import.excel.data may later
be exported as a text file so that other programs can use the data.
}
\usage{
export.text.data(filename, data)
}

\arguments{
  \item{filename}{Name of the text file to be written.}
  \item{data}{A vector of experimental data to be exported.}
}

\author{Qi Zheng <qzheng@sprh.tamhsc.edu>}

\seealso{import.text.data}

\examples{export.text.data('testing.txt',demerec.data)}


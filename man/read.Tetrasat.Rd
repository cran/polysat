\name{read.Tetrasat}
\alias{read.Tetrasat}
\title{
Read Data from a TETRASAT Input File
}
\description{
Given a file containing genotypes in the TETRASAT format,
\code{read.Tetrasat} produces a vector containing the population
identity of each individual as well as a two-dimensional list of vectors
containing the genotypes.
}
\usage{
read.Tetrasat(infile, missing = -9)
}
\arguments{
  \item{infile}{
A character string of the file path to be read.
}
  \item{missing}{
A numerical value to be used to indicate missing data in the object
produced by the function.
}
}
\details{
\code{read.Tetrasat} reads text files that are in the exact format
specified by the software TETRASAT  and TETRA (see references for more
information).  This is similar to the file format for GenePop but allows
for up to four alleles per locus.  All alleles must be coded by two
digits.  Another diference between the TETRASAT and GenePop formats is
that in TETRASAT the sample name and genotypes are not separated by a
comma, because the columns of data have fixed widths.

Since TETRASAT files also contain information about which samples belong
to which populations, a vector is also produced containing this
information.  The elements of the vector are numbers representing the
populations, and the names of the elements are the sample names.  This
can be used as the popinfo argument for \code{estimate.freq}, for
example, or converted to a column in an array for use in the
\code{extracols} argument of \code{write.Structure}.
}
\value{
\item{PopData}{An integer vector containing the population ID for each
  sample.  Sample names are used as the names of the elements.}
\item{Genotypes}{A two-dimensional list of integer vectors, in the
  standard polysat genotype format.  Samples are the first dimension,
  and loci are the second dimension.  All white space is stripped from
  sample names.  Each vector contains all unique alleles for that sample
  at that locus.}
}
\references{
  \url{http://markwith.freehomepage.com/tetrasat.html}

  Markwith, S. H., Stewart, D. J. and Dyer, J. L.  (2006)
  TETRASAT: a program for the population analysis of
  allotetraploid microsatellite data. \emph{Molecular Ecology Notes}
  \bold{6}, 586-589.

  \url{http://ecology.bnu.edu.cn/zhangdy/TETRA/TETRA.htm}

  Liao, W. J., Zhu, B. R., Zeng, Y. F. and Zhang, D. Y.  (2008)
  TETRA: an improved program for population genetic analysis of
  allotetraploid microsatellite data.  \emph{Molecular Ecology Resources}
  \bold{8}, 1260-1262.
}
\author{
Lindsay V. Clark
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{read.GeneMapper}}, \code{\link{write.Tetrasat}},
\code{\link{read.ATetra}}, \code{\link{read.GenoDive}},
\code{\link{read.Structure}}, \code{\link{dominant.to.codominant}},
\code{\link{read.SPAGeDi}}
}
\examples{
\dontrun{
# example from the Tetrasat website
mydata <- read.Tetrasat("http://markwith.freehomepage.com/sample.txt")
mydata$PopData
mydata$Genotypes["BCRHE3",]
mydata$Genotypes[["BR6","B1_Gtype"]]
}

# example with defined data:
cat("Sample Data", "A1_Gtype", "A10_Gtype", "B1_Gtype", "D7_Gtype",
"D9_Gtype", "D12_Gtype", "Pop",
"BCRHE 1             0406     04040404 0208     02020202 03030303 0710",
"BCRHE 10            0406     04040404 07070707 02020202 0304     0710",
"BCRHE 2             04040404 04040404 0708     02020202 010305   0710",
"BCRHE 3             04040404 04040404 02020202 0203     03030303 0809",
"BCRHE 4             04040404 04040404 0608     0203     03030303 070910",
"BCRHE 5             04040404 04040404 0208     02020202 03030303 050710",
"BCRHE 6             0304     04040404 0207     02020202 03030303 07070707",
"BCRHE 7             0406     04040404 0708     02020202 03030303 07070707",
"BCRHE 8             0304     04040404 0203     0203     03030303 0709",
"BCRHE 9             0406     04040404 0708     02020202 03030303 0710",
"Pop",
"BR 1                0406     04040404 05050505 02020202 03030303 1012",
"BR 10               030406   04040404 0607     02020202 03030303 1011",
"BR 2                030406   04040404 07070707 02020202 03030303 09090909",
"BR 3                010304   04040404 07070707 02020202 03030303 09090909",
"BR 4                030406   04040404 07070707 0203     03030303 10101010",
"BR 5                030406   04040404 07070707 02020202 03030303 10101010",
"BR 6                0406     04040404 0507     0203     03030303 10101010",
"BR 7                0304     04040404 0809     02020202 03030303 070910",
"BR 8                030406   04040404 07070707 02020202 03030303 070910",
"BR 9                0406     04040404 07070707 02020202 03030303 07070707",
sep="\n", file="TetrasatExample.txt")
mydata2 <- read.Tetrasat("TetrasatExample.txt")
}
\keyword{file}

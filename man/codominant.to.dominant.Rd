\name{codominant.to.dominant}
\alias{codominant.to.dominant}
\title{
Convert a Genotype Object to Binary Data
}
\description{
Given a genotype object in the form of a two-dimensional list of
vectors, \code{codomiant.to.dominant} produces a matrix of binary data
indicating the presence and absence of alleles in samples.
}
\usage{
codominant.to.dominant(gendata, makecolinfo = FALSE, allelepresent = 1,
alleleabsent = 0, missingin = -9, missingout = -9,
loci = dimnames(gendata)[[2]], samples = dimnames(gendata)[[1]])
}
\arguments{
  \item{gendata}{
    Genotypes.  A two dimensional list of vectors as produced by other
    functions in polysat.  Samples are represented in the first
    dimension of the list and loci in the second dimension.  Each vector
    contains all unique alleles for a given sample at a given locus.
}
  \item{makecolinfo}{
Boolean.  If TRUE, \code{codominant.to.dominant} also produces a data frame
indicating which loci and alleles are represented in which columns.
}
  \item{allelepresent}{
The symbol to be used in the output to indicate that a sample has an allele.
}
  \item{alleleabsent}{
The symbol to be used in the output to indicate that a sample does not
have an allele.
}
  \item{missingin}{
The symbol used to indicate missing data in \code{gendata}.
}
  \item{missingout}{
The symbol used to indicate missing data in the output.
}
  \item{loci}{
A character vector that is a subset of \code{dimnames(gendata)[[2]]},
indicating which loci to process.
}
  \item{samples}{
A character vector that is a subset of \code{dimnames(gendata)[[1]]},
indicating which samples to process.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  If \code{makecolinfo=FALSE}, a matrix is returned containing the binary
  genotype data.  The row names of the matrix are the sample names.  The
  column names of the matrix are locus and allele names concatenated
  together, separated by a period.  The symbols \code{allelepresent},
  \code{alleleabsent}, and \code{missingout} (by default 1, 0, and -9) are
  used in the matrix to indicate the presence or absence of alleles,
  or missing data, respectively.  Loci are in the same order as in the
  input, and alleles are arranged in numerical order.

  If \code{makecolinfo=TRUE}, the matrix described above is returned as the
  first element of a list.  The second element is a data frame with loci
  in the first column and alleles in the second column, with data in the
  same order as the columns in the matrix.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Lindsay V. Clark
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dominant.to.codominant}}, \code{\link{write.table}}
}
\examples{
# Create a data set to convert
mygendata <-
  array(list(c(124,128,138),c(122,130,140,142),c(122,132,136),c(122,134,140),
             c(203,212,218),c(197,206,221),c(215),c(200,218),
             c(140,144,148,150),c(-9),c(146,150),c(152,154,158))
        , dim=c(4,3), dimnames=list(c("ind1","ind2","ind3","ind4"),
  c("locus1","locus2","locus3")))

# Convert the data
mybinary <- codominant.to.dominant(mygendata, makecolinfo=TRUE)

# Inspect the results
mybinary[[1]]
mybinary[[2]]

}
\keyword{manip}
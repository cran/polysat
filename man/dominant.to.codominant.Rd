\name{dominant.to.codominant}
\alias{dominant.to.codominant}
\title{Convert Genotypes from Dominant Format to Codominant Format}
\description{
  This function takes an array or matrix of genotypes where each allele
  is represented in one column and symbols indicate the presence or
  absence of that allele in a sample.  It produces a two-dimensional
  list of vectors representing genotypes, indexed by sample and locus.
}
\usage{dominant.to.codominant(domdata, colinfo = NULL,
samples = dimnames(domdata)[[1]], missing = -9, allelepresent = 1, split
= ".")
}
\arguments{
  \item{domdata}{A two-dimensional array or matrix, in which samples are
    represented in the first dimension (and named accordingly) and
    alleles are represented in the second dimension.  The
    symbol specified by \code{allelepresent} indicates that a
    sample has a particular allele, and any other symbol indicates that it
    does not.  If \code{colinfo} is not provided by the user, the second
    dimension names of the array should be the locus name and allele number
  connected by a period or by another character as
  specified in \code{split}.}
\item{colinfo}{A data frame, indexed by column number from
  \code{domdata}, containing locus names as the first column and
  allele numbers as the second column.}
\item{samples}{A character vector containing the names of samples to be
  converted, if only a subset of samples in \code{domdata} are to be used.}
\item{missing}{The symbol to use to represent missing data in the
  output.}
\item{allelepresent}{The symbol used in \code{domdata} to indicate that a
  particular sample has a particular allele.}
\item{split}{If \code{colinfo=NULL}, the character used to separate the locus
  name and allele number in the column names of \code{domdata}.}
}
\value{
  A two-dimensional list of integer vectors, in the standard polysat
  genotype format.  Samples are represented in the first
  dimension and loci in the second dimension, and both are named
  accordingly.  Each vector contains all unique alleles for a given
  sample at a given locus.
}
\details{
  Because allele copy number is often unknown, many researchers who work
  with microsatellites in polyploids record genotype data in a dominant
  format, such as a matrix of 1's and 0's to represent the presence and
  absence of peaks as is done with AFLPs.  \code{dominant.to.codominant} is
  written to convert that data back to a semi-codominant format so that
  other analyses or data conversion can be performed.

  The default symbol to indicate the presence of an allele is 1, but
  this can be set to any other symbol using the \code{allelepresent}
  argument.  It does not matter which symbols are used to indicate that
  an allele is absent or that there is missing data.  If
  \code{dominant.to.codominant} does not find any alleles present for a
  given sample and locus, it fills in a missing data symbol in that
  position in the two-dimensional genotype list.

  This function does not read or write files.  Since the user would
  already have dominant data in an array-like format in a spreadsheet or
  text document, it should be easily read by \code{read.table} and
  converted to a matrix by \code{as.matrix}.

  There are two options for indicating which locus and allele is
  represented by each column:

  1) These can be specified in the second
  dimension names of the array or matrix.  The name of each column
  should be a concatenation of the locus name followed by the allele
  number, and these should be separated by a period or other character
  as specified in \code{split} (e.g. \dQuote{locus1.204}).  Note that with
  \code{check.names=TRUE}, \code{read.table} will convert a lot of symbols (like
  hyphens or spaces) to periods.  It is probably a good idea to inspect
  the column names of \code{domdata} before setting \code{split}.

  2) Create a data frame containing locus and allele information.  The
  rows should be in the same order as the columns of \code{domdata}.  The
  first vector in the data frame should contain the locus names, and the
  second vector in the data frame should contain the numerical alleles.
  Use this data frame as \code{colinfo}.
}
\references{
}
\seealso{
  \code{\link{codominant.to.dominant}}, \code{\link{read.table}}, \code{\link{as.matrix}}
}
\examples{
# Create a matrix of dominant data (usually read from a file instead)
mysamples <- c("ind1","ind2","ind3")
myalleles <- c("loc1.100","loc1.102","loc1.104","loc1.106",
"loc2.141","loc2.144","loc2.147","loc2.150")
mydomdata <- matrix(nrow = length(mysamples), ncol = length(myalleles),
                    dimnames = list(mysamples, myalleles))
mydomdata["ind1",] <- c(1,1,1,0,0,1,1,0)
mydomdata["ind2",] <- c(1,0,0,1,0,0,1,1)
mydomdata["ind3",] <- c(-9,-9,-9,-9,1,1,0,1)

# inspect the matrix
mydomdata

# convert to codominant data
mycodomdata <- dominant.to.codominant(mydomdata)

# view the list created
mycodomdata
# view genotypes by individual
mycodomdata["ind1",]
mycodomdata["ind2",]
mycodomdata["ind3",]

# Alternately, use a matrix without alleles labeled in the colunn names
dimnames(mydomdata)[[2]] <- NULL
mydomdata

# Make a data frame for a locus and allele index
# (Under normal circumstances you would read this from a file)
laindex <- data.frame(Loci = c(rep("loc1",4), rep("loc2",4)),
Alleles = c(100, 102, 104, 106, 141, 144, 147, 150))
laindex

# convert to codominant data
mycodomdata2 <- dominant.to.codominant(mydomdata, colinfo=laindex)
# look at the results
mycodomdata2["ind1",]
# etc.
}
\author{Lindsay V. Clark}
\keyword{manip}
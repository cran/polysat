\name{write.GeneMapper}
\alias{write.GeneMapper}
\title{
Write Genotypes to a Table Similarly to ABI GeneMapper
}
\description{
Given a genotype object in the standard polysat format,
\code{write.GeneMapper} writes a text file of a table containing columns
for sample name, locus, and alleles.
}
\usage{
write.GeneMapper(gendata, file = "", samples = dimnames(gendata)[[1]],
loci = dimnames(gendata)[[2]])
}
\arguments{
  \item{gendata}{
Genotypes in the standard polysat format: a two-dimensional list of
vectors.  Samples are represented in the first dimension of the list and
loci in the second dimension of the list, and both dimensions are named
accordingly.  Each vector contains all unique alleles for a given sample
and locus.
}
  \item{file}{
Character string.  The path to which to write the file.
}
  \item{samples}{
Character vector.  Samples to write to the file.  This should be a
subset of \code{dimnames(gendata)[[1]]}.
}
  \item{loci}{
Character vector.  Loci to write to the file.  This should be a subset
of \code{dimnames(gendata)[[2]]}.
}
}
\details{
  Although I do not know of any population genetic software that will
  read this data format directly, the ABI GeneMapper Genotypes Table
  format is a convenient way for the user to store microsatellite
  genotype data and view it in a text editor or spreadsheet software.
  Each row contains the sample name, locus name, and alleles separated
  by tabs.

  The number of allele columns needed is automatically detected by the
  function by using \code{estimate.ploidy} to find the maximum number of
  alleles per genotype.

  \code{write.GeneMapper} handles missing data in a very simple way, in
  that it writes the missing data symbol directly to the table as though
  it were an allele.  If you want missing data to be represented
  differently in the table, you can open it in spreadsheet software and
  use find/replace or conditional formatting to locate missing data.

  The file that is produced can be read back into R directly by
  \code{read.GeneMapper}, and therfore may be a convenient way to back up
  genotype data for future analysis and manipulation in polysat.
  (\code{save} can also be used to back up an R object more directly.)  If
  you need to edit the genotypes and would rather do it using
  spreadsheet software than in polysat
  (e.g. \code{mygendata[["ind1","loc3"]]<-c(122,126,134)}) writing the file
  in GeneMapper format, editing it, then reading it back in is an option.
}
\value{
A file is written but no value is returned.
}
\references{
  \url{http://www.appliedbiosystems.com/genemapper}
}
\author{
Lindsay V. Clark
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{read.GeneMapper}}, \code{\link{write.Structure}},
\code{\link{write.GenoDive}}, \code{\link{write.Tetrasat}},
\code{\link{write.ATetra}}, \code{\link{codominant.to.dominant}},
\code{\link{write.SPAGeDi}}
}
\examples{
# create a genotype object (usually done by reading a file)
mysamples <- c("ind1", "ind2", "ind3", "ind4")
myloci <- c("loc1", "loc2")
mygendata <- array(list(-9), dim=c(4,2), dimnames=list(mysamples, myloci))
mygendata[,"loc1"] <- list(c(202,204), c(204), c(200,206,208,212),
                           c(198,204,208))
mygendata[,"loc2"] <- list(c(78,81,84), c(75,90,93,96,99), c(87), c(-9))

# write a GeneMapper file
write.GeneMapper(mygendata, "exampleGMoutput.txt")

# view the file with read.table
read.table("exampleGMoutput.txt", sep="\t", header=TRUE)
}
\keyword{ file }
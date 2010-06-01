\name{write.SPAGeDi}
\alias{write.SPAGeDi}
\title{
Write Genotypes in SPAGeDi Format
}
\description{
\code{write.SPAGeDi} takes a genotype object in the standard polysat
format and creates a file that can be read by the software SPAGeDi.  The
user controls how the genotypes are formatted, and can specify the
ploidy, population, and spatial coordinates of each sample.
}
\usage{
write.SPAGeDi(gendata, samples = dimnames(gendata)[[1]],
              loci = dimnames(gendata)[[2]],
              indploidies = rep(4, length(samples)),
              popinfo = rep(1, length(samples)), allelesep = "/",
              digits = 2, file = "",
              spatcoord = data.frame(X = rep(1, length(samples)),
                                     Y = rep(1, length(samples)),
                                     row.names = samples),
              usatnts = rep(2, length(loci)), missing = -9)
}
\arguments{
  \item{gendata}{
A genotype object in the standard polysat format.  A two dimensional
list of integer vectors, where samples are represented and named in the first
dimension, and loci are represented and named in the second dimension.
Each vector contains all unique allele for a given sample and locus.
}
  \item{samples}{
Character vector.  Samples to write to the file.  Must be a subset of
dimnames(gendata)[[1]].
}
  \item{loci}{
Character vector.  Loci to write to the file.  Must be a subset of
dimnames(gendata)[[2]].
}
  \item{indploidies}{
Integer vector.  Ploidy of each sample.  This can either be named by
sample, or be in the same order as \code{samples}.
}
  \item{popinfo}{
Vector.  Population identity (or category) of each individual.  This can
either be named by sample, or be in the same order as \code{samples}.
}
  \item{allelesep}{
The character that will be used to separate alleles within a genotype.
If each allele should instead be a fixed number of digits, with no
characters to delimit alleles, set \code{allelesep = ""}.
}
  \item{digits}{
Integer.  The number of digits used to represent each allele.
}
  \item{file}{
Character string.  The file path to write to.
}
  \item{spatcoord}{
Data frame.  Spatial coordinates of each sample.  Column names are used
for column names in the file.  Row names indicate sample, or if absent
it is assumed that the rows are in the same order as \code{samples}.
}
  \item{usatnts}{
Integer vector.  Repeat length of each locus (eg. 2 to indicate
dinucleotide repeats, or 3 to indicate trinucleotide repeats).  If the
alleles in \code{gendata} are already in repeat lengths rather than
nucleotides, the value should be 1.  The vector should be named by loci,
or else is assumed to be in the same order as \code{loci}.
}
  \item{missing}{
The symbol used in \code{gendata} to indicate missing data.
}
}
\details{
\code{popinfo}, \code{indploidies}, \code{spatcoord}, and \code{usatnts}
can have more samples or loci than are to be written to the file, as
long as they are named.

The first line of the file contains the number of individuals, number of
categories, number of spatial coordinates, number of loci, number of
digits for coding alleles, and maximum ploidy, and is generated
automatically from the data provided.

The function does not write distance intervals to the file, but instead
writes \code{0} to the second line.

All alleles for a given locus are divided by the \code{usatnts} value
for that locus, after all missing data symbols have been replaced with
zeros.  If necessary, a multiple of 10 is subtracted from all
alleles at a locus in order to get the alleles down to the right number
of digits.

If a genotype has fewer alleles than the \code{indploidies} value for
that sample, zeros are added up to the ploidy.  If the genotype has more
alleles than the ploidy, a random subset of alleles is used and a
warning is printed.  If the genotype has only one allele (is fully
heterozygous), then that allele is replicated to the ploidy of the
individual.  Genotypes are then concatenated into strings,
delimited by \code{allelesep}.  If \code{allelesep=""}, leading zeros
are first added to alleles as necessary to make them the right number of
digits.
}
\value{
A file is written but no value is returned.
}
\references{
  \url{http://ebe.ulb.ac.be/ebe/Software_files/manual_SPAGeDi_1-3.pdf}

  Hardy, O. J. and Vekemans, X. (2002) SPAGeDi: a versatile computer
  program to analyse spatial genetic structure at the individual or
  population levels.  \emph{Molecular Ecology Notes} \bold{2}, 618-620.
}
\author{
Lindsay V. Clark
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{read.SPAGeDi}}, \code{\link{write.GenoDive}},
\code{\link{write.Structure}}, \code{\link{write.GeneMapper}},
\code{\link{write.ATetra}}, \code{\link{write.Tetrasat}},
\code{\link{codominant.to.dominant}}
}
\examples{
# set up data to write (usually read from a file)
mygendata <- array(list(-9), dim=c(4,2),
                   dimnames=list(c("ind1","ind2","ind3","ind4"),
                                 c("loc1", "loc2")))
mygendata["ind1",] <- list(c(102,106,108),c(207,210))
mygendata["ind2",] <- list(c(104),c(204,210))
mygendata["ind3",] <- list(c(100,102,108),c(201,213))
mygendata["ind4",] <- list(c(102,112),c(-9))
myploidies <- c(3,2,2,2)
names(myploidies) <- c("ind1","ind2","ind3","ind4")
myusatnts <- c(2,3)
names(myusatnts) <- c("loc1","loc2")
myspatcoord <- data.frame(X=c(27,29,24,30), Y=c(44,41,45,46),
                          row.names=c("ind1","ind2","ind3","ind4"))

# write a file
write.SPAGeDi(mygendata, indploidies = myploidies, usatnts = myusatnts,
              spatcoord = myspatcoord, file="SpagOutExample.txt")
}
\keyword{ file }


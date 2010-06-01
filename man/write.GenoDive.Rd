\name{write.GenoDive}
\alias{write.GenoDive}
\title{
Write a File in GenoDive Format
}
\description{
\code{write.GenoDive} takes the standard genotype object used by polysat
and creates a file formatted for the software GenoDive.
}
\usage{
write.GenoDive(gendata, popnames = "onebigpop",
commentline = "file description goes here", digits = 2, file = "",
samples = dimnames(gendata)[[1]], loci = dimnames(gendata)[[2]],
popinfo = rep(1, times = length(samples)),
usatnts = rep(2, times = length(loci)), missing=-9)
}
\arguments{
  \item{gendata}{
A two-dimensional list of vectors representing genotypes, in the
standard polysat format.  Samples are the first dimension of the list
and loci are the second dimension.  Dimensions are named accordingly.
Each vector is numerical and contains all unique alleles for a given
sample and locus.
}
  \item{popnames}{
A character vector of population names, ordered by the population
numbers used in \code{popinfo}.
}
  \item{commentline}{
Character string.  The first line of the file, containing a description
of the data.
}
  \item{digits}{
Integer.  How many digits to use to represent each allele (usually 2 or 3).
}
  \item{file}{
A character string of the file path to which to write.
}
  \item{samples}{
    A character vector of samples to include in the file.  A subset of
    \code{dimnames(gendata)[[1]]}.
}
  \item{loci}{
    A character vector of loci to include in the file.  A subset of
    \code{dimnames(gendata)[[2]]}.
}
  \item{popinfo}{
A numeric or integer vector containing the population number of each
individual.  \code{names(popinfo)} should include all of \code{samples}.
If the vector is unnamed, it is assumed to be in the same order as
\code{samples}.
}
  \item{usatnts}{
A numeric or integer vector, with one element for each locus.
\code{names(usatnts)} should include all of \code{loci}, or if the vector is
unnamed it is assumed to be in the same order as \code{loci}.  Each element
represents the length of the nucleotide repeat for the locus, such as 2
for dinucleotide repeats or 3 for trinucleotide repeats.  If the alleles
are already expressed in \code{gendata} as repeat number rather than
fragment length in nucleotides, 1 should be the entry for that locus.
}
\item{missing}{
  The symbol used to indicate missing data in the input.
  }
}
\details{
The number of individuals, number of populations, number of loci, and
maximum ploidy of the sample are calculated automatically and entered in
the second line of the file.  If the maximum ploidy needs to be reduced
by random removal of alleles, it is possible to do this in the software
GenoDive after importing the data.

Population names and population identies of individuals can optionally
be entered as arguments in the file.  If you do not already have this
data readily available in R, you may prefer to enter it manually using a
text editor or spreadsheet software after the file is written by
\code{write.GenoDive}.

Several steps happen in order to convert alleles to the right format.
First, all instances of the missing data symbol are replaced with zero.
Alleles are then divided by the numbers provided in \code{usatnts} (and
rounded down if necessary) in order to
convert them from nucleotides to repeat numbers.  If the alleles are
still too big to be represented by the number of digits specified,
\code{write.GenoDive} repeatedly subtracts a number (10^(digits-1), so 10
if \code{digits=2}) from all alleles at a locus until the alleles are small
enough.  Alleles are then converted to characters, and a leading zero is
added to an allele if it does not have enough digits.  These alleles are
concatenated at each locus so that each sample*locus genotype is an
uninterrupted string of numbers.
}
\value{
A file is written but no value is returned.
}
\references{
  Meirmans, P. G. and Van Tienderen P. H. (2004) GENOTYPE and GENODIVE:
  two programs for the analysis of genetic diversity of asexual
  organisms \emph{Molecular Ecology Notes} \bold{4}, 792-794.

  \url{http://www.bentleydrummer.nl/software/software/GenoDive.html}
}
\author{
Lindsay V. Clark
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{read.GenoDive}}, \code{\link{write.Structure}},
\code{\link{write.ATetra}}, \code{\link{write.Tetrasat}},
\code{\link{write.GeneMapper}}, \code{\link{codominant.to.dominant}},
\code{\link{write.SPAGeDi}}
}
\examples{
# set up the genotype object (usually done by reading a file)
mysamples <- c("Mal", "Inara", "Kaylee", "Simon", "River", "Zoe",
               "Wash", "Jayne", "Book")
myloci <- c("loc1", "loc2")
mygendata <- array(list(-9), dim=c(9,2), dimnames=list(mysamples,
                                                       myloci))
mygendata[,"loc1"] <- list(c(304,306), c(302,310), c(306), c(312,314),
c(312,314), c(308,310), c(312), c(302,308,310), c(-9))
mygendata[,"loc2"] <- list(c(118,133), c(121,130), c(122,139),
c(124,133), c(118,124), c(121,127), c(124,136), c(124,127,136),
c(121,130))

# make the vector of nucleotide repeat lengths
myusatnts <- c(2,3)
names(myusatnts) <- myloci

# set up population info (some or all of this might also be read from a file)
mypopnames <- c("Core","Outer Rim")
mypopinfo <- c(rep(1,4) ,rep(2,6))
names(mypopinfo) <- c("Simon", "River", "Inara", "Book", "Kaylee",
                      "Mal", "Wash", "Zoe", "Badger", "Jayne")

# write files (use file="" to display in console instead)
write.GenoDive(mygendata, usatnts=myusatnts, popnames=mypopnames,
  popinfo=mypopinfo, digits=2, commentline="Serenity crew",
  file="testGenoDive2.txt")
write.GenoDive(mygendata, usatnts=myusatnts, popnames=mypopnames,
  popinfo=mypopinfo, digits=3, commentline="Serenity crew",
  file="testGenoDive3.txt")
}
\keyword{file}
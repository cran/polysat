\name{read.GenoDive}
\alias{read.GenoDive}
\title{
Import Genotype Data from GenoDive File
}
\description{
\code{read.GenoDive} takes a text file in the format for the software
GenoDive and produces a vector indicating which samples are in which
populations and a two-dimensional list of vectors containing the genotypes.
}
\usage{
read.GenoDive(infile, missing = -9)
}
\arguments{
  \item{infile}{
A character string.  The path to the file to be read.
}
  \item{missing}{
The symbol used to represent missing data in the output.
}
}
\details{
GenoDive is a Mac-only program for population genetic analysis that
allows for polyploid data.  \code{read.GenoDive} imports data from text
files formatted for this program.

The first line of the file is a comment line, ignored by the function.
On the second line, separated by tabs, are the number of individuals,
number of populations, number of loci, maximum ploidy
(ignored), and number of digits used to code alleles.

The following lines contain the names of populations, which are ignored
by the function although the number of lines should be the same as the
number of populations as specified in the first line.  After that is a
header line for the genotype data.  This line contains, separated by
tabs, column headers for populations, clones (optional), and
individuals, followed by the name of each locus.  The loci names for the
genotype object are derived from this line.

Each individual is on one line following the genotype header line.
Separated by tabs are the population number, the clone number
(optional), the individual name (used as the sample name in the output)
and the genotypes at each locus.  Alleles at one locus are concatenated
together in one string without any characters to separate them.  Each
allele must have the same number of digits, although leading zeros can
be omitted.

If the only alleles listed for a particular individual and locus are
zeros, this is interpreted by \code{read.GenoDive} as missing data, and
\code{missing} (-9 by default) is recorded in the list.
GenoDive allows for a genotype to be partially missing but polysat does
not; therefore, if an allele is coded as zero but other alleles are
recorded for that sample and locus, the output genotype will just
contain the alleles that are present, with the zeros thrown out.
}
\value{
  \item{PopData}{A numeric vector containing the population number of
    each individual.  The names of the vector are the names of the individuals.}
  \item{Genotypes}{A two-dimensional list of integer vectors, containing
  genotypes in the same format that is produced and read by other
  polysat functions.  The first dimension is indexed by sample and the
  second dimension is indexed by locus.  Each vector contains all unique
  alleles for a given sample and locus.}
}
\references{
  Meirmans, P. G. and Van Tienderen, P. H. (2004) GENOTYPE and GENODIVE:
  two programs for the analysis of genetic diversity of asexual
  organisms. \emph{Molecular Ecology Notes} \bold{4}, 792-794.

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
\code{\link{read.GeneMapper}}, \code{\link{write.GenoDive}},
\code{\link{read.Tetrasat}}, \code{\link{read.ATetra}},
\code{\link{read.Structure}}, \code{\link{dominant.to.codominant}},
\code{\link{read.SPAGeDi}}
}
\examples{
# create data file (normally done in a text editor or spreadsheet software)
cat(c("example comment line", "5\t2\t2\t3\t2", "pop1", "pop2",
      "pop\tind\tloc1\tloc2", "1\tJohn\t102\t1214",
      "1\tPaul\t202\t0", "2\tGeorge\t101\t121213",
      "2\tRingo\t10304\t131414","1\tYoko\t10303\t120014"),
    file = "genodiveExample.txt", sep = "\n")

# import file data
exampledata <- read.GenoDive("genodiveExample.txt")

# view data
exampledata$PopData
exampledata$Genotypes[,"loc1"]
exampledata$Genotypes[,"loc2"]
}
\keyword{file}

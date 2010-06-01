\name{read.SPAGeDi}
\alias{read.SPAGeDi}
\title{
Read Genotypes in SPAGeDi format
}
\description{
\code{read.SPAGeDi} can read a text file formatted for the SPAGeDi
software and return a genotype object in the standard polysat format,
as well as optionally returning a vector of individual ploidies and/or a
data frame of categories and spatial coordinates.
}
\usage{
read.SPAGeDi(infile, allelesep = "/", returncatspatcoord = FALSE,
             returnploidies = FALSE, missing = -9)
}
\arguments{
  \item{infile}{
Character string.  Path of the file to read.
}
  \item{allelesep}{
The character that is used to delimit alleles within genotypes, or
\code{""} if alleles have a fixed number of digits and are not delimited
by any character.  Other examples shown in section 3.2.1 of the SPAGeDi
1.3 manual include \code{"/"}, \code{" "}, \code{", "}, \code{"."}, and
\code{"--"}.
}
  \item{returncatspatcoord}{
Boolean.  Indicates whether a data frame should be returned containing
the category and spatial coordinates columns.
}
  \item{returnploidies}{
Boolean.  Indicates whether a vector should be returned containing the
ploidy of each individual.
}
  \item{missing}{
The symbol to be used to specify missing data in the genotype object
that is returned.
}
}
\details{
SPAGeDi offers a lot of flexibility in how data files are formatted.
\code{read.SPAGeDi} accomodates most of that flexibility.  The primary
exception is that alleles must be delimited in the same way across all
genotypes, as specified by \code{allelesep}.  Comment lines beginning
with \code{//}, as well as blank lines, are ignored by
\code{read.SPAGeDi} just as they are by SPAGeDi.

\code{read.SPAGeDi} is not designed to read dominant data (see section
3.2.2 of the SPAGeDi 1.3 manual).  However, see
\code{dominant.to.codominant} for a way to read this type
of data after some simple manipulation in a spreadsheet program.

The first line of a SPAGeDi file contains information that is used by
\code{read.SPAGeDi}.  The ploidy as specified in the 6th position of the
first line is ignored, and is instead calculated by counting alleles for
each individual (including zeros on the right, but not the left, side of
the genotype).  The number of digits specified in the 5th position of
the first line is only used if \code{allelesep=""}.  All other values
in the first line are important for the function.

If the only alleles found for a particular individual and locus are
zeros, the genotype is interpreted as missing.  Otherwise, zeros on the
left side of a genotype are ignored, and zeros on the right side of a
genotype are used in calculating the ploidy but are not included in the
genotype object that is returned.  If \code{allelesep=""},
\code{read.SPAGeDi} checks that the number of characters in the genotype
can be evenly divided by the number of digits per allele.  If not, zeros
are added to the left of the genotype string before splitting it into
alleles.
}
\value{
Under the default where \code{returncatspatcoord=FALSE} and
\code{returnploidies=FALSE}, a genotype object in the standard polysat
format is returned.  This is a two-dimensional list of integer vectors,
where the first dimension of the list represents samples and the second
dimension of the list represents loci.  Both dimensions are named by the
individual and locus names found in the file.  Each vector contains all
unique alleles, formatted as integers.  

Otherwise, a list of two or three objects is returned:
\item{CatSpatCoord }{A data frame of categories and spatial coordinates,
unchanged from the file.  The format of each column is determined under
the default \code{read.table} settings.  Row names are individual names
from the file.  Column names are the same as in the file.}
\item{Indploidies }{A vector containing the ploidy of each individual,
  as determined by the (maximum) number of alleles per genotype,
  including zeros on the right.  The vector is named by the individual
  names from the file.}
\item{Genotypes }{A genotype object as described above.}
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

\seealso{
\code{\link{write.SPAGeDi}}, \code{\link{dominant.to.codominant}},
\code{\link{read.table}}, \code{\link{read.GeneMapper}},
\code{\link{read.GenoDive}}, \code{\link{read.Structure}},
\code{\link{read.ATetra}}, \code{\link{read.Tetrasat}}
}
\examples{
# create a file to read (usually done with spreadsheet software or a
# text editor):
cat("// here's a comment line at the beginning of the file",
"5\t0\t-2\t2\t2\t4",
"4\t5\t10\t50\t100",
"Ind\tLat\tLong\tloc1\tloc2",
"ind1\t39.5\t-120.8\t00003133\t00004040",
"ind2\t39.5\t-120.8\t3537\t4246",
"ind3\t42.6\t-121.1\t5083332\t40414500",
"ind4\t38.2\t-120.3\t00000000\t41430000",
"ind5\t38.2\t-120.3\t00053137\t00414200",
"END",
sep="\n", file="SpagInputExample.txt")

# display the file
cat(readLines("SpagInputExample.txt"), sep="\n")

# read the file
mydata <- read.SPAGeDi("SpagInputExample.txt", allelesep = "",
returncatspatcoord = TRUE, returnploidies = TRUE)

# view the data
mydata

# view genotypes by locus
mydata$Genotypes[,"loc1"]
mydata$Genotypes[,"loc2"]
}
\keyword{ file }


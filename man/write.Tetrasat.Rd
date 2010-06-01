\name{write.Tetrasat}
\alias{write.Tetrasat}
\title{
Write Genotype Data in Tetrasat Format
}
\description{
  Given a genotype object in the format used by polysat,
  \code{write.Tetrasat} creates a file that can be read by the software
  Tetrasat and Tetra.
}
\usage{
write.Tetrasat(gendata, commentline = "insert data description here",
samples = dimnames(gendata)[[1]], loci = dimnames(gendata)[[2]],
popinfo = rep(1, length(samples)), usatnts = rep(2, length(loci)),
file = "", missing = -9)
}
\arguments{
  \item{gendata}{
Genotype data in the format created and used by other polysat
functions.  A two dimensional list of numerical vectors, where
samples are represented in the first dimension and loci are
represented in the second dimension, and dimensions are named
accordingly.  Each vector contains each unique allele once.
}
  \item{commentline}{
A character string to be written as the first line of the file.
}
  \item{samples}{
A character vector of samples to write to the file.  Should be a subset
of \code{dimnames(gendata)[[1]]}.
}
  \item{loci}{
A character vector of loci to write to the file.  Should be a subset of
\code{dimnames(gendata)[[2]]}.
}
  \item{popinfo}{
An integer vector indicating the population number of each sample.  The
vector should be named using sample names, or should be in the same
order as \code{samples}.
}
  \item{usatnts}{
An integer vector indicating the length of the nucleotide repeats for
each locus.  The vector should be named using locus names, or should be
in the same order as \code{loci}.  The value indicating dinucleotide
repeats is 2, trinucleotide repeats would be represented by 3, etc.
Note that if the alleles in \code{gendata} are already in terms of repeat
number rather than fragment length, the value should be 1.
}
  \item{file}{
A character string indicating the file to which to write.
}
  \item{missing}{
The symbol used to indicate missing data in \code{gendata}.
  }
}
\details{
Tetrasat files are space-delimited text files in which all alleles at a
locus are concatenated into a string eight characters long.  Population
names or numbers are not used in the file, but samples are ordered by
population, with the line \dQuote{Pop} delimiting populations.

\code{write.Tetrasat} divides each allele by the length of the repeat and
rounds down in order to convert alleles to repeat numbers.  If
necessary, it subtracts a multiple of 10 from all alleles at a locus to
make all allele values less than 100, or puts a zero in front of the
number if it only has one digit.  If the individual is fully homozygous
at a locus, the single allele is repeated four times.  If any genotype
has more than four alleles, \code{write.Tetrasat} picks a random sample of
four alleles without replacement, and prints a warning.  Missing data
are represented by blank spaces.

Sample names should be a maximum of 20 characters long in order for the
file to be read correctly by Tetrasat or Tetra.
}
\value{
A file is written but no value is returned.
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
\code{\link{read.Tetrasat}}, \code{\link{write.Structure}},
\code{\link{write.GeneMapper}}, \code{\link{write.ATetra}},
\code{\link{write.GenoDive}}, \code{\link{codominant.to.dominant}},
\code{\link{write.SPAGeDi}}
}
\examples{
# set up sample data (usually done by reading files)
mysamples <- c("ind1", "ind2", "ind3", "ind4")
myloci <- c("loc1", "loc2")
myusatnts <- c(2, 3)
names(myusatnts) <- myloci
mygendata <- array(list(-9), dim=c(4,2), dimnames=list(mysamples, myloci))
mygendata[,"loc1"] <- list(c(202,204), c(204), c(200,206,208,212),
                           c(198,204,208))
mygendata[,"loc2"] <- list(c(78,81,84), c(75,90,93,96,99), c(87), c(-9))
mypopinfo <- c(1,2,1,2)
names(mypopinfo) <- mysamples

# write a Tetrasat file
write.Tetrasat(mygendata, popinfo=mypopinfo, usatnts=myusatnts,
               commentline="sample data", file="tetrasattest.txt")

# view the file
cat(readLines("tetrasattest.txt"),sep="\n")
}
\keyword{file}
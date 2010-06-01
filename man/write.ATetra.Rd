\name{write.ATetra}
\alias{write.ATetra}
\title{
Write Genotypes in ATetra Format
}
\description{
Given a genotype object in the standard polysat format and an index of
which samples belong to which populations, \code{write.ATetra} creates a
text file of genotypes in the ATetra format.
}
\usage{
write.ATetra(gendata, samples = dimnames(gendata)[[1]],
loci = dimnames(gendata)[[2]], popinfo = rep(1, length(samples)),
popnames = "onebigpop", commentline = "insert data info here",
missing = -9, file = "")
}
\arguments{
  \item{gendata}{
A two-dimensional list of integer vectors, in the standard polysat
format.  Samples are represented in the first dimension and loci are
represented in the second dimension, and both dimensions are named
accordingly.  Each vector contains each unique allele once.
}
  \item{samples}{
    A character vector of samples to write to the file.  This is a subset
    of \code{dimnames(gendata)[[1]]}.
}
  \item{loci}{
    A character vector of loci to write to the file.  This is a subset of
    \code{dimnames(gendata)[[2]]}.
}
  \item{popinfo}{
An integer vector where each number represents to which population a sample
belongs.  The vector should have names corresponding to \code{samples},
or if unnamed is assumed to be in the same order as \code{samples}.
}
  \item{popnames}{
A character vector containing the name of each population.  The
populations should be ordered by the numbers used to represent them in
\code{popinfo}.
}
  \item{commentline}{
A character string to be put into the title line of the file.
}
  \item{missing}{
The symbol used to represent missing data in \code{gendata}.
}
  \item{file}{
A character string indicating the path and name to which to write the file.
}
}
\details{
Note that missing data are not allowed in ATetra, although
\code{write.ATetra} will still process missing data.  When it does so, it
leaves all alleles blank in the file for that particular sample and
locus, and also prints a warning indicating which sample and locus had
missing data.

\code{popinfo} can contain information for additional samples that are not
to be used in the file (are not in \code{samples}).  If the number of
poplations in \code{popinfo} is greater than the length of \code{popnames}, a
warning will be printed and only the populations up to the length of
\code{popnames} will be used.  If \code{popnames} is longer than the number of
populations in \code{popinfo}, the same warning will be printed and there
will be a subscripting error.
}
\value{
A file is written but no value is returned.
}
\references{
  \url{http://www.vub.ac.be/APNA/ATetra_Manual-1-1.pdf}

  van Puyvelde, K., van Geert, A. and Triest, L.  (2010)  ATETRA, a new
  software program to analyze tetraploid microsatellite data: comparison
  with TETRA and TETRASAT.  \emph{Molecular Ecology Resources}
  \bold{10}, 331-334.
}
\author{
Lindsay V. Clark
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{read.ATetra}}, \code{\link{write.Tetrasat}},
\code{\link{write.GenoDive}}, \code{\link{write.Structure}},
\code{\link{write.GeneMapper}}, \code{\link{codominant.to.dominant}},
\code{\link{write.SPAGeDi}}
}
\examples{
# set up sample data (usually done by reading files)
mysamples <- c("ind1", "ind2", "ind3", "ind4")
myloci <- c("loc1", "loc2")
mygendata <- array(list(-9), dim=c(4,2), dimnames=list(mysamples, myloci))
mygendata[,"loc1"] <- list(c(202,204), c(204), c(200,206,208,212),
                           c(198,204,208))
mygendata[,"loc2"] <- list(c(78,81,84), c(75,90,93,96,99), c(87), c(-9))
mypopinfo <- c(1,2,1,2)
names(mypopinfo) <- mysamples
mypopnames <- c("this pop", "that pop")

# write an ATetra file
write.ATetra(mygendata, popinfo=mypopinfo, popnames=mypopnames,
             commentline="sample data", file="atetratest.txt")

# view the file
cat(readLines("atetratest.txt"),sep="\n")

}
\keyword{ file }
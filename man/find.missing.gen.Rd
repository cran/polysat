\name{find.missing.gen}
\alias{find.missing.gen}
\title{
Find Missing Genotypes
}
\description{
This function returns a data frame listing the locus and sample names of
all genotypes with missing data.
}
\usage{
find.missing.gen(gendata, missing = -9, samples =
  dimnames(gendata)[[1]], loci = dimnames(gendata)[[2]])
}
\arguments{
  \item{gendata}{
A genotype object in the standard polysat format.  A two dimensional
list of vectors, where samples are represented and named in the first
dimension and loci in the second dimension.  Each vector contains all
unique alleles for a given sample and locus.  If data is missing, the
missing data symbol is the first element of the vector.
}
  \item{missing}{
The symbol used to indicate missing data.
}
  \item{samples}{
A character vector of all samples to be searched.  Must be a subset of
\code{dimnames(gendata)[[1]]}.
}
  \item{loci}{
A character vector of all loci to be searched.  Must be a subset of
\code{dimnames(gendata)[[2]]}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame with no row names.  The first column is named
\dQuote{Locus} and the second column is named \dQuote{Sample}.  Each
row represents one missing genotype, and gives the locus and sample
of that genotype.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# set up the genotype data
samples <- paste("ind", 1:4, sep="")
samples
loci <- paste("loc", 1:3, sep="")
loci
testgen <- array(list(-9), dim = c(4,3), dimnames = list(samples,loci))
testgen[,"loc1"] <- list(c(-9), c(102,104), c(100,106,108,110,114),
                         c(102,104,106,110,112))
testgen[,"loc2"] <- list(c(77,79,83), c(79,85), c(-9), c(83,85,87,91))
testgen[,"loc3"] <- list(c(122,128), c(124,126,128,132), c(120,126),
                         c(124,128,130))

# look up which samples*loci have missing genotypes
find.missing.gen(testgen)

}
\keyword{ manip }


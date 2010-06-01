\name{estimate.freq}
\alias{estimate.freq}
\title{
Estimate Allele Frequencies in Populations
}
\description{
Given genotypes, population identity, and ploidy of each individual,
\code{estimate.freq} produces a data frame showing the estimated frequency
of each allele in each population, as well as the number of genomes
in each population.
}
\usage{
estimate.freq(gendata, missing = -9, samples = dimnames(gendata)[[1]],
loci = dimnames(gendata)[[2]], popinfo = rep(1, length(samples)),
indploidies = rep(4, length(samples)))
}
\arguments{
  \item{gendata}{
A genotype object in the standard polysat format.  A two-dimensional
list of vectors, where samples are represented and named in the first
dimension and loci in the second dimension.  Each vector contains all
unique alleles for a given sample and locus.
}
  \item{missing}{
The symbol used to represent missing data in \code{gendata}.
}
  \item{samples}{
Character vector.  The samples to be used in analysis.  This should be a
subset of \code{dimnames(gendata)[[1]]}.
}
  \item{loci}{
Character vector.  The loci to be used in analysis.  This should be a
subset of \code{dimnames(gendata)[[2]]}.
}
  \item{popinfo}{
Integer or character vector.  The population identity (population number
or name) of each sample.  The names of the vector should correspond to
\code{samples}.  If the vector is unnamed, it is assumed to be in the same
order as \code{samples}.
}
  \item{indploidies}{
Integer vector.  The ploidy of each sample.  Should be named similarly
to \code{popinfo}, or if unnamed is assumed to be in the same order as
\code{samples}.
}
}
\details{
This function estimates allele frequencies rather than calculating them
exactly from the sample, because if there are any partially heterozygous
genotypes present then allele copy number cannot be known exactly.  For
each sample*locus, a conversion factor is generated that is the ploidy
of the sample as specified in \code{indploidies} divided by the number of
alleles that the sample has at that locus.  Each allele is then
considered to be present in as many copies as the the conversion factor
(note that this is not necessarily an integer).  The number of copies of
an allele is totaled for the population and is divided by the total
number of genomes in the population (minus missing data at the locus)
in order to calculate allele frequency.

A major assumption of this calculation method is that each allele in a
partially heterozygous genotype has an equal chance of being present in
more than one copy.  This is almost never true, because common alleles
in a population are more likely to be partially homozygous in an
individual.  The result is that the frequency of common alleles is
underestimated and the frequency of rare alleles is overestimated.
}
\value{
Data frame, where each population is in one row.  The first column is
called \code{Genomes} and contains the number of genomes in each
population.  Each remaining column contains frequencies for one allele.
Columns are named by locus and allele, separated by a period.
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
\code{\link{calcFst}}
}
\examples{
# create a data set (typically done by reading files)
mygenotypes <- array(list(-9), dim = c(6,2), dimnames =
                     list(paste("ind",1:6, sep=""),c("loc1","loc2")))
mygenotypes[,"loc1"] <- list(c(206),c(208,210),c(204,206,210),
    c(196,198,202,208),c(196,200),c(198,200,202,204))
mygenotypes[,"loc2"] <- list(c(130,134),c(138,140),c(130,136,140),
    c(138),c(136,140),c(130,132,136))

mypopinfo <- c(1,1,1,2,2,2)
names(mypopinfo) <- dimnames(mygenotypes)[[1]]

myploidies <- c(2,2,4,4,2,4)
names(myploidies) <- dimnames(mygenotypes)[[1]]

# calculate allele frequencies
myfreq <- estimate.freq(mygenotypes, popinfo=mypopinfo,
indploidies=myploidies)

# look at the results
myfreq
}
\keyword{ arith }


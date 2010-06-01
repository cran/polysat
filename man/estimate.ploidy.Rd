\name{estimate.ploidy}
\alias{estimate.ploidy}
\title{
Maximum and Mean Allele Count for Estimation of Ploidy
}
\description{
Given genotypes in the form of a two-dimensional list of vectors,
estimate.ploidy produces a two-dimensional array containing the maximum
number of alleles and the mean number of alleles for each sample, across
all loci.
}
\usage{
estimate.ploidy(gendata, samples = dimnames(gendata)[[1]],
loci = dimnames(gendata)[[2]])
}
\arguments{
  \item{gendata}{
Genotypes in the standard polysat format.  A two-dimensional list of
vectors, where samples are represented and named in the first dimension
and loci in the second dimension.  Each vector contains all unique
alleles for a given sample and locus.
}
  \item{samples}{
    An optional character vector of samples to evaluate, which is a subset
    of \code{dimnames(gendata)[[1]]}.
}
\item{loci}{
  An optional character vector of loci to use in the calculation, which
  is a subset of \code{dimnames(gendata)[[2]]}.
}
}
\details{
  To assist the user in determining the ploidy of each sample,
  \code{estimate.ploidy} looks at the genotype of the sample across
  all loci and returns the maximum number of alleles per locus.  The
  mean number of alleles is also returned to assist with checking for
  errors (for example, if an octoploid genotype was accidentally scored
  at one locus for a diploid sample).  Both of these are calculated
  using the \code{length} function on the genotype vectors.

  The user may want to extract the vector containing the maximum number
  of alleles (for example, \code{myploidies<-ploidyinfo[,1]}) and then
  manually edit the values based on other knowledge of the organism.
  This vector can then be used as the \code{indploidies} argument for
  \code{write.Structure} or \code{estimate.freq}.
}
\value{
  An array with the second dimension of length 2 and the first dimension
  as long as \code{samples}.  The rows are labeled by sample name and
  the columns are labeled \code{max.alleles} and \code{mean.alleles}.
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
  \code{\link{estimate.freq}}, \code{\link{write.Structure}}
}
\examples{
# Create a data set to analyze
mygendata <-
  array(list(c(124,128,138),c(122,130,140,142),c(122,132,136),
             c(122,134,140),
             c(203,212,218),c(197,206,221),c(215),c(200,218),
             c(140,144,148,150),c(-9),c(146,150),c(152,154,158))
        , dim=c(4,3), dimnames=list(c("ind1","ind2","ind3","ind4"),
        c("locus1","locus2","locus3")))

# Run the function
estimate.ploidy(mygendata)
}

\keyword{arith}
\keyword{array}
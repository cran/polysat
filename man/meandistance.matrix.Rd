\name{meandistance.matrix}
\alias{meandistance.matrix}
\title{Mean Pairwise Distance Matrix}
\description{
  Given a two-dimensional list of genotypes, indexed by sample and
  locus, \code{meandistance.matrix} produces a symmetrical matrix of
  pairwise distances between samples, averaged across all loci.  An
  array of all distances prior to averaging may also be produced.
}
\usage{
meandistance.matrix(gendata, samples=dimnames(gendata)[[1]],
loci=dimnames(gendata)[[2]], all.distances=FALSE, usatnts=NULL, ...)
}
\arguments{
  \item{gendata}{A two-dimensional list of genotypes in the standard
    polysat format.  The first dimension is an index of samples, the
    second dimension is an index of loci, and the elements are numerical
    vectors containing the alleles as elements.}
  \item{samples}{A character vector of samples to be analyzed.  These
    should be all or a subset of the sample names used in \code{gendata}.}
  \item{loci}{A character vector of loci to be analyzed.  These should
    be all or a subset of the loci names used in \code{gendata}.}
  \item{all.distances}{If FALSE, only the mean distance matrix will be
    returned.  If TRUE, a list will be returned containing an array of
    all distances by locus and sample as well as the mean distance matrix.}
  \item{usatnts}{A numerical vector that contains the length of
    nucleotide repeats for each locus.  For example, 3 would be used to
    indicate a locus with trinucleotide repeats.  1 should be used if
    alleles are written in terms of repeat number, not fragment length
    in nucleotides. \code{names(usatnts)} should be the same as those used in
    \code{dimnames(gendata)[[2]]} (the names of the loci).  This argument can be
    omitted if repeat length is irrelevant to the distance metric.}
  \item{...}{If \code{distmetric} or \code{progress} are given here they will be
    passed to \code{distance.matrix.1locus}.  Any other arguments will be
    passed to \code{distmetric}.}
}
\value{
  A symmetrical matrix containing pairwise distances between all
  samples, averaged across all loci.  Row and column names of the matrix
  will be the sample names provided in the \sQuote{samples} argument.  If
  \code{all.distances=TRUE}, a list will be produced containing the above
  matrix as well as a three-dimensional array containing all distances
  by locus and sample.  The array is the first item in the list, and the
  mean matrix is the second.
}
\details{
  \code{meandistance.matrix} uses \code{distance.matrix.1locus} once
  for each locus to be analyzed, then averages values across these
  matrices.  Any arguments that need to be passed to
  \code{distance.matrix.1locus} may be given to
  \code{meandistance.matrix}.  If the loci are of different repeat
  types and the type of repeat is important for the distance metric
  being used (e.g. \code{Bruvo.distance}), the \code{usatnts} argument
  can be used to pass a different \code{usatnt} argument to
  \code{distmetric} depending on the locus.

  Because the user may want to omit samples or loci, the \code{samples}
  and \code{loci} arguments are given for convenient indexing of the
  data to be analyzed.  If \code{gendata} contains only the data that
  the user wants to analyze, the user can simply omit these arguments.

  Missing data must be represented by the missing data symbol, rather
  than \code{NA}.
}
\references{
}
\seealso{
  \code{\link{distance.matrix.1locus}}, \code{\link{Bruvo.distance}},
  \code{\link{Lynch.distance}}
}
\examples{
# create a list of genotype data
mygendata <-
  array(list(c(124,128,138),c(122,130,140,142),c(122,132,136),c(122,134,140),
             c(203,212,218),c(197,206,221),c(215),c(200,218),
             c(140,144,148,150),c(-9),c(146,150),c(152,154,158),
             c(233,236,280),c(-9),c(-9),c(-9))
        ,dim=c(4,4),dimnames=list(c("ind1","ind2","ind3","ind4"),
                                  c("locus1","locus2","locus3","locus4")))

# make index vectors of data to use
myloci <- c("locus1","locus2","locus3")
mysamples <- c("ind1","ind2","ind4")

# locus1 and locus3 have dinucleotide repeats, and locus2 has
# trinucleotide repeats
myusatnts <- c(2,3,2)
names(myusatnts) <- myloci

meandistance.matrix(mygendata, mysamples, myloci, all.distances=TRUE,
                     usatnts=myusatnts)
}
\author{Lindsay V. Clark}
\keyword{array}
\keyword{arith}
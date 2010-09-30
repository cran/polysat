\name{meandistance.matrix}
\alias{meandistance.matrix}
\title{Mean Pairwise Distance Matrix}
\description{
  Given a \code{genambig} object, \code{meandistance.matrix}
  produces a symmetrical matrix of
  pairwise distances between samples, averaged across all loci.  An
  array of all distances prior to averaging may also be produced.
}
\usage{
meandistance.matrix(object, samples = Samples(object),
                    loci = Loci(object), all.distances=FALSE,
                    distmetric = Bruvo.distance, progress = TRUE,
                    ...)
}
\arguments{
  \item{object}{A \code{genambig} object containing the genotypes to be
    analyzed.  If \code{distmetric = Bruvo.distance}, the \code{Usatnts}
    slot should be filled in.}
  \item{samples}{A character vector of samples to be analyzed.  These
    should be all or a subset of the sample names used in \code{object}.}
  \item{loci}{A character vector of loci to be analyzed.  These should
    be all or a subset of the loci names used in \code{object}.}
  \item{all.distances}{If \code{FALSE}, only the mean distance matrix will be
    returned.  If \code{TRUE}, a list will be returned containing an array of
    all distances by locus and sample as well as the mean distance matrix.}
  \item{distmetric}{The function to be used to calculate distances
    between genotypes.  \code{Bruvo.distance}, \code{Lynch.distance}, or
  a distance function written by the user.}
  \item{progress}{If \code{TRUE}, loci and samples will be printed to
    the console as distances are calculated, so that the user can
    monitor the progress of the computation.}
  \item{...}{Additional arguments (such as \code{maxl}) to
    pass to \code{distmetric}.}
}
\value{
  A symmetrical matrix containing pairwise distances between all
  samples, averaged across all loci.  Row and column names of the matrix
  will be the sample names provided in the \code{samples} argument.  If
  \code{all.distances=TRUE}, a list will be produced containing the above
  matrix as well as a three-dimensional array containing all distances
  by locus and sample.  The array is the first item in the list, and the
  mean matrix is the second.
}
\details{
  Each distance for the three-dimensional array is calculated only once,
  to save computation time.  Since the array (and resulting mean matrix)
  is symmetrical, the distance is written to two positions in the array
  at once.
}
\references{
}
\seealso{
  \code{\link{Bruvo.distance}}, \code{\link{Lynch.distance}},
  \code{\link{meandist.from.array}}
}
\examples{
# create a list of genotype data
mygendata <- new("genambig", samples = c("ind1","ind2","ind3","ind4"),
                 loci = c("locus1","locus2","locus3","locus4"))
Genotypes(mygendata) <-
  array(list(c(124,128,138),c(122,130,140,142),c(122,132,136),c(122,134,140),
             c(203,212,218),c(197,206,221),c(215),c(200,218),
             c(140,144,148,150),c(-9),c(146,150),c(152,154,158),
             c(233,236,280),c(-9),c(-9),c(-9)))
Usatnts(mygendata) <- c(2,3,2,1)

# make index vectors of data to use
myloci <- c("locus1","locus2","locus3")
mysamples <- c("ind1","ind2","ind4")

# calculate array and matrix
mymat <- meandistance.matrix(mygendata, mysamples, myloci,
                             all.distances=TRUE)
# view the results
mymat[[1]]["locus1",,]
mymat[[1]]["locus2",,]
mymat[[1]]["locus3",,]
mymat[[2]]

}
\author{Lindsay V. Clark}
\keyword{array}
\keyword{arith}
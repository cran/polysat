\name{distance.matrix.1locus}
\alias{distance.matrix.1locus}
\title{Pairwise Genetic Distances at One Locus}
\description{
  Given all genotypes for one locus, create a pairwise genetic distance
  matrix.
}
\usage{
distance.matrix.1locus(gendata, distmetric = Bruvo.distance,
progress = TRUE, ...)
}
\arguments{
  \item{gendata}{
  A sublist from a genotype object in the standard polysat format,
  representing all samples for one locus.  A list of vectors, where
  each vector contains all unique alleles for a given sample at this
  locus.  Sample names are used to name the list.}
  \item{distmetric}{This is the function that will be used to calculate
    each pairwise distance.  This should be a function that, given two
    vectors of alleles, returns a numerical distance.}
  \item{progress}{If \code{TRUE}, \code{distance.matrix.1locus} will
    print the names of sample pairs as it finishes each calculation
    with \code{distmetric}.  For large datasets, this is intended so
    that the user can monitor the progress of the calculations.}
  \item{...}{These arguments will be passed to \code{distmetric}.  For
    example, with \code{Bruvo.distance}, \code{maxl}, \code{usatnt}, or
    \code{missing} may be used.}
}
\value{
  A symmetrical matrix of distances, with the names of samples used as
  row and column names.
}
\details{
  Given a list of genotypes at one locus, \code{distance.matrix.1locus}
  produces a symmetrical matrix of pairwise distances between
  genotypes.  If using a polysat genotype object such as that produced
  by \code{read.GeneMapper}, the \code{gendata} argument should be one
  of the sublists, for example \code{mygenotypedata[,"locus1"]}.
  The measure of distance can be any that is provided with polysat, or
  any function written by the user, so long as it takes genotypes as
  vectors  of alleles (or any other type of object that is given as
  elements of the list \code{gendata}) and returns a numerical distance.
  Any arguments that need to be passed to the \code{distmetric} function
  can be given to \code{distance.matrix.1locus}.  To save processing
  time, each pairwise distance is only calculated once and then written
  to both locations in the matrix simultaneously.  The user also has the
  option to have each pair of sample names printed after the distance is
  calculated, so that progress can be monitored if evaluation is
  expected to take a long time.
}
\references{
}
\seealso{
  \code{\link{Bruvo.distance}}, \code{\link{Lynch.distance}}, \code{\link{meandistance.matrix}}
}
\examples{
mygenotypes <- list(IND1=c(124,127,133), IND2=c(130,139,145,151),
    IND3=c(118,127,133,154))
distance.matrix.1locus(mygenotypes,usatnt=3)
}
\author{Lindsay V. Clark}
\keyword{array}
\name{write.Structure}
\alias{write.Structure}
\title{Write Genotypes in Structure 2.3 Format}
\description{
  Given genotypes in the form of a two-dimensional list of vectors,
  \code{write.Structure} produces a text file of the genotypes in a format
  readable by Structure 2.2 and higher.  The user specifies the overall
  ploidy of the file as well as the ploidy of each sample.
}
\usage{
write.Structure(gendata, ploidy, file="",
samples=dimnames(gendata)[[1]], loci=dimnames(gendata)[[2]],
indploidies=rep(ploidy,times=length(samples)),
extracols=NULL, missing=-9)
}
\arguments{
  \item{gendata}{Genotypes stored as a two-dimensional list of vectors,
    in the standard polysat format.}
  \item{ploidy}{PLOIDY for Structure, i.e. how many rows per individual
    to write.}
  \item{file}{A character string of the file path where the file
    should be written to.}
  \item{samples}{An optional character vector listing the names of samples to be
    written to the file.}
  \item{loci}{An optional character vector listing the names of the loci to be
    written to the file.}
  \item{indploidies}{An integer vector containing the ploidy of each
    sample. \code{names(indploidies)} should correspond to
    \code{samples}, or if the vector is unnamed it is assumed to be in
    the same order as \code{samples}.}
  \item{extracols}{An array, with the first dimension names
    corresponding to \code{samples}, of PopData, PopFlag, LocData, Phenotype,
    or other values to be included in the extra columns in the file.}
  \item{missing}{The number used to indicate missing data.}
}
\value{
  No value is returned, but instead a file is written at the path specified.
}
\details{
  Structure 2.2 and higher can process polyploid microsatellite data,
  although 2.3.3 or higher is recommended for its improvements on
  polyploid handling.  The input format of Structure requires that
  each locus take up one column and that each individual take up as
  many rows as the parameter PLOIDY.  Because of the multiple rows per
  sample, each sample name must be duplicated, as well as any
  population, location, or phenotype data.  Partially heterozygous
  genotypes also must have one arbitrary allele duplicated up to the
  ploidy of the sample, and samples that have a lower ploidy than that
  used in the file (for mixed polyploid data sets) must have a missing
  data symbol inserted up to fill in the extra rows.  Additionally, if
  some samples have more alleles than PLOIDY (if you are using a lower
  PLOIDY to save processing time, or if there are extra alleles from
  scoring errors), some alleles must be randomly removed from the data.
  \code{write.Structure} performs this duplication, insertion, and random
  deletion of data.

  The argument \code{samples} contains all of the sample names to be
  written to the file, and is used to index \code{gendata}, \code{indploidies},
  and \code{extracols}.  These sample names will also be used as row
  names in the Structure file.  Each sample name should only be in the
  vector sample once, because \code{write.Structure} will duplicate
  the sample names a number of times as dictated by \code{ploidy}.
  Likewise, \code{indploidies} and \code{extracols} only need to contain
  data for each sample once.  If \code{samples} isn't specified by the user it
  will be extracted from \code{gendata}.

  In writing genotypes to the file, \code{write.Structure} compares the number
  of alleles in the genotype, the ploidy of the sample as stored in
  \code{indploidies}, and the ploidy of the file as stored in
  \code{ploidy}, and does one of six things (for a given sample x and
  locus loc):

  1) If indploidies[x] is greater than or equal to ploidy, and
  length(gendata[[x,loc]]) is equal to ploidy, the genotype data is
  used as is.
  
  2) If indploidies[x] is greater than or equal to ploidy, and
  length(gendata[[x,loc]]) is less than ploidy, the first allele is
  duplicated as many times as necessary for there to be as many alleles
  as ploidy.
  
  3) If indploidies[x] is greater than or equal to ploidy, and
  length(gendata[[x,loc]]) is greater than ploidy, a random sample of
  the alleles, without replacement, is used as the genotype.
  
  4) If indploidies[x] is less than ploidy, and
  length(gendata[[x,loc]]) is equal to indploidies[x], the genotype
  data is used as is and missing data symbols are inserted in the extra
  rows.
  
  5) If indploidies[x] is less than ploidy, and
  length(gendata[[x,loc]]) is less than indploidies[x], the first
  allele is duplicated as many times as necessary for there to be as
  many alleles as indploidies[x], and missing data symbols are inserted
  in the extra rows.
  
  6) If indploidies[x] is less than ploidy, and
  length(gendata[[x,loc]]) is greater than indploidies[x], a random
  sample, without replacement, of indploidies[x] alleles is used, and
  missing data symbols are inserted in the extra rows. (Alleles are
  removed even though there is room for them in the file.)

  Two of the header rows that are optional for Structure are written by
  \code{write.Structure}.  These are \sQuote{Marker Names}, containing the
  names of loci supplied in \code{gendata}, and \sQuote{Recessive
    Alleles}, which contains the missing data symbol once for each locus.
  This indicates to the program that all alleles are codominant with
  copy number ambiguity.

  The output file requires a few small modifications, done in a text
  editor or spreadsheet software, in order to be read
  by Structure.  In the upper left corner the words \dQuote{rowlabel} and
  \dQuote{missing} should be deleted.  Likewise the first and second rows for
  any non-locus columns should be deleted if the extracols argument was
  used.  These should include the second dimension names used in
  \code{extracols}, and zeros, respectively.
}
\references{
  \url{http://pritch.bsd.uchicago.edu/structure_software/release_versions/v2.3.3/structure_doc.pdf}

  Hubisz, M. J., Falush, D., Stephens, M. and Pritchard, J. K. (2009)
  Inferring weak population structure with the assistance of sample
  group information.  \emph{Molecular Ecology Resources} \bold{9}, 1322-1332.

  Falush, D., Stephens, M. and Pritchard, J. K.  (2007)
  Inferences of population structure using multilocus genotype data:
  dominant markers and null alleles.  \emph{Molecular Ecology Notes}
  \bold{7}, 574-578.
}
\note{If \code{extracols} is a character array, make sure none of the
  elements contain white space.}
\seealso{
  \code{\link{read.Structure}}, \code{\link{write.GeneMapper}},
  \code{\link{write.GenoDive}}, \code{\link{write.ATetra}},
  \code{\link{write.Tetrasat}}, \code{\link{codominant.to.dominant}},
  \code{\link{write.SPAGeDi}}
}
\examples{
# input genotype data (this is usually done by reading a file)
mygendata <- array(list(c(100,102,106,108,114,118),c(102,110),
                      c(98,100,104,108,110,112,116),c(102,106,112,118),
                      c(104,108,110),c(-9),
                      c(204),c(206,208,210,212,220,224,226),
                      c(202,206,208,212,214,218),c(200,204,206,208,212),
                      c(-9),c(202,206)),
                 dim=c(6,2), dimnames=list(c("ind1","ind2","ind3",
                                             "ind4","ind5","ind6"),
                                           c("locus1","locus2")))

# create a vector of sample names to be used.  Note that this excludes
#  ind6.
# Also note that this could be obtained as names(mygendata[[1]]).
mysamples <- c("ind1","ind2","ind3","ind4","ind5")

# create a vector of the ploidy of each sample.
# Note that some of the above genotypes have more or fewer alleles than
# the ploidy of the sample.
myploidies <- c(6,6,6,4,4)
names(myploidies) <- mysamples

# Create an array containing data for additional columns to be written
# to the file.  You might also prefer to just read this and the ploidies
# in from a file.
myexcols <- array(data=c(1,2,1,2,1,1,1,0,0,0),dim=c(5,2),
                  dimnames=list(mysamples, c("PopData","PopFlag")))

# Write the Structure file, with six rows per individual.
# Since outfile="", the data will be written to the console instead of a file.
write.Structure(mygendata, 6, "", samples=mysamples, indploidies=myploidies,
                extracols=myexcols)
}
\author{Lindsay V. Clark}
\keyword{file}
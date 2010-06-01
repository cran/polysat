\name{read.GeneMapper}
\alias{read.GeneMapper}
\title{Read GeneMapper Genotypes Tables}
\description{
  Given a vector of filepaths to tab-delimited text files containing
  genotype data in the ABI GeneMapper Genotypes Table format,
  \code{read.GeneMapper} produces a two-dimensional list of genotypes that can
  be read by other functions in the polysat package.
}
\usage{
  read.GeneMapper(infiles, missing = -9)
}
\arguments{
  \item{infiles}{A character vector of paths to the files to be read.}
  \item{missing}{A numerical value used to indicate missing data for a
    given sample and locus.}
}
\value{
  The object produced is a two-dimensional list of vectors representing
  the genotypes.  Samples are represented by the first dimension and
  loci by the second dimension.  The names of the samples and loci are
  taken from the Sample Names and Markers columns, respectively, of the
  GeneMapper files.  The vectors at each position of the list are
  numerical and are only as long as needed to contain each allele (for
  that sample and locus) as one element.
}
\details{
  \code{read.GeneMapper} can read the genotypes tables that are exported
  by the Applied Biosystems GeneMapper software.  The only alterations
  to the files that the user may have to make are 1) delete
  any rows with missing data or fill in a numerical missing
  data symbol of your choice (such as -9) in the first allele slot for
  that row, 2) make sure that all allele names are numeric
  representations of fragment length (no question marks or dashes), and
  3) put sample names into the Sample Name column, if the names that you
  wish to use in analysis are not already there.  Each file should have
  the standard header row produced by the software.  If any sample has
  more than one genotype listed for a given locus, only the last
  genotype listed will be used.

  The file format is simple enough that the user can easily create files
  manually if GeneMapper is not the software used in allele calling.
  The files are tab-delimited text files.  There should be a header row
  with column names.  The column labeled \dQuote{Sample Name} should contain
  the names of the samples, and the column labeled \dQuote{Marker} should
  contain the names of the loci.  You can have as many or as few columns as
  needed to contain the alleles, and each of these columns should be
  labeled \dQuote{Allele X} where X is a number unique to each column.  Row
  labels and any other columns are ignored.  For any given sample, each
  allele is listed only once and is given as an integer that is the
  length of the fragment in nucleotides.  Alleles are separated by
  tabs.  If you have more allele columns than alleles for any given
  sample, leave the extra cells blank so that \code{read.table} will
  read them as \code{NA}.  Example data files in this format are
  included in the package.

  \code{read.GeneMapper} will read all of your data at once.  It takes
  as its first argument a character vector containing paths to all of
  the files to be read.  How the data are distributed over these files
  does not matter.  The function finds all unique sample names and all
  unique markers across all the files, and automatically puts a missing
  data symbol into the list if a particular sample and locus combination
  is not found.  Rows in which all allele cells are blank should NOT be
  included in the input files; either delete these rows or put the
  missing data symbol into the first allele cell.

  Sample and locus names must be consistent within and across the
  files.  The list that is produced is indexed by these names.  For
  example, if the object produced was called \code{mygenotypes},
  \code{mygenotypes[["AB1","ABC5"]]} would be a vector containing alleles for
  sample AB1 at locus ABC5.  \code{mygenotypes[,"ABC5"]} would be a list of all
  genotypes at locus ABC5, and \code{mygenotypes["AB1",]} would be a list of
  all genotypes for sample AB1.
}
\references{
  \url{http://www.appliedbiosystems.com/genemapper}
}
\note{
  A \sQuote{subscript out of bounds} error may mean that a sample name
  or marker was left blank in one of the input files.
  }
\seealso{
  \code{\link{read.Tetrasat}}, \code{\link{read.ATetra}},
  \code{\link{read.Structure}}, \code{\link{read.GenoDive}},
  \code{\link{write.GeneMapper}}, \code{\link{dominant.to.codominant}},
  \code{\link{read.SPAGeDi}}
}
\examples{
  \dontrun{
    myinfiles<-c("data\\sample CBA15.txt","data\\sample
    CBA23.txt","data\\sample CBA28.txt")
    mygenotypes<-read.GeneMapper(myinfiles)

    #Look at the object produced.  Alleles are not listed but you can
    #see that the array was filled.
    mygenotypes

    #Look at the genotype of individual FCR5.
    mygenotypes["FCR5",]

    #Correct one of the genotypes
    mygenotypes[["FCR5","RhCBA15"]]<-c(208)
  }

# an example with defined data:
# create a table of data
gentable <- data.frame(Sample.Name=rep(c("ind1","ind2","ind3"),2),
                       Marker=rep(c("loc1","loc2"), each=3),
                       Allele.1=c(202,200,204,133,133,130),
                       Allele.2=c(206,202,208,136,142,136),
                       Allele.3=c(NA,208,212,145,148,NA),
                       Allele.4=c(NA,216,NA,151,157,NA)
                       )
# create a file (inspect this file in a text editor or spreadsheet
# software to see the required format)
write.table(gentable, file="readGMtest.txt", quote=FALSE, sep="\t",
            na="", row.names=FALSE, col.names=TRUE)

# read the file
mygenotypes <- read.GeneMapper("readGMtest.txt")

# inspect the results
mygenotypes[,"loc1"]
mygenotypes[,"loc2"]
}
\author{Lindsay V. Clark}
\keyword{file}
\name{read.ATetra}
\alias{read.ATetra}
\title{
Read File in ATetra Format
}
\description{
Given a file formatted for the software ATetra, \code{read.ATetra}
produces a vector of the population identities of each sample, as well
as a two-dimensional list of vectors representing the genotypes.
}
\usage{
read.ATetra(infile)
}
\arguments{
  \item{infile}{
Character string.  A file path to the file to be read.
}
}
\details{
\code{read.ATetra} reads text files in the exact format specified by
the ATetra documentation.  Note that this format only allows tetraploid
(or lower) data and that there can be no missing data.  The genotype
object produced is in the same format that is produced or read by other
functions in the polysat package.
}
\value{
\item{PopData}{A numerical vector containing the population number for
  each sample (individual).  The sample names are used as the names of
  this vector.}
\item{Genotypes}{A two dimensional list of integer vectors in the
  standard polysat genotype format.  The first dimension
represents samples and the second represents loci, and both are named
accordingly.  Each vector contains as elements the alleles for that
sample at that locus.}
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
\code{\link{read.Tetrasat}}, \code{\link{read.GeneMapper}},
\code{\link{read.Structure}}, \code{\link{read.GenoDive}},
\code{\link{dominant.to.codominant}}, \code{\link{write.ATetra}},
\code{\link{read.SPAGeDi}}
}
\examples{
# create a file to be read
# (this would normally be done in a text editor or with ATetra's Excel template)
cat("TIT,Sample Rubus Data for ATetra", "LOC,1,CBA15",
"POP,1,1,Commonwealth", "IND,1,1,1,CMW1,197,208,211,213",
"IND,1,1,2,CMW2,197,207,211,212", "IND,1,1,3,CMW3,197,208,212,219",
"IND,1,1,4,CMW4,197,208,212,219", "IND,1,1,5,CMW5,197,208,211,212",
"POP,1,2,Fall Creek Lake", "IND,1,2,6,FCR4,197,207,211,212",
"IND,1,2,7,FCR7,197,208,212,218", "IND,1,2,8,FCR14,197,207,212,218",
"IND,1,2,9,FCR15,197,208,211,212", "IND,1,2,10,FCR16,197,208,211,212",
"IND,1,2,11,FCR17,197,207,212,218","LOC,2,CBA23","POP,2,1,Commonwealth",
"IND,2,1,1,CMW1,98,100,106,125","IND,2,1,2,CMW2,98,125,,",
"IND,2,1,3,CMW3,98,126,,","IND,2,1,4,CMW4,98,106,119,127",
"IND,2,1,5,CMW5,98,106,125,","POP,2,2,Fall Creek Lake",
"IND,2,2,6,FCR4,98,125,,","IND,2,2,7,FCR7,98,106,126,",
"IND,2,2,8,FCR14,98,127,,","IND,2,2,9,FCR15,98,108,117,",
"IND,2,2,10,FCR16,98,125,,","IND,2,2,11,FCR17,98,126,,","END",
file = "atetraexample.txt", sep = "\n")

# Read the file and examine the data
exampledata <- read.ATetra("atetraexample.txt")
exampledata$PopData
exampledata$Genotypes
exampledata$Genotypes[["CMW1","CBA23"]]
}
\keyword{file}
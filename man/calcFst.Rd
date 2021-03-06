\name{calcPopDiff}
\alias{calcFst}
\alias{calcPopDiff}
\title{
Estimate Population Differentiation Statistics
}
\description{
Given a data frame of allele frequencies and population sizes,
\code{calcPopDiff} calculates a matrix of pairwise \eqn{F_{ST}}{FST}, 
\eqn{G_{ST}}{GST}, Jost's \eqn{D}{D}, or \eqn{R_{ST}}{RST} values, or a single global
value for any of these four statistics.
\code{calcFst} is a wrapper for \code{calcPopDiff} to allow backwards
compatibility with previous versions of \pkg{polysat}.
}
\usage{
calcPopDiff(freqs, metric, pops = row.names(freqs), 
            loci = unique(gsub("\\\\..*$", "", names(freqs))), global = FALSE,
            bootstrap = FALSE, n.bootstraps = 1000, object = NULL)

calcFst(freqs, pops = row.names(freqs),
        loci = unique(gsub("\\\\..*$", "", names(freqs))), ...)
}
\arguments{
  \item{freqs}{
A data frame of allele frequencies and population sizes such as that
produced by \code{simpleFreq} or \code{deSilvaFreq}.  Each population 
is in one row, and a
column called \code{Genomes} (or multiple columns containing the locus
names and \dQuote{Genomes} seperated by a period) contains the relative
size of each
population.  All other columns contain allele frequencies.  The names of
these columns are the locus name and allele name, separated by a period.
}
\item{metric}{
The population differentiation statistic to estimate.  Can be \dQuote{Fst},
\dQuote{Gst}, \dQuote{Jost's D}, or \dQuote{Rst}.
}
  \item{pops}{
A character vector.  Populations to analyze, which should be
a subset of \code{row.names(freqs)}.
}
  \item{loci}{
A character vector indicating which loci to analyze.  These should be a
subset of the locus names as used in the column names of \code{freqs}.
}
\item{global}{
Boolean.  If \code{TRUE}, a single global statistic will be estimated across all
populations.  If \code{FALSE}, pairwise statistics will be estimated between 
populations.
}
\item{bootstrap}{
Boolean.  If \code{TRUE}, a set of replicates bootstrapped across loci will be returned.  If
\code{FALSE}, a single value will be returned for each pair of populations (if 
\code{global = FALSE}) or for the whole set (if \code{global = TRUE}).
}
\item{n.bootstraps}{
Integer.  The number of bootstrap replicates to perform.  Ignored if 
\code{bootstrap = FALSE}.
}
\item{object}{
A \code{"genambig"} or \code{"genbinary"} object with the \code{Usatnts} slot filled in.
Required for \code{metric = "Rst"} only.
}
\item{...}{
Additional arguments to be passed to \code{calcPopDiff} (\emph{i.e.} \code{global},
\code{bootstrap}, and/or \code{n.bootstraps}).
}
}
\details{
For \code{metric = "Fst"} or \code{calcFst}:

\eqn{H_S}{HS} and \eqn{H_T}{HT} are estimate directly from allele frequencies
for each locus for each pair of populations, then averaged across loci.  
Wright's \eqn{F_{ST}}{FST} is then calculated for each pair of populations as 
\eqn{\frac{H_T - H_S}{H_T}}{(HT-HS)/HT}.

\eqn{H}{H} values (expected heterozygosities for populations and combined
populations) are calculated as one minus the sum of all squared allele
frequencies at a locus.  To calculte \eqn{H_T}{HT}, allele frequencies between two
populations are averaged before the calculation.  To calculate \eqn{H_S}{HS}, \eqn{H}{H}
values are averaged after the calculation.  In both cases, the averages
are weighted by the relative sizes of the two populations (as indicated
by \code{freqs$Genomes}).

For \code{metric = "Gst"}:

This metric is similar to \eqn{F_{ST}}{FST}, but mean allele frequencies and 
mean expected heterozygosities are not weighted by population size.  Additionally,
unbiased estimators of \eqn{H_S}{HS} and \eqn{H_T}{HT} are used according to Nei
and Chesser (1983; equations 15 and 16, reproduced also in Jost (2008)).
Instead of using twice the harmonic mean of the number of individuals in the two
subpopulations as described by Nei and Chesser (1983), the harmonic mean of the
number of allele copies in the two subpopulations (taken from \code{freq$Genomes})
is used, in order to allow for polyploidy.
\eqn{G_{ST}}{GST} is estimated for each locus and then averaged across loci.

For \code{metric = "Jost's D"}:

The unbiased estimators of \eqn{H_S}{HS} and \eqn{H_T}{HT} and calculated as 
with \eqn{G_{ST}}{GST}, without weighting by population size.  They are then used 
to estimate \eqn{D}{D} at each locus according to Jost (2008; equation 12):

\deqn{2 * \frac{H_T - H_S}{1 - H_S}}{2 * (HT - HS)/(1 - HS)}

Values of \eqn{D}{D} are then averaged across loci.

For \code{metric = "Rst"}:

\eqn{R_{ST}}{RST} is estimated on a per-locus basis according to Slatkin (1995),
but with populations weighted equally regardless of size. Values are then averaged
across loci.

For each locus:

\deqn{S_w = \frac{1}{d} * \sum_j^d{\frac{\sum_i{\sum_{i' < i}{p_{ij} * p_{i'j} * n_j^2 * (a_i - a_{i'})^2}}}{n_j * (n_j - 1)}}}{S_w = 1/d * \sum_j^d (\sum_i \sum_{i' < i} p_ij * p_i'j * n_j^2 * (a_i - a_i')^2)/(n_j * (n_j - 1)))}

\deqn{\bar{S} = \frac{\sum_i{\sum_{i' < i}{\bar{p}_i * \bar{p}_{i'} * n^2 * (a_i - a_{i'})^2}}}{n * (n-1)}}{S = (\sum_i \sum_{i' < i} p_i * p_i' * n^2 * (a_i - a_i')^2)/(n * (n - 1))}

\deqn{R_{ST} = \frac{\bar{S} - S_w}{\bar{S}}}{RST = (S - S_w)/S}

where \eqn{d}{d} is the number of populations, \eqn{j}{j} is an individual population, \eqn{i}{i} and 
\eqn{i'}{i'} are individual alleles, \eqn{p_{ij}}{p_ij} is the frequency of an allele in a population,
\eqn{n_j}{n_j} is the number of allele copies in a population, \eqn{a_i}{a_i} is the size of an allele
expressed in number of repeat units, \eqn{\bar{p}_i}{p_i} is an allele frequency averaged across
populations (with populations weighted equally), and \eqn{n}{n} is the total number of allele copies
across all populations.

}
\value{
If \code{global = FALSE} and \code{bootstrap = FALSE}, a square matrix containing \eqn{F_{ST}}{FST}, \eqn{G_{ST}}{GST}, \eqn{D}{D}, or \eqn{R_{ST}}{RST} values.  
The rows and columns of the matrix are both named by population.

If \code{global = TRUE} and \code{bootstrap = FALSE}, a single value indicating the \eqn{F_{ST}}{FST}, \eqn{G_{ST}}{GST}, \eqn{D}{D}, or \eqn{R_{ST}}{RST} value.

If \code{global = TRUE} and \code{bootstrap = TRUE}, a vector of bootstrapped replicates of \eqn{F_{ST}}{FST}, \eqn{G_{ST}}{GST}, \eqn{D}{D}, or \eqn{R_{ST}}{RST}.

If \code{global = FALSE} and \code{bootstrap = TRUE}, a square two-dimensional list, with rows and columns named by population.  Each item is
a vector of bootstrapped values for \eqn{F_{ST}}{FST}, \eqn{G_{ST}}{GST}, \eqn{D}{D}, or \eqn{R_{ST}}{RST} for that pair of populations.
}
\references{
Nei, M. (1973) Analysis of gene diversity in subdivided populations.
\emph{Proceedings of the National Academy of Sciences of the United
  States of America} \bold{70}, 3321--3323.
  
Nei, M. and Chesser, R. (1983) Estimation of fixation indices and gene diversities.
\emph{Annals of Human Genetics} \bold{47}, 253--259.
  
Jost, L. (2008) \eqn{G_{ST}}{GST} and its relatives to not measure differentiation.
\emph{Molecular Ecology} \bold{17}, 4015--4026.

Slatkin, M. (1995) A measure of population subdivision based on microsatellite 
allele frequencies.  \emph{Genetics} \bold{139}, 457--462.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{simpleFreq}}, \code{\link{deSilvaFreq}}
}
\examples{
# create a data set (typically done by reading files)
mygenotypes <- new("genambig", samples = paste("ind", 1:6, sep=""),
                   loci = c("loc1", "loc2"))
Genotypes(mygenotypes, loci = "loc1") <- list(c(206), c(208,210),
                                              c(204,206,210),
    c(196,198,202,208), c(196,200), c(198,200,202,204))
Genotypes(mygenotypes, loci = "loc2") <- list(c(130,134), c(138,140),
                                              c(130,136,140),
    c(138), c(136,140), c(130,132,136))
PopInfo(mygenotypes) <- c(1,1,1,2,2,2)
mygenotypes <- reformatPloidies(mygenotypes, output="sample")
Ploidies(mygenotypes) <- c(2,2,4,4,2,4)
Usatnts(mygenotypes) <- c(2,2)

# calculate allele frequencies
myfreq <- simpleFreq(mygenotypes)

# calculate pairwise differentiation statistics
myfst <- calcPopDiff(myfreq, metric = "Fst")
mygst <- calcPopDiff(myfreq, metric = "Gst")
myD <- calcPopDiff(myfreq, metric = "Jost's D")
myrst <- calcPopDiff(myfreq, metric = "Rst", object = mygenotypes)

# examine the results
myfst
mygst
myD
myrst

# get global statistics
calcPopDiff(myfreq, metric = "Fst", global = TRUE)
calcPopDiff(myfreq, metric = "Gst", global = TRUE)
calcPopDiff(myfreq, metric = "Jost's D", global = TRUE)
calcPopDiff(myfreq, metric = "Rst", global = TRUE, object = mygenotypes)
}
\keyword{ arith }
\keyword{ array }

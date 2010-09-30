Bruvo.distance <- function(genotype1, genotype2, maxl=9, usatnt=2, missing=-9) {

require(combinat) # for permn

if(identical(genotype1, genotype2)) {dist <- 0} else {
    if((length(genotype1)>maxl & length(genotype2)>maxl)|
       genotype1[1]==missing|genotype2[1]==missing){
        dist <- NA} else {
        if(length(genotype1) >= length(genotype2)) {
            genotypeL <- genotype1/usatnt; genotypeS <- genotype2/usatnt} else {
                genotypeL <- genotype2/usatnt; genotypeS <- genotype1/usatnt}

# if genotypes are identical, just return a distance of zero without doing the rest
        # of the calculation
# if genotypes are both longer than 9, skip this calculation because it will take hours
# whichever genotype has more alleles, make this genotypeL (long) and the other
        # genotypeS (short)
# convert alleles into repeat counts by dividing by usatnt

kl <- length(genotypeL) # sets the ploidy level for this genotype comparison
ks <- length(genotypeS) # number of alleles in the shorter genotype

allele.distances <- array(0 , c(kl,ks))
# Create an empty matrix to contain the raw distances between alleles

for(n in 1:kl) { for(m in 1:ks) {allele.distances[n,m] <- genotypeL[n] - genotypeS[m]}}
# fills the array with the differences in allele repeat count

geometric.distances <- array(1 - 2^-abs(allele.distances) , c(kl,ks))
# geometric transformation based on mutation probabilities

#Next, find the minimum distance sum among all permutations

column <- 1:ks # an index of all columns (genotypeS alleles)
row <- 1:kl # an index of all rows (genotypeL alleles)

combinations <- combn(row, ks, FUN = NULL, simplify=FALSE)
# all combinations of alleles in genotypeL that can be matched to non-virtual
        # alleles in genotypeS

permutations <- permn(ks)
# all possible orders that alleles within these combinations can go in

mindist <- Inf # this variable will store the minimum sum encountered so far.

for(i in 1:length(combinations)) {
# the loop to go through every possible sum of compatible allele comparisons

rowcomb <- combinations[[i]] # choose one combination of rows for this round

for(l in 1:length(permutations)){ # go through all orders of this combinations of rows

sum <- 0 # this is si, the sum of allele comparisons

for(j in 1:ks){
    sum <- sum + geometric.distances[rowcomb[permutations[[l]][j]],column[j]]}
# the loop to calculate the sum for this permutation

if(sum < mindist) {mindist <- sum} # is this the minimum sum found so far?

}}

dist <- (mindist+kl-ks)/kl
# add 1 for each infinite virtual allele, then divide by the ploidy

}}
return(dist)
}

Lynch.distance<-function(genotype1,genotype2,usatnt=NA,missing=-9){
    if(genotype1[1]==missing || genotype2[1]==missing){
        # return NA if there is any missing data
        distance<-NA
    } else {
        # get the average number of bands for the two genotypes
        meanbands<-(length(genotype1)+length(genotype2))/2
        # find how many bands the genotypes have in common
        commonbands<-length(genotype1[genotype1 %in% genotype2])
        # calculate the distance
        distance<- 1-(commonbands/meanbands)
    }

    # return distance
    return(distance)
}


meandistance.matrix <- function(object, samples=Samples(object),
                                 loci=Loci(object), all.distances=FALSE,
                                 distmetric=Bruvo.distance,
                                progress=TRUE, ...){
    # subset the object so that samples can be numbered
    object <- object[samples, loci]

# create an array containing all distances by locus and sample
    loci.matrices<-array(dim=c(length(loci),length(samples),length(samples)),
                         dimnames=list(loci,samples,samples))
    for(L in loci){
       for(m in 1:length(samples)){
           for(n in m:length(samples)){
               thisdistance <- distmetric(Genotype(object, m, L),
                                          Genotype(object, n, L),
                                          usatnt = Usatnts(object)[L],
                                          missing = Missing(object),
                                          ...)
               loci.matrices[L,m,n] <- thisdistance
               loci.matrices[L,n,m] <- thisdistance
               if(progress) print(c(L, samples[m], samples[n]))
           }
       }
    }

    # calculate the mean matrix across all loci
    mean.matrix<-matrix(nrow=length(samples),ncol=length(samples),
                        dimnames=list(samples,samples))
    for(j in samples){
        for(k in samples){
            mean.matrix[j,k]<-mean(loci.matrices[,j,k][!is.na(loci.matrices[,j,k])])
        }
    }

    # return either the mean matrix and possibly the array as well
    if(all.distances){
        return(list(DistByLoc=loci.matrices, MeanMatrix=mean.matrix))
    } else {
        return(mean.matrix)
    }
}

meandist.from.array<-function(distarray, samples=dimnames(distarray)[[2]],
                              loci=dimnames(distarray)[[1]]){
    # get the array to be averaged
    subarray<-distarray[loci,samples,samples]
    # make a matrix to put the means into
    mean.matrix<-matrix(nrow=length(samples),ncol=length(samples),
                        dimnames=list(samples,samples))
    for(j in samples){
        for(k in samples){
            mean.matrix[j,k]<-mean(subarray[,j,k][!is.na(subarray[,j,k])])
        }
    }
    return(mean.matrix)
}

find.na.dist<-function(distarray, samples=dimnames(distarray)[[2]],
                              loci=dimnames(distarray)[[1]]){
    # set up vectors for data frame to contain info on where missing data is
    Locus<-""
    Sample1<-""
    Sample2<-""
    # current row in the data frame
    currrow<-1
    # go through the array, find NA, and put the index into the data frame
    for(L in loci){
        for(s1 in samples){
            for(s2 in samples){
                if(is.na(distarray[L,s1,s2])){
                    Locus[currrow]<-L
                    Sample1[currrow]<-s1
                    Sample2[currrow]<-s2
                    currrow<-currrow+1
                }
            }
        }
    }
    #return data frame
    return(data.frame(Locus=Locus, Sample1=Sample1, Sample2=Sample2, stringsAsFactors=FALSE))
}

find.missing.gen<-function(object, samples=Samples(object),
                           loci=Loci(object)){
    # set up vectors to contain the indices to put into the data frame
    Locus <- c("")
    Sample <- c("")
    # current row in the data frame
    currrow<-1
    # find which data are missing
    for(L in loci){
        for(s in samples){
            if(isMissing(object, s, L)){
                Locus[currrow]<-L
                Sample[currrow]<-s
                currrow<-currrow+1
            }
        }
    }
    # return the data frame
    return(data.frame(Locus=Locus,Sample=Sample, stringsAsFactors=FALSE))
}

# find NA distances that aren't the result of missing data
find.na.dist.not.missing<-function(object, distarray,
                                   samples=dimnames(distarray)[[2]],
                                   loci=dimnames(distarray)[[1]]){
    # get the data frames of NA distances and missing genotypes
    na.dist<-find.na.dist(distarray,samples=samples,loci=loci)
    missing.gen<-find.missing.gen(object,samples=samples,loci=loci)
    # set up vectors for data frame to contain results
    Locus<-""
    Sample1<-""
    Sample2<-""
    currrow<-1
    # for each row of na.dist, look for that locus and samples in missing.gen
    for(i in 1:length(na.dist$Locus)){
        L<-na.dist$Locus[i]
        s1<-na.dist$Sample1[i]
        s2<-na.dist$Sample2[i]
        missthislocus<-missing.gen[missing.gen$Locus==L,]
        if(identical(c(s1,s2) %in% missthislocus$Sample, c(FALSE,FALSE))){
            Locus[currrow]<-L
            Sample1[currrow]<-s1
            Sample2[currrow]<-s2
            currrow<-currrow+1
        }
    }
    # return data frame
    return(data.frame(Locus=Locus, Sample1=Sample1, Sample2=Sample2,
                      stringsAsFactors=FALSE))
}

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

read.GeneMapper<-function(infiles, missing=-9){
    samples<-c()
    loci<-c()
    locusdata<-list()
    length(locusdata)<-length(infiles)
    for(i in 1:length(infiles)){
        locusdata[[i]]<-read.table(infiles[i],sep="\t",header=TRUE,
                                   as.is=c("Sample.Name","Marker"))
        samples<-c(samples,locusdata[[i]][["Sample.Name"]])
        loci<-c(loci,locusdata[[i]][["Marker"]])
    }
    samples<-unique(samples)
    loci<-unique(loci)
    genotypedata<-array(list(missing),c(length(samples), length(loci)),
                            list(samples,loci))
    for(m in 1:length(locusdata)){
        alleleindex<-grep("Allele",names(locusdata[[m]]),value=FALSE)
        for(j in 1:length(locusdata[[m]][["Sample.Name"]])){
                untrimmedgenotype<-locusdata[[m]][j,alleleindex]
                #Gives an error if there are loci or samples not in the arguments
                genotypedata[[locusdata[[m]][["Sample.Name"]][j],
                              locusdata[[m]][["Marker"]][j]]]<-
                    untrimmedgenotype[!is.na(untrimmedgenotype)]
            }
        }
        return(genotypedata)
}

distance.matrix.1locus<-function(gendata, distmetric=Bruvo.distance,
                                 progress=TRUE, ...){
    count<-length(gendata) # how many genotypes there are
    distances<-matrix(nrow=count, ncol=count,
                      dimnames=list(names(gendata),names(gendata)))
    for(m in 1:count){for(n in m:count){
        thisdistance<-distmetric(gendata[[m]],gendata[[n]], ...)
        distances[m,n]<-thisdistance
        distances[n,m]<-thisdistance
        if(progress){print(c(names(gendata)[[m]],names(gendata)[[n]]))}
    }}
    return(distances)
}

meandistance.matrix <- function(gendata, samples=dimnames(gendata)[[1]],
                                 loci=dimnames(gendata)[[2]], all.distances=FALSE,
                                 usatnts=NULL, ...){

# create an array containing all distances by locus and sample
    loci.matrices<-array(dim=c(length(loci),length(samples),length(samples)),
                         dimnames=list(loci,samples,samples))
    for(i in loci){
       if(is.null(usatnts)){
           loci.matrices[i,,]<-distance.matrix.1locus(gendata[samples,i],...)
       } else {
           loci.matrices[i,,]<-distance.matrix.1locus(gendata[samples,i],
                                                      usatnt=usatnts[i],...)
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

write.Structure <- function(gendata, ploidy, file="",
                            samples=dimnames(gendata)[[1]], loci=dimnames(gendata)[[2]],
                            indploidies=rep(ploidy,times=length(samples)),
                            extracols=NULL, missing=-9){
    if(is.null(names(indploidies))) {names(indploidies)<-samples}
    structdata <- data.frame(rowlabel=c("missing",rep(samples, each=ploidy)))
    thiscol = 1 # which column we are working with

    # Fill in extra columns.  Ignored if extracols=NULL.
    for(excol in dimnames(extracols)[[2]]){
        thiscol <- thiscol+1
        structdata[[thiscol]]<-c(0,rep(extracols[samples, excol], each=ploidy))
        names(structdata)[thiscol]<-excol
    }

    # Fill in the genotypes
    for(locus in loci){
        thiscol<-thiscol+1
        alleles<-missing #set up the vector of alleles
        for(s in samples){
            #If missing data must be inserted to reflect a lower ploidy level:
            if(indploidies[s] < ploidy){
                if(length(gendata[[s,locus]]) == indploidies[s]){
                    #fully heterozygous genotype
                    thesealleles <- c(gendata[[s,locus]],
                                    rep(missing, times=ploidy-indploidies[s]))
              } else {
                  if(length(gendata[[s,locus]]) < indploidies[s]){
                      #duplicate the first allele to get to the right ploidy
                      thesealleles <- c(gendata[[s,locus]],
                                        rep(gendata[[s,locus]][1],
                                            times=indploidies[s]-
                                            length(gendata[[s,locus]])),
                                        rep(missing, times=ploidy-indploidies[s]))
                  } else {
                      #randomly choose alleles to use if there are too many
                      thesealleles <- c(sample(gendata
                                       [[s,locus]],indploidies[s],replace=FALSE),
                                        rep(missing, times=ploidy-indploidies[s]))
                  }
              }
            #If the individual has equal or greater ploidy to that
                #being used in the file:
            } else {
                if(length(gendata[[s,locus]]) == ploidy){
                    #genotype fills available spaces
                    thesealleles <- gendata[[s,locus]]
                } else {
                    if(length(gendata[[s,locus]]) < ploidy){
                        #duplicate the first allele to get to the right ploidy
                        thesealleles<-c(gendata[[s,locus]],
                                        rep(gendata[[s,locus]][1],
                                            times=ploidy-length(gendata[[s,locus]])))
                    } else {
                        #randomly choose alleles to use if there are too many
                        thesealleles<-sample(gendata[[s,locus]],ploidy,replace=FALSE)
                    }
                }
            }
            alleles<-c(alleles,thesealleles)
        }
        structdata[[thiscol]]<-alleles
        names(structdata)[thiscol]<-locus
    }
   write.table(structdata, file=file, sep="\t", row.names=FALSE,
               quote=FALSE, col.names=TRUE)
}

estimate.ploidy <- function(gendata, samples=dimnames(gendata)[[1]],
                            loci=dimnames(gendata)[[2]]){
    # set up array to contain the maximum and average number of alleles
    ploidyinfo <- array(dim=c(length(samples),2),
                        dimnames=list(samples, c("max.alleles","mean.alleles")))

    # fill the array
    for(s in samples){
        numalleles<-mapply(length,gendata[s,loci])
        ploidyinfo[s,"max.alleles"] <- max(numalleles)
        ploidyinfo[s,"mean.alleles"] <- mean(numalleles)
    }

    return(ploidyinfo)
}

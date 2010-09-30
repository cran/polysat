genambig.to.genbinary <- function(object, samples=Samples(object), loci=Loci(object)){
    # set up the object that will ultimately be returned
    objectN <- new("genbinary", samples, loci)

    # fill in the slots that will be identical
    PopNames(objectN) <- PopNames(object)
    PopInfo(objectN) <- PopInfo(object)[samples]
    Ploidies(objectN) <- Ploidies(object)[samples]
    Usatnts(objectN) <- Usatnts(object)[loci]
    Description(objectN) <- Description(object)

    # find all unique alleles for each locus
    locvector<-c()
    allelevector<-c()
    for(L in loci){
        thesealleles <- c()
        for(s in samples){
            if(!isMissing(object, s, L)) thesealleles <-
                c(thesealleles, Genotype(object, s, L))
        }
        thesealleles<-sort(unique(thesealleles))
        locvector <- c(locvector, rep(L, length(thesealleles)))
        allelevector <- c(allelevector, thesealleles)
    }

    # Build data frame of locus and allele information
    colinfo <- data.frame(Loci=locvector, Alleles=allelevector,
                          stringsAsFactors=FALSE)
    # Set up matrix to ultimately put in the Genotypes slot
    domdata <- matrix(nrow = length(samples), ncol=dim(colinfo)[1],
                      dimnames=list(samples, paste(locvector, allelevector,
                      sep=".")))
    # Fill the matrix
    for(m in 1:dim(colinfo)[1]){
        for(s in samples){
            if(isMissing(object, s, colinfo[m,1])){
                domdata[s,m] <- Missing(objectN)
            } else {
                if(colinfo[[m,2]] %in% Genotype(object, s, colinfo[m,1])){
                    domdata[s,m] <- Present(objectN)
                } else { domdata[s,m] <- Absent(objectN) }
            }
        }
    }

    # Fill in genotypes slot
    Genotypes(objectN) <- domdata

    # return the converted object
    return(objectN)
}

genbinary.to.genambig <- function(object, samples = Samples(object),
                                  loci = Loci(object)){
    # set up new genambig object
    objectN <- new("genambig", samples, loci)

    # fill in the slots that will be identical
    PopNames(objectN) <- PopNames(object)
    PopInfo(objectN) <- PopInfo(object)[samples]
    Ploidies(objectN) <- Ploidies(object)[samples]
    Usatnts(objectN) <- Usatnts(object)[loci]
    Description(objectN) <- Description(object)

    # Go through matrix one locus at a time to get genotypes
    for(L in loci){
        thesegen <- as.matrix(Genotypes(object, samples, L))
        # Go through each allele
        for(n in 1:dim(thesegen)[2]){
            allele <- strsplit(dimnames(thesegen)[[2]][n], split=".",
                               fixed=TRUE)[[1]][2]
            allele <- as.integer(allele)
            # Look for allele in each sample and write to new object
            for(s in samples){
                if(thesegen[s,n] == Present(object)){
                    if(isMissing(objectN, s, L)){
                        Genotype(objectN, s, L) <- allele
                    } else {
                        Genotype(objectN, s, L) <- c(Genotype(objectN, s, L),
                                                     allele)
                    }
                }
            }
        }
    }

    # return the new genambig object
    return(objectN)
}

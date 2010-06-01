Lynch.distance<-function(genotype1,genotype2,missing=-9){
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

find.missing.gen<-function(gendata, missing=-9, samples=dimnames(gendata)[[1]],
                           loci=dimnames(gendata)[[2]]){
    # set up vectors to contain the indices to put into the data frame
    Locus <- c("")
    Sample <- c("")
    # current row in the data frame
    currrow<-1
    # find which data are missing
    for(L in loci){
        for(s in samples){
            if(gendata[[s,L]][1]==missing){
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
find.na.dist.not.missing<-function(gendata, distarray, missing=-9,
                                   samples=dimnames(distarray)[[2]],
                                   loci=dimnames(distarray)[[1]]){
    # get the data frames of NA distances and missing genotypes
    na.dist<-find.na.dist(distarray,samples=samples,loci=loci)
    missing.gen<-find.missing.gen(gendata,missing=missing,samples=samples,loci=loci)
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
    return(data.frame(Locus=Locus, Sample1=Sample1, Sample2=Sample2, stringsAsFactors=FALSE))
}

estimate.freq<-function(gendata, missing=-9, samples=dimnames(gendata)[[1]],
                        loci=dimnames(gendata)[[2]], popinfo=rep(1,length(samples)),
                        indploidies=rep(4,length(samples))){
    # name popinfo and indploidies if necessary
    if(is.null(names(popinfo))){
        names(popinfo)<-samples
    }
    if(is.null(names(indploidies))){
        names(indploidies)<-samples
    }
    # get presence/absence data for alleles, and the table of which allele and locus
    # is in which column
    gentable<-codominant.to.dominant(gendata[samples,loci], makecolinfo=TRUE,
                                     missingin=missing)
    colinfo<-gentable[[2]]
    gentable<-gentable[[1]]
    # get a list of populations
    pops<-sort(unique(popinfo[samples]))
    # get total number of genomes per population
    totgenomes<-rep(0,length(pops))
    names(totgenomes)<-pops
    for(p in pops){
        totgenomes[p]<-sum(indploidies[samples[popinfo[samples]==p]])
    }

    # set up data frame to contain allele frequencies
    freqtable<-data.frame(Genomes=totgenomes,row.names=pops)

    # loop to get frequency data
    for(L in loci){
        # get all samples without missing data at this locus
        xsamples<-samples[gentable[,match(L,colinfo[[1]])]!=-9]
        # get the total number of genomes per population
        totgenomes<-rep(0,length(pops))
        names(totgenomes)<-pops
        for(p in pops){
            totgenomes[p]<-sum(indploidies[xsamples[popinfo[xsamples]==p]])
        }
        # make a conversion factor to weight allele presence based on ploidy
        # of each individual and number of alleles at this locus
        numalleles<-sapply(gendata[xsamples,L],length)
        names(numalleles)<-xsamples
        convf<-indploidies[xsamples]/numalleles
        # convert alleles in the table to estimated copy number
        loctable<-gentable[xsamples,colinfo[[1]]==L]*convf
        # loop through alleles at this locus
        for(al in colinfo[[2]][colinfo[[1]]==L]){
            theseallelefreqs<-rep(0,length(pops))
            names(theseallelefreqs)<-pops
            # loop through populations
            for(p in pops){
                theseallelefreqs[p]<-sum(loctable[xsamples[popinfo[xsamples]==p],
                                                  paste(L,al,sep=".")])/totgenomes[p]
            }
            freqtable<-cbind(freqtable,theseallelefreqs)
        }
    }

    # return data frame
    names(freqtable)<-c("Genomes",dimnames(gentable)[[2]])
    return(freqtable)
}

calcFst<-function(freqs, pops=row.names(freqs),
                  loci=unique(as.matrix(as.data.frame(strsplit(names(freqs), split=".",
                                                      fixed=TRUE),
                                             stringsAsFactors=FALSE))[1,])){
    # Clean up loci
    loci<-loci[loci!="Genomes"]
    # Set up matrix for Fst values
    fsts<-matrix(0,nrow=length(pops),ncol=length(pops),dimnames=list(pops,pops))
    # Get genome number from the table
    genomes<-freqs$Genomes
    names(genomes)<-pops
    for(m in 1:length(pops)){
        for(n in m:length(pops)){
            # set up array for HT and HS values
            hets<-array(0,dim=c(length(loci),2),dimnames=list(loci,c("HT","HS")))
            for(L in loci){
                # get just the frequencies for these pops and this locus
                thesefreqs<-freqs[c(pops[m],pops[n]),grep(L,names(freqs),fixed=TRUE)]
                # get average allele frequencies weighted by genomes/pop
                avgfreq<-(thesefreqs[1,]*genomes[pops[m]] + thesefreqs[2,]*genomes[pops[n]])/
                    (genomes[pops[m]] + genomes[pops[n]])
                # estimate H by 1 - sum of squared allele frequencies
                # put the heterozygositites in the array
                hets[L,"HT"]<-1-sum(avgfreq^2)
                hets[L,"HS"]<-((1-sum(thesefreqs[1,]^2))*genomes[pops[m]] +
                               (1-sum(thesefreqs[2,]^2))*genomes[pops[n]])/
                                   (genomes[pops[m]] + genomes[pops[n]])
            }
            HT<-mean(hets[,"HT"])
            HS<-mean(hets[,"HS"])
            fsts[m,n]<-(HT-HS)/HT
            fsts[n,m]<-(HT-HS)/HT
        }
    }
    # return matrix of Fst values
    return(fsts)
}


# join two distance arrays


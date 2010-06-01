read.SPAGeDi<-function(infile, allelesep="/", returncatspatcoord=FALSE,
                       returnploidies=FALSE, missing=-9){
    # get all lines from the file, and eliminate comment lines
    Lines<-readLines(infile)
    Lines<-Lines[Lines != ""]
    first2<-sapply(Lines, substring, first=1, last=2)
    Lines<-Lines[first2 != "//"]
    # get data from the first line
    fileinfo<-as.integer(strsplit(Lines[1], "\t")[[1]])
    numind<-fileinfo[1] # number of samples
    numcat<-fileinfo[2] # number of categories
    numsc<-fileinfo[3] # number of spatial coordinates
    numloc<-fileinfo[4] # number of loci
    digits<-fileinfo[5] # number of digits to represent alleles

    # Is there a column of categories?
    if(numcat==0){
        catpres<-0
    } else {
        catpres<-1
    }
    # Is latitude and longitude used instead of Cartesian coordinates?
    if(numsc==-2){
        numsc<-2
    }

    # read the rest of the file as a table
    cat(Lines[3:(3+numind)], sep="\n", file="SpagTemp.txt")
    gentable <- read.table("SpagTemp.txt", sep="\t", header=TRUE,
                           row.names=1,
                           colClasses=c("character",rep(NA,catpres+numsc),
                           rep("character",numloc)))
    # get sample and locus names
    samples<-row.names(gentable)
    loci<-names(gentable)[(length(gentable)-numloc+1):length(gentable)]

    # set up list to contain genotypes
    gendata <- array(list(missing), dim=c(length(samples), length(loci)),
                     dimnames=list(samples,loci))
    # set up vector to contain ploidies
    indploidies <- rep(4, length(samples))
    names(indploidies) <- samples

    # If there is no separation of alleles, count digits off with substring
    if(allelesep==""){
        for(s in samples){
            #set up a list to contain genotypes
            thesegenotypes<-list(0)
            length(thesegenotypes)<-length(loci)
            names(thesegenotypes)<-loci

            for(L in loci){
                # add leading zeros if necessary
                while(nchar(gentable[s,L])%%digits !=0){
                    gentable[s,L]<-paste("0",gentable[s,L],sep="")
                }
                # split into alleles and convert to integers
                thesegenotypes[[L]]<-as.integer(substring(gentable[s,L],
                                               first=seq(1,
                                               nchar(gentable[s,L])-digits+1,
                                               by=digits),
                                               last=seq(digits,
                                               nchar(gentable[s,L]),by=digits)))
                # if genotype only has zeros, write missing data symbol
                if(length(unique(thesegenotypes[[L]]))==1 && thesegenotypes[[L]][1]==0){
                    thesegenotypes[[L]]<-missing
                }
                # otherwise remove zeros on left
                while(thesegenotypes[[L]][1]==0){
                    thesegenotypes[[L]]<-thesegenotypes[[L]][-1]
                }
            }
            # get ploidy of sample
            indploidies[s] <- max(sapply(thesegenotypes,length))

            for(L in loci){
                # remove zeros on the right
                thesegenotypes[[L]]<-thesegenotypes[[L]][thesegenotypes[[L]] != 0]
                # get unique alleles
                thesegenotypes[[L]]<-unique(thesegenotypes[[L]])
            }
            # add genotypes to list
            gendata[s,]<-thesegenotypes
        }

    # get alleles by strsplit
    } else {
        for(s in samples){
            # get alleles by splitting the strings
            thesegenotypes<-sapply(gentable[s,loci],strsplit,
                                   split=allelesep,fixed=TRUE)
            names(thesegenotypes) <- loci

            for(L in loci){
                # convert to integer
                thesegenotypes[[L]]<-as.integer(thesegenotypes[[L]])
                # if genotype only has zeros, write missing data symbol
                if(length(unique(thesegenotypes[[L]]))==1 && thesegenotypes[[L]][1]==0){
                    thesegenotypes[[L]]<-missing
                }
                # otherwise remove zeros on left
                while(thesegenotypes[[L]][1]==0){
                    thesegenotypes[[L]]<-thesegenotypes[[L]][-1]
                }
            }
            # get ploidy of sample
            indploidies[s] <- max(sapply(thesegenotypes,length))

            for(L in loci){
                # remove zeros on the right
                thesegenotypes[[L]]<-thesegenotypes[[L]][thesegenotypes[[L]] != 0]
                # get unique alleles
                thesegenotypes[[L]]<-unique(thesegenotypes[[L]])
            }
            # add genotypes to list
            gendata[s,]<-thesegenotypes
        }
    }

    # return the genotypes, and other data as applicable
    if(identical(c(returncatspatcoord,returnploidies), c(FALSE,FALSE))){
        return(gendata)
    } else {
        return(list(CatSpatCoord=gentable[,!names(gentable) %in% loci],
                    Indploidies=indploidies,
                    Genotypes=gendata)[c(returncatspatcoord,
                    returnploidies,TRUE)])
    }
}

write.SPAGeDi<-function(gendata,samples=dimnames(gendata)[[1]],
                        loci=dimnames(gendata)[[2]],
                        indploidies=rep(4,length(samples)),
                        popinfo=rep(1,length(samples)),
                        allelesep="/", digits=2, file="",
                        spatcoord=data.frame(X=rep(1,length(samples)),
                                             Y=rep(1,length(samples)),
                                             row.names=samples),
                        usatnts=rep(2, length(loci)),
                        missing=-9){
    # name indploidies, popinfo, usatnts, and if not already done
    if(is.null(names(indploidies))){names(indploidies)<-samples}
    if(is.null(names(popinfo))){names(popinfo)<-samples}
    if(is.null(names(usatnts))){names(usatnts)<-loci}
    if(identical(row.names(spatcoord), as.character(1:dim(spatcoord)[1]))){
        row.names(spatcoord)<-samples}

    # set up data frame to contain genotypes
    gentable<-data.frame(Ind=samples, Cat=popinfo[samples], spatcoord[samples,])
    # find missing data
    misstable<-find.missing.gen(gendata[samples,loci], missing=missing)
    # get genotype data by locus (column)
    for(L in loci){
        genotypesL<-gendata[samples,L]
        missingsamples<-misstable$Sample[misstable$Locus==L]
        # replace missing data with zeros
        genotypesL[missingsamples]<-0
        # divide and subtract to convert to repeats
        divallele<-function(value){floor(value/usatnts[L])}
        genotypesL<-mapply(divallele,genotypesL,SIMPLIFY=FALSE)
        suballele<-function(value){if(value[1]!=0){
            value-(10^(digits-1))} else{
            0}}
        while(max(mapply(max,genotypesL)) >= 10^digits){
            genotypesL<-mapply(suballele, genotypesL, SIMPLIFY=FALSE)
        }
        # add zeros up to correct ploidy, delete alleles if necessary
        zerostoadd<-indploidies[samples] - mapply(length,genotypesL)
        names(zerostoadd)<-samples
        for(s in samples){
            if(length(genotypesL[[s]])==1){
                # replicate allele if totally homozygous
                genotypesL[[s]]<-rep(genotypesL[[s]],indploidies[s])
            } else {
                if(zerostoadd[s] < 0){
                    # randomly remove alleles if there are too many
                    genotypesL[[s]]<-sample(genotypesL[[s]],indploidies[s],
                                            replace=FALSE)
                    cat("Alleles randomly removed to get to ploidy:",L,s,"\n")
                } else {
                    # add zeros for partial heterozygotes
                    genotypesL[[s]]<-c(genotypesL[[s]],rep(0,zerostoadd[s]))
                }
            }
            # also make each allele the right number of digits if necessary
            if(allelesep==""){
                genotypesL[[s]]<-as.character(genotypesL[[s]])
                for(a in 1:length(genotypesL[[s]])){
                    while(nchar(genotypesL[[s]][a]) < digits){
                        genotypesL[[s]][a]<-paste(0,genotypesL[[s]][a],sep="")
                    }
                }
            }
        }
        # concatenate into strings
        genvect<-mapply(paste,genotypesL,collapse=allelesep)
        # add the vector to the data frame
        gentable<-data.frame(gentable,genvect)
        names(gentable)[dim(gentable)[2]]<-L
    }
    # write file
    write.table(gentable,file="SpagTemp.txt",sep="\t",row.names=FALSE,
                col.names=TRUE,quote=FALSE)
    cat(paste(length(samples),length(unique(popinfo[samples])),dim(spatcoord)[2],
                length(loci),digits,max(indploidies[samples]),sep="\t"),"0",
        readLines("SpagTemp.txt"),"END",sep="\n",file=file)
}

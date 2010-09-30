read.GeneMapper<-function(infiles){
    # read files and get names of samples and loci
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

    # set up genambig object to put data into
    object <- new("genambig", samples, loci)

    # extract and insert genotypes
    for(m in 1:length(locusdata)){
        alleleindex<-grep("Allele",names(locusdata[[m]]),value=FALSE)
        for(j in 1:length(locusdata[[m]][["Sample.Name"]])){
                untrimmedgenotype<-locusdata[[m]][j,alleleindex]
                #Gives an error if there are loci or samples not in the arguments
                Genotype(object,locusdata[[m]][j, "Sample.Name"],
                              locusdata[[m]][j, "Marker"])<-
                    unique(untrimmedgenotype[!is.na(untrimmedgenotype)])
        }
    }

    return(object)
}

read.GenoDive <- function(infile){
    rawdata<-readLines(infile)
    #information about number of samples, number of loci, etc. is in second line
    datainfo<-as.integer(strsplit(rawdata[[2]],"\t")[[1]])
    numsam<-datainfo[1] # number of individuals
    numpop<-datainfo[2] # number of populations
    numloc<-datainfo[3] # number of loci
    digits<-datainfo[5] # number of digits used to code each allele

    # get the data header with column names, including loci
    colheader<-strsplit(rawdata[[3+numpop]],"\t")[[1]]
    lastloc<-length(colheader) # index of the last locus column
    firstloc<-lastloc-numloc+1 # index of the first locus column (3 or 4)
    loci<-colheader[firstloc:lastloc] # get the locus names

    # set up genotype object
    object <- new("genambig", samples=1:numsam, loci = loci)
    PopNames(object) <- rawdata[3:(2+numpop)]
    Description(object) <- rawdata[1]
    # fill in sample data
    for(s in 1:numsam){
        # get data from the row for this sample
        sampledata<-strsplit(rawdata[[3+numpop+s]],"\t")[[1]]
        # extract name and population
        PopInfo(object)[s]<-as.integer(sampledata[1])
        Samples(object)[s]<-sampledata[2]
        # extract alleles
        for(L in 1:numloc){
            rawalleles<-sampledata[(firstloc:lastloc)[L]]
            while(nchar(rawalleles)%%digits !=0) {
                # If leading zeros were removed
                rawalleles<-paste("0",rawalleles,sep="")
            }
            # convert the character string to an integer vector
            myfirst<-seq(length=nchar(rawalleles)/digits,from=1,by=digits)
            mylast<-seq(length=nchar(rawalleles)/digits,from=digits,by=digits)
            thesealleles<-as.integer(substring(rawalleles,myfirst,mylast))
            # get rid of duplicate alleles
            thesealleles<-unique(thesealleles)
            # insert the missing data symbol as appropriate
            if(length(thesealleles)==1 && thesealleles==0){
                thesealleles<-Missing(object)
            }
            # get rid of "missing alleles" if the whole genotype is not missing
            thesealleles<-thesealleles[thesealleles !=0]
            # write the allele vector to the list
            Genotype(object, s, L) <-thesealleles
        }
    }
    # return the object
    return(object)
}

# make extracols include the popinfo column
# popinfocol is column number not including labels, or NA
read.Structure<-function(infile,ploidy,missingin=-9,sep="\t",markernames=TRUE,
                         labels=TRUE, extrarows=1, popinfocol=1, extracols=1,
                         getexcols=FALSE
                         ){
    # read the file
    rawdata<-read.table(infile,header=markernames,sep=sep)
    # get an index of samples and a column labeling the samples
    if(labels){
        # if row labels are used, get those as the sample names
        samples<-unique(rawdata[[1]])
        samples<-samples[!is.na(samples)]
        samples<-as.character(samples)
        samples<-samples[samples != ""]
        names(rawdata)[1]<-"Samples"
    } else {
        # make an integer vector to represent samples
        samples<-1:((dim(rawdata)[1]-extrarows)/ploidy)
        samindex<-c(rep(0,times=extrarows),rep(samples, each=ploidy))
        rawdata[length(rawdata)+1]<-samindex
        names(rawdata)[length(rawdata)]<-"Samples"
    }
    # get an index of loci (will be V2 etc. if loci not named)
    loci<-names(rawdata)
    loci<-loci[loci != "Samples"]
    loci<-loci[(extracols+1):length(loci)]
    # mark the column containing popinfo
    if(!is.na(popinfocol)){
        names(rawdata)[popinfocol+ifelse(labels,1,0)] <- "PopInfo"
    }

    # set up the object to store genotypes
    object <- new("genambig", samples, loci)

    # fill the genotypes and popinfo
    for(s in samples){
        for(L in loci){
            rawalleles<-rawdata[rawdata$Samples==s,L]
            # process missing data and get unique alleles
            thesealleles<-unique(rawalleles)
            if(length(thesealleles)==1 && thesealleles[1]==missingin){
                thesealleles<-Missing(object)
            } else {
                thesealleles<-thesealleles[thesealleles != missingin]
                if(is.na(Ploidies(object)[s]) || Ploidies(object)[s] <
                   length(rawalleles[rawalleles != missingin])){
                    # get ploidy
                    Ploidies(object)[s] <-
                        length(rawalleles[rawalleles != missingin])
                }
            }
            Genotype(object,s,L)<-thesealleles
        }
        # get popinfo
        if(!is.na(popinfocol)){
            PopInfo(object)[s] <- rawdata[match(s, rawdata$Samples),"PopInfo"]
        }
    }

    # extract the extra columns, if needed
    if(getexcols){
        Extracol<-data.frame(row.names=samples)
        excolindex <- 1:extracols
        if(!is.na(popinfocol)){
            excolindex <- excolindex[-popinfocol]
        }
        if(labels){
            excolindex <- excolindex + 1
        }
        for(x in excolindex){
            # look up values by sample in this column
            # add the column to the data frame
            Extracol<-data.frame(Extracol,rawdata[[x]][seq(1+extrarows,
                                                        length(rawdata[[1]])+
                                                           1-ploidy,
                                                        by=ploidy)])
            names(Extracol)[length(Extracol)] <- x
        }
    }

    # return extra columns and genotypes
    if(getexcols){
        return(list(ExtraCol=Extracol, Dataset=object))
    } else {
        return(object)
    }
}

read.SPAGeDi<-function(infile, allelesep="/", returnspatcoord=FALSE){
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

    # set up genambig object
    object <- new("genambig", samples=samples, loci=loci)
    # take PopInfo from categories if present
    if(catpres==1){
        PopNames(object) <- unique(as.character(gentable[[1]]))
        PopInfo(object) <- match(gentable[[1]], PopNames(object))
    }

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
                    thesegenotypes[[L]]<-Missing(object)
                }
                # otherwise remove zeros on left
                while(thesegenotypes[[L]][1]==0){
                    thesegenotypes[[L]]<-thesegenotypes[[L]][-1]
                }
            }
            # get ploidy of sample
            Ploidies(object)[s] <- max(sapply(thesegenotypes,length))

            for(L in loci){
                # remove zeros on the right
                thesegenotypes[[L]]<-thesegenotypes[[L]][thesegenotypes[[L]] != 0]
                # get unique alleles
                thesegenotypes[[L]]<-unique(thesegenotypes[[L]])
            }
            # add genotypes to object
            Genotypes(object, samples=s)<-thesegenotypes
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
                    thesegenotypes[[L]]<-Missing(object)
                }
                # otherwise remove zeros on left
                while(thesegenotypes[[L]][1]==0){
                    thesegenotypes[[L]]<-thesegenotypes[[L]][-1]
                }
            }
            # get ploidy of sample
            Ploidies(object)[s] <- max(sapply(thesegenotypes,length))

            for(L in loci){
                # remove zeros on the right
                thesegenotypes[[L]]<-thesegenotypes[[L]][thesegenotypes[[L]] != 0]
                # get unique alleles
                thesegenotypes[[L]]<-unique(thesegenotypes[[L]])
            }
            # add genotypes to list
            Genotypes(object, samples=s)<-thesegenotypes
        }
    }

    # return the genotypes, popinfo, ploidies, and spatial coordinates
    if(!returnspatcoord){
        return(object)
    } else {
        spatcoord <- gentable[,!names(gentable) %in% loci]
        if(catpres == 1){
            spatcoord <- spatcoord[,-1]
        }
        return(list(SpatCoord=spatcoord,
                    Dataset=object))
    }
}

read.Tetrasat <- function(infile){
    #read the file into a character vector, containing all the lines
    rawdata<-readLines(infile)
    #find which lines delimit populations
    popindex<-grep("pop",rawdata,ignore.case=TRUE,value=FALSE)
    #get a character vector of loci, whether they were stored in one or
    #several lines
    loci<-rawdata[2:(popindex[1]-1)]
    if(length(loci) == 1){loci<-strsplit(loci,",")[[1]]}
    #find which lines contain genotype data
    samindex<-c(popindex[1]:length(rawdata))
    samindex<-samindex[!samindex %in% popindex]
    # set up the genambig object
    object<-new("genambig", samples=1:length(samindex) ,loci=loci)
    Ploidies(object) <- rep(4, length(samindex))
    Description(object) <- rawdata[1]

    #Extract the data out of the lines
    for(i in 1:length(samindex)){
        #Find the largest value in popindex that is smaller than samindex[i]
        #Use the position of that value in popindex for the pop id
        PopInfo(object)[i]<-match(max(popindex[popindex < samindex[i]]),
                                  popindex)
        #extract sample names and genotypes
        samname<-gsub(" ","",substring(rawdata[samindex[i]],1,20),fixed=TRUE)
        Samples(object)[i] <- samname
        thesegenotypes<-substring(rawdata[samindex[i]],
                                  seq(length=length(loci),from=21,by=9),
                                  seq(length=length(loci),from=28,by=9))
        for(j in 1:length(loci)){
            thesealleles<-gsub(" ","",thesegenotypes[j],fixed=TRUE)
            if(nchar(thesealleles) != 0){
              thesealleles<-substring(thesealleles,
                                      seq(length=nchar(thesealleles)/2,from=1,by=2),
                                      seq(length=nchar(thesealleles)/2,from=2,by=2))
              thesealleles<-unique(as.integer(thesealleles))
            } else {thesealleles<-Missing(object)}
            Genotype(object, i,j)<-thesealleles
        }
    }
    #return genotype and pop data
    return(object)
}

read.ATetra<-function(infile){
    #read the file into a list of vectors, one for each line
    rawdata<-strsplit(readLines(infile),",")
    #make an index of which lines contain locus and sample info
    #also make an index of all loci, and vector of population names
    locindex<-c()
    samindex<-c()
    popnames<-character(0)
    for(i in 1:length(rawdata)){
        if(rawdata[[i]][1] == "TIT"){
            description <- rawdata[[i]][2]
        }
        if(rawdata[[i]][1] == "LOC"){
            locindex<-c(locindex,i)
            names(locindex)[length(locindex)]<-rawdata[[i]][3]
        }
        if(rawdata[[i]][1] == "IND"){
            samindex<-c(samindex,i)
        }
        if(rawdata[[i]][1] == "POP" && length(locindex) == 1){
            popnames <- c(popnames, rawdata[[i]][4])
        }
    }
    #make a vector of all loci
    loci<-names(locindex)
    #make a vector of all samples, and a vector of which samples go in which pops
    samples<-c()
    popdata<-c()
    samindex1loc<-samindex[samindex < locindex[2]]
    for(j in 1:length(samindex1loc)){
        samples[j]<-rawdata[[samindex1loc[j]]][5]
        popdata[j]<-as.integer(rawdata[[samindex1loc[j]]][3])
    }
    #set up the object to contain genotypes
    object <- new("genambig", samples=samples, loci=loci)
    PopInfo(object) <- popdata
    PopNames(object) <- popnames
    Ploidies(object) <- rep(4, length(samples))
    Description(object) <- description

    #fill the array of genotypes
    for(m in samindex){
        thesealleles<-rawdata[[m]][6:9]
        thesealleles<-as.integer(thesealleles[thesealleles !=""])
        thesealleles<-thesealleles[!is.na(thesealleles)]
        Genotype(object,as.integer(rawdata[[m]][4]),
                 as.integer(rawdata[[m]][2]))<-thesealleles
    }
    #return population data and genotypes
    return(object)
}

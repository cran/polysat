read.Tetrasat <- function(infile,missing=-9){
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
    #make a two-dimensional list to hold the genotypes
    gendata<-array(list(missing),dim=c(length(samindex),length(loci)),
                   dimnames=list(c(1:length(samindex)),loci))
    #make a vector containing population data
    popdata<-c()
    length(popdata)<-length(samindex)

    #Extract the data out of the lines
    for(i in 1:length(samindex)){
        #Find the largest value in popindex that is smaller than samindex[i]
        #Use the position of that value in popindex for the pop id
        popdata[i]<-match(max(popindex[popindex < samindex[i]]),popindex)
        #extract sample names and genotypes
        samname<-gsub(" ","",substring(rawdata[samindex[i]],1,20),fixed=TRUE)
        names(popdata)[i]<-samname
        dimnames(gendata)[[1]][i]<-samname
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
            } else {thesealleles<-missing}
            gendata[[i,j]]<-thesealleles
        }
    }
    #return genotype and pop data
    return(list(PopData=popdata,Genotypes=gendata))
}

read.ATetra<-function(infile){
    #read the file into a list of vectors, one for each line
    rawdata<-strsplit(readLines(infile),",")
    #make an index of which lines contain locus and sample info
    #also make an index of all loci
    locindex<-c()
    samindex<-c()
    for(i in 1:length(rawdata)){
        if(rawdata[[i]][1] == "LOC"){
            locindex<-c(locindex,i)
            names(locindex)[length(locindex)]<-rawdata[[i]][3]
        }
        if(rawdata[[i]][1] == "IND"){
            samindex<-c(samindex,i)
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
    names(popdata)<-samples
    #set up the array to contain genotypes
    gendata<-array(list(missing),dim=c(length(samples),length(loci)),
                   dimnames=list(samples,loci))
    #fill the array of genotypes
    for(m in samindex){
        thesealleles<-rawdata[[m]][6:9]
        thesealleles<-as.integer(thesealleles[thesealleles !=""])
        thesealleles<-thesealleles[!is.na(thesealleles)]
        gendata[[as.integer(rawdata[[m]][4]),as.integer(rawdata[[m]][2])]]<-thesealleles
    }
    #return population data and genotypes
    return(list(PopData=popdata,Genotypes=gendata))
}

dominant.to.codominant<-function(domdata,colinfo=NULL,samples=dimnames(domdata)[[1]],
                                 missing=-9,allelepresent=1,split="."){
    #get locus and allele information
    if(is.null(colinfo)){
        rawcol<-strsplit(dimnames(domdata)[[2]],split,fixed=TRUE)
        loc<-c()
        alleles<-c()
        for(i in 1:length(rawcol)){
            loc[i]<-rawcol[[i]][1]
            alleles[i]<-as.integer(rawcol[[i]][2])
        }
        colinfo<-data.frame(Loci=loc, Alleles=alleles)
    }

    #extract a list of loci
    loci<-unique(colinfo[[1]])
    #set up the genotype list
    gendata<-array(list(missing),dim=c(length(samples),length(loci)),
                   dimnames=list(samples,loci))
    #fill the genotype vectors
    for(n in 1:dim(domdata)[2]){
        for(m in samples){
            if(domdata[m,n] == allelepresent){
                gendata[[m,colinfo[n,1]]]<-c(gendata[[m,colinfo[n,1]]],colinfo[n,2])
                gendata[[m,colinfo[n,1]]]<-
                    gendata[[m,colinfo[n,1]]][gendata[[m,colinfo[n,1]]] != missing]
            }
        }
    }
    #return codominant genotypes
    return(gendata)
}

codominant.to.dominant<-function(gendata, makecolinfo=FALSE, allelepresent=1, alleleabsent=0,
                                 missingin=-9, missingout=-9,
                                 loci=dimnames(gendata)[[2]], samples=dimnames(gendata)[[1]]){
    # Find all unique alleles for each locus
    locvector<-c()
    allelevector<-c()
    for(l in loci){
        thesealleles<-c()
        for(s in samples){
            thesealleles<-c(thesealleles,gendata[[s,l]])
        }
        thesealleles<-unique(thesealleles)
        thesealleles<-thesealleles[thesealleles != missingin]
        thesealleles<-sort(thesealleles)
        locvector<-c(locvector,rep(l,length(thesealleles)))
        allelevector<-c(allelevector,thesealleles)
    }
    # Build data frame of locus and allele information.
    colinfo<-data.frame(Loci=locvector,Alleles=allelevector)
    # Set up matrix of dominant data
    domdata<-matrix(nrow=length(samples), ncol=dim(colinfo)[1],
                    dimnames=list(samples,1:dim(colinfo)[1]))
    # Search for alleles and fill matrix with presence and absence symbols
    for(m in 1:dim(colinfo)[1]){
        for(s in samples){
            if(gendata[[s,colinfo[m,1]]][1] == missingin){
                domdata[s,m]<-missingout
            } else {
            if(colinfo[[m,2]] %in% gendata[[s,colinfo[m,1]]]){
                domdata[s,m]<-allelepresent
            } else { domdata[s,m]<-alleleabsent } }
        }
        # Create the locus.allele name for the column
        dimnames(domdata)[[2]][m]<-paste(colinfo[m,1],colinfo[m,2], sep=".")
    }

    if(makecolinfo) return(list(Domdata=domdata, Colinfo=colinfo))
    else return(domdata)
}

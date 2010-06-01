read.GenoDive <- function(infile,missing=-9){
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

    # set up the list to contain the genotype data and the vector to contain the
    # population info
    gendata<-array(list(missing),dim=c(numsam,numloc),dimnames=list(1:numsam,loci))
    popinfo<-rep(0,times=numsam)
    names(popinfo)<-1:numsam

    # fill the list and the vector!
    for(s in 1:numsam){
        # get data from the row for this sample
        sampledata<-strsplit(rawdata[[3+numpop+s]],"\t")[[1]]
        # extract name and population
        popinfo[s]<-as.integer(sampledata[1])
        names(popinfo)[s]<-sampledata[2]
        dimnames(gendata)[[1]][s]<-sampledata[2]
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
            if(length(thesealleles)==1 && thesealleles==0){thesealleles<-missing}
            # get rid of "missing alleles" if the whole genotype is not missing
            thesealleles<-thesealleles[thesealleles !=0]
            # write the allele vector to the list
            gendata[[s,L]]<-thesealleles
        }
    }
    # return the vector and list
    return(list(PopData=popinfo,Genotypes=gendata))
}

write.GenoDive<-function(gendata, popnames="onebigpop",
                         commentline="file description goes here",
                         digits=2, file="", samples=dimnames(gendata)[[1]],
                         loci=dimnames(gendata)[[2]], popinfo=rep(1,times=length(samples)),
                         usatnts=rep(2,times=length(loci)), missing=-9){
    # fill in names if popinfo and usatnts are default
    if(is.null(names(popinfo))){names(popinfo)<-samples}
    if(is.null(names(usatnts))){names(usatnts)<-loci}
    # get some info for the second line of the file
    numsam<-length(samples)
    numloc<-length(loci)
    numpop<-length(popnames)
    if(numpop != length(unique(popinfo[samples]))){
        cat("Warning: number of populations in popnames and popinfo is different.\n")
    }
    maxploidy<-max(estimate.ploidy(gendata[samples,loci])[,"max.alleles"])

    # start a character vector to contain all the lines of the file
    lines<-c(commentline,paste(numsam,numpop,numloc,maxploidy,digits, sep="\t"))
    lines[3:(numpop+2)]<-popnames
    lines[numpop+3]<-paste("Population\tIndividual",paste(loci,sep="",collapse="\t"),sep="\t")

    # enter population and sample names
    for(s in 1:numsam){
        lines[numpop+3+s]<-paste(popinfo[samples[s]], samples[s], sep="\t")
    }

    # process alleles and write them to lines
    for(L in loci){
        # replace missing data with zeros
        repmiss<-function(value){
            if(value[1]==missing){ 0 } else {value}
        }
        convertedalleles<-mapply(repmiss, gendata[samples,L],SIMPLIFY=FALSE)
        # convert alleles to repeat number
        divallele<-function(value){floor(value/usatnts[L])}
        convertedalleles<-mapply(divallele,convertedalleles,SIMPLIFY=FALSE)
        # subtract if necessary to get them to the right number of digits
        suballele<-function(value){if(value[1]!=0){
            value-(10^(digits-1))} else{
            0}}
        while(max(mapply(max,convertedalleles)) >= 10^digits){
            convertedalleles<-mapply(suballele, convertedalleles, SIMPLIFY=FALSE)
        }
        # For each sample, concatenate into strings and add to lines
        for(s in 1:numsam){
            # Convert alleles to character and set up string for concatenation
            charalleles<-as.character(convertedalleles[[s]])
            allelestring<-""
            # For each allele:
            for(ca in charalleles){
                allele<-ca
                # Add zeros if necessary
                while(nchar(allele) < digits){
                    allele<-paste("0",allele,sep="")
                }
                # add this allele to the string
                allelestring<-paste(allelestring,allele,sep="")
            }
            # Add the allele string to the line for that sample.
            lines[numpop+3+s]<-paste(lines[numpop+3+s],allelestring,sep="\t")
        }
    }

    # write the file
    cat(lines, file=file, sep="\n")
}

write.Tetrasat<-function(gendata, commentline="insert data description here",
                         samples=dimnames(gendata)[[1]], loci=dimnames(gendata)[[2]],
                         popinfo=rep(1, length(samples)), usatnts=rep(2, length(loci)),
                         file="", missing=-9){
    # index popinfo and usatnts if not already done
    if(is.null(names(popinfo))){names(popinfo)<-samples}
    if(is.null(names(usatnts))){names(usatnts)<-loci}
    allpops<-unique(popinfo) # a vector of populations to cycle through
    # set up vector of lines to write to file
    lines<-c(commentline,loci)
    # where does population/genotype data begin?
    datastart<-length(lines) + 1

    # make lines with "Pop" and sample names, and an index of where samples are
    sampleindex<-as.integer(c())
    currentline<-datastart
    for(pop in allpops){
        lines[currentline]<-"Pop"
        currentline<-currentline+1
        # make a vector of all samples that belong to this population
        thesesamples<-samples[popinfo[samples]==pop]
        # make an index of the lines that they will go on
        indextoadd<-currentline:(currentline+length(thesesamples)-1)
        names(indextoadd)<-thesesamples
        # add this to the sample index
        sampleindex<-c(sampleindex,indextoadd)
        # write each sample name (plus spaces up to 20 characters) to the
        # appropriate line
        for(s in thesesamples){
            samname<-s
            while(nchar(samname)<20){
                samname<-paste(samname," ", sep="")
            }
            lines[currentline]<-samname
            currentline<-currentline+1
        }
    }

    # now go through the loci, convert the genotypes, and fill them in
    for(L in loci){
        # convert alleles to repeat number
        divallele<-function(value){
            if(value[1]==missing){
                value
            } else {
                floor(value/usatnts[L])
            }
        }
        convertedalleles<-mapply(divallele, gendata[samples,L],SIMPLIFY=FALSE)
        # make sure alleles have two digits or fewer
        suballele<-function(value){
            if(value[1]==missing){
                value
            } else {
                value-10
            }
        }
        while(max(mapply(max,convertedalleles)) >= 100){
            convertedalleles<-mapply(suballele, convertedalleles, SIMPLIFY=FALSE)
        }
        # go through the genotypes by sample
        for(s in samples){
            # convert missing data to blank spaces, convert to alleles to characters
            if(convertedalleles[[s]][1]==missing){
                charalleles<-"         "
            } else {
                charalleles<-as.character(convertedalleles[[s]])
                # duplicate allele if fully homozygous
                if(length(charalleles)==1){
                    charalleles<-rep(charalleles,4)
                }
                # randomly remove alleles if there are more than 4
                if(length(charalleles)>4){
                    charalleles<-sample(charalleles,4,replace=FALSE)
                    cat(c("Alleles randomly removed:",s,L,"\n"),sep=" ")
                }
            }
            # concatenate all alleles into one string
            allelestring<-""
            for(ca in charalleles){
                allele<-ca
                # Add zeros if necessary
                while(nchar(allele) < 2){
                    allele<-paste("0",allele,sep="")
                }
                # add this allele to the string
                allelestring<-paste(allelestring,allele,sep="")
            }
            # add spaces up to nine characters
            while(nchar(allelestring) < 9){
                allelestring<-paste(allelestring," ", sep="")
            }
            # add the concatenated genotype to the appropriate line
            lines[sampleindex[s]]<-paste(lines[sampleindex[s]],allelestring,sep="")
        }
    }

    cat(lines, sep="\n", file=file)
}

write.ATetra<-function(gendata, samples=dimnames(gendata)[[1]], loci=dimnames(gendata)[[2]],
                       popinfo=rep(1,length(samples)), popnames="onebigpop",
                       commentline="insert data info here", missing=-9, file=""){
    # name popinfo if necessary
    if(is.null(names(popinfo))){names(popinfo)<-samples}
    # set up a character vector to hold the lines for the file
    lines<-paste("TIT",commentline, sep=",")
    currentline<-2 # a variable to say which line to write to next
    # make numbers to go with loci, pops, and samples
    locnums<-1:length(loci)
    names(locnums)<-loci
    samnums<-1:length(samples)
    names(samnums)<-samples
    popnums<-1:length(popnames)
    names(popnums)<-popnames

    if(length(popnames) != length(unique(popinfo))){
      cat("Warning: number of populations in popnames and popinfo is different.\n")
    }

    # fill in data using a loop
    for(L in loci){
        # write a line describing the locus
        lines[currentline]<-paste("LOC",locnums[L],L,sep=",")
        currentline<-currentline+1
        for(p in popnames){
            # write a line describing the population
            lines[currentline]<-paste("POP",locnums[L],popnums[p],p, sep=",")
            currentline<-currentline+1
            # get a vector of individuals in this population
            thesesamples<-names(popinfo)[popinfo==popnums[p]]
            thesesamples<-thesesamples[thesesamples %in% samples]
            # write lines for individuals
            for(s in thesesamples){
                # first put sample info into the line
                lines[currentline]<-paste("IND",locnums[L],popnums[p],samnums[s],s,sep=",")
                # get the alleles and print a warning if there is missing data
                thesealleles<-gendata[[s,L]]
                if(thesealleles[1]==missing){
                    thesealleles<-""
                    cat("Missing data:",s,L,"\n",sep=" ")
                }
                # take a random sample if there are more than 4
                if(length(thesealleles)>4){
                    thesealleles<-sample(thesealleles,4,replace=FALSE)
                    cat("More than 4 alleles:",s,L,"\n",sep=" ")
                }
                # add alleles to the line
                for(a in 1:4){
                    if(length(thesealleles)>=a ){
                      lines[currentline]<-paste(lines[currentline],thesealleles[a],sep=",")
                  } else {
                      lines[currentline]<-paste(lines[currentline],"",sep=",")
                  }
                }
                currentline<-currentline+1
            }
        }
    }

    # write the "END-record"
    lines[currentline]<-"END"

    # write the file
    cat(lines, sep="\n", file=file)
}

read.Structure<-function(infile,missingin=-9,missingout=-9,sep="\t",markernames=TRUE,
                         labels=TRUE, extrarows=1, extracols=0, ploidy=4,
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

    # set up the list to store the genotypes
    gendata<-array(list(missingout),dim=c(length(samples),length(loci)),
                   dimnames=list(samples,loci))

    # fill the list
    for(L in loci){
        for(s in samples){
            rawalleles<-rawdata[rawdata$Samples==s,L]
            # process missing data and get unique alleles
            thesealleles<-unique(rawalleles)
            if(length(thesealleles)==1 && thesealleles[1]==missingin){
                thesealleles<-missingout
            } else {
                thesealleles<-thesealleles[thesealleles != missingin]
            }
            gendata[[s,L]]<-thesealleles
        }
    }

    # extract the extra columns, if needed
    if(getexcols){
        Extracol<-data.frame(row.names=samples)
        if(labels){
            colsbeforex<-1
        } else {
            colsbeforex<-0
        }
        for(x in extracols){
            # look up values by sample in this column
            # add the column to the data frame
            Extracol[[x]]<-rawdata[[x+colsbeforex]][seq(1+extrarows,
                                                        length(rawdata[[1]])+1-ploidy,
                                                        by=ploidy)]
        }
    }

    # return extra columns and genotypes
    if(getexcols){
        return(list(ExtraCol=Extracol, Genotypes=gendata))
    } else {
        return(gendata)
    }
}

write.GeneMapper<-function(gendata,file="",samples=dimnames(gendata)[[1]],
                           loci=dimnames(gendata)[[2]]){
    # figure out how many allele columns are needed
    numallelecol<-max(estimate.ploidy(gendata[samples,loci])[,"max.alleles"])
    # figure out how many rows are needed
    numrows<-length(samples)*length(loci)
    # set up data frame
    gentable<-data.frame(Sample.Name=rep(samples, times=length(loci)),
                         Marker=rep(loci, each=length(samples)))
    # put empty allele columns into data frame
    alcollabels<-paste("Allele.",1:numallelecol,sep="")
    for(ac in 1:numallelecol){
        gentable[[ac+2]]<-rep("",times=numrows)
        names(gentable)[ac+2]<-alcollabels[ac]
    }

    # put alleles into their respective cells
    currentrow<-1
    for(L in loci){
        for(s in samples){
            thesealleles<-as.character(gendata[[s,L]])
            gentable[currentrow,3:(length(thesealleles)+2)]<-thesealleles
            currentrow<-currentrow+1
        }
    }

    # write the table to file
    write.table(gentable, file=file, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

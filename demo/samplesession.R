# load the data from the package into R
data(testgenotypes)
data(FCRinfo)

# commands from "How genotypes are stored in polysat"
testgenotypes
testgenotypes[,"RhCBA23"]
testgenotypes["FCR14",]
myloci <- c("RhCBA23","RhCBA28")
mysamples <- c("FCR1","FCR2","FCR3","FCR4","FCR5",
               "FCR6","FCR7","FCR8","FCR9","FCR10")
subgenotypes <- testgenotypes[mysamples, myloci]
subgenotypes
testdist <- meandistance.matrix(testgenotypes[mysamples, myloci])
testdist <- meandistance.matrix(testgenotypes, samples = mysamples,
                                 loci = myloci)
to.exclude <- c("FCR11","FCR12")
all.samples <- dimnames(testgenotypes)[[1]]
to.analyze <- all.samples[!all.samples %in% to.exclude]
to.analyze
testgenotypes[to.analyze,]
find.missing.gen(testgenotypes)

# "Examples of how to edit genotype data in R"
dimnames(testgenotypes)[[2]] <- c("C15","C23","C28")
testgenotypes[["FCR5","C15"]] <- 208
testgenotypes[["FCR19","C23"]]
testgenotypes[["FCR19","C23"]] <- c(98, 125)
testgenotypes[["FCR2","C23"]]
testgenotypes[["FCR2","C23"]][4] <- 112
testgenotypes["FCR5",]
testgenotypes["FCR19",]
testgenotypes["FCR2",]
testgenotypes["FCR7","C28"]
testgenotypes[["FCR7","C28"]]
#testgenotypes[["FCR7",]]

# "Editing data using spreadsheet software"
write.GeneMapper(testgenotypes, "genotypestoedit.txt")
testgenotypes <- read.GeneMapper("genotypestoedit.txt")

# "Merging genotype objects"
mygenotypes1 <- array(list(1,2,3,4), dim=c(2,2),
                      dimnames=list(c("ind1","ind2"),c("loc1","loc2")))
mygenotypes2 <- array(list(5,6,7,8), dim=c(2,2),
                      dimnames=list(c("ind1","ind2"),c("loc3","loc4")))
mygenotypes1
mygenotypes2
mysamples <- c("ind1","ind2")
mygenotypes <- cbind(mygenotypes1[mysamples,], mygenotypes2[mysamples,])
mygenotypes
mygenotypes2 <- array(list(9,10,11,12), dim=c(2,2),
                      dimnames=list(c("ind3","ind4"),c("loc1","loc2")))
mygenotypes2
myloci <- c("loc1","loc2")
mygenotypes <- rbind(mygenotypes1[,myloci], mygenotypes2[,myloci])
mygenotypes

write.GeneMapper(mygenotypes1, "mygenotypes1.txt")
write.GeneMapper(mygenotypes2, "mygenotypes2.txt")
mygenotypes <- read.GeneMapper(c("mygenotypes1.txt", "mygenotypes2.txt"))


# "Creating a genotype object from scratch"
missing <- -9
samples <- c("ind1","ind2","ind3")
loci <- c("loc1","loc2")
gendata <- array(list(missing), dim = c(length(samples), length(loci)),
                 dimnames = list(samples, loci))
gendata
gendata[["ind1","loc1"]] <- c(100,102,104)
#for(L in loci){
#	for(s in samples){
#
#		# Insert code here that would find the genotype of sample
#		# s at locus L in your data structure and convert it to a
#		# vector called thesealleles.
#
#		gendata[[s,L]] <- thesealleles
#	}
#}

# "Importing data from files"
#folderpath <- "C:\\Users\\lvclark\\R\\win-library\\2.11\\polysat\\extdata\\"
#ATdata <- read.ATetra(paste(folderpath,"ATetraExample.txt",sep=""))
#ATdata
#Tetdata <- read.Tetrasat(paste(folderpath,"tetrasatExample.txt",sep=""))
#Tetdata
#GDdata <- read.GenoDive(paste(folderpath,"genodiveExample.txt",sep=""))
#GDdata
#GMdata <- read.GeneMapper(paste(folderpath,"GeneMapperCBA",
#                                c("15.txt","23.txt","28.txt"), sep=""))
#GMdata
#Structdata <- read.Structure(paste(folderpath, "structureExample.txt", sep=""),
#                             extracols = 1, ploidy = 8, getexcols = TRUE)
#Structdata
#Spagdata <- read.SPAGeDi(paste(folderpath, "spagediExample.txt", sep=""),
#                         returnploidies = TRUE)
#Spagdata
#Spagdata$Genotypes[,"locA"]
#Domdata <- as.matrix(read.table(paste(folderpath, "dominantExample.txt", sep=""),
#                               header=TRUE,sep="\t",row.names=1))
#Domdata
#ConvDomdata <- dominant.to.codominant(Domdata)
#ConvDomdata

# "Exporting genotype data to files"
mypopinfo <- FCRinfo$Species
names(mypopinfo) <- row.names(FCRinfo)
mypopinfo
mypopnames <- c("A","B","C")

write.GeneMapper(testgenotypes, file="GMout.txt")
write.Structure(testgenotypes, ploidy=8, file="Structout.txt",
                indploidies=c(8,8,8,4,8,8,rep(4,14)),
                extracols=array(mypopinfo, dim=c(20,1),
                                dimnames=list(names(mypopinfo),"PopData")))
write.GenoDive(testgenotypes, mypopnames, file="GDout.txt",popinfo=mypopinfo)
tetrasamples <- names(mypopinfo)[mypopinfo != 1]
tetrasamples
write.ATetra(testgenotypes, popinfo=mypopinfo,
             popnames=mypopnames,
             file="ATout.txt", samples=tetrasamples)
write.Tetrasat(testgenotypes, popinfo=mypopinfo[tetrasamples],
               file="TSout.txt", samples=tetrasamples)
write.SPAGeDi(testgenotypes, file="SpagOut.txt", popinfo=mypopinfo,
              spatcoord=data.frame(Lat=c(rep(43.943,6), rep(43.957,14)),
              Long=c(rep(-122.768,6), rep(-122.755,14))),
              indploidies=c(8,8,8,4,8,8,rep(4,14)))
Domdata <- codominant.to.dominant(testgenotypes)
write.table(Domdata, file="Domout.txt")

# "Individual level statistics"
# calculate distance matrices from the genotypes,
# using two different distance metrics
Bmatrix <- meandistance.matrix(testgenotypes, progress=FALSE)
Lmatrix <- meandistance.matrix(testgenotypes, distmetric=Lynch.distance,
                               progress=FALSE)
# do Principal Coordinate Analysis
Bprcoord <- cmdscale(Bmatrix)
Lprcoord <- cmdscale(Lmatrix)
plot(Bprcoord[,1],Bprcoord[,2], col=FCRinfo$Plot.color,
     pch=FCRinfo$Plot.symbol)
plot(Lprcoord[,1],Lprcoord[,2], col=FCRinfo$Plot.color,
     pch=FCRinfo$Plot.symbol)
# histogram
hist(as.vector(Bmatrix))
hist(as.vector(Lmatrix))
# export
write.table(Bmatrix, file="Bmatrix.txt")
# experiment with removing loci
Larray <- meandistance.matrix(testgenotypes, progress=FALSE,
                              distmetric=Lynch.distance, all.distances=TRUE)[[1]]
mdist15.23 <- meandist.from.array(Larray, loci=c("C15","C23"))
mdist23.28 <- meandist.from.array(Larray, loci=c("C23","C28"))
mdist15.28 <- meandist.from.array(Larray, loci=c("C15","C28"))

# "Estimate ploidy of samples"
myploidies <- as.data.frame(estimate.ploidy(testgenotypes))
myploidies
myploidies[[3]]<-myploidies$max.alleles
names(myploidies)[3] <- "ploidy"
myploidies <- edit(myploidies)

# "Population-level statistics"
# "Estimate allele frequencies"
mypopinfo <- FCRinfo$Species
names(mypopinfo) <- row.names(FCRinfo)
myploidies <- c(8,8,8,4,8,8,rep(4,14))
names(myploidies) <- row.names(FCRinfo)
freqtable <- estimate.freq(testgenotypes, popinfo=mypopinfo,
                           indploidies=myploidies)
freqtable[,1:10]

# "Calculating pairwise FST"
testfsts<-calcFst(freqtable)
testfsts

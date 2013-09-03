### R code from vignette source 'polysattutorial.Rnw'

###################################################
### code chunk number 1: polysattutorial.Rnw:25-27
###################################################
library(polysat)
options(width=60)


###################################################
### code chunk number 2: polysattutorial.Rnw:209-210
###################################################
getwd()


###################################################
### code chunk number 3: polysattutorial.Rnw:217-218
###################################################
simgen <- read.GeneMapper("GeneMapperExample.txt")


###################################################
### code chunk number 4: polysattutorial.Rnw:226-231
###################################################
summary(simgen)
Samples(simgen)
Loci(simgen)
viewGenotypes(simgen, samples=paste("A", 1:20, sep=""), loci="loc1")
find.missing.gen(simgen)


###################################################
### code chunk number 5: polysattutorial.Rnw:239-243
###################################################
Description(simgen) <- "Dataset for the tutorial"
PopNames(simgen) <- c("PopA", "PopB", "PopC")
PopInfo(simgen) <- rep(1:3, each = 100)
Usatnts(simgen) <- c(2, 3, 2)


###################################################
### code chunk number 6: polysattutorial.Rnw:250-252
###################################################
rep(1:3, each = 100)
PopInfo(simgen)


###################################################
### code chunk number 7: polysattutorial.Rnw:257-258
###################################################
Samples(simgen, populations = "PopA")


###################################################
### code chunk number 8: polysattutorial.Rnw:269-270
###################################################
Usatnts(simgen)


###################################################
### code chunk number 9: polysattutorial.Rnw:305-307
###################################################
simgen <- reformatPloidies(simgen, output="sample")
Ploidies(simgen) <- read.table("vignettebuild/simgenPloidies2.txt")[[1]]


###################################################
### code chunk number 10: polysattutorial.Rnw:312-313
###################################################
summary(simgen)


###################################################
### code chunk number 11: polysattutorial.Rnw:324-325
###################################################
save(simgen, file="simgen.RData")


###################################################
### code chunk number 12: polysattutorial.Rnw:344-345
###################################################
load("vignettebuild/testmat.RData")


###################################################
### code chunk number 13: polysattutorial.Rnw:348-351
###################################################
pca <- cmdscale(testmat)
plot(pca[,1], pca[,2], col=rep(c("red", "green", "blue"), each=100),
     main = "PCA with Bruvo distance")


###################################################
### code chunk number 14: polysattutorial.Rnw:362-363
###################################################
load("vignettebuild/testmat2.RData")


###################################################
### code chunk number 15: polysattutorial.Rnw:366-369
###################################################
pca2 <- cmdscale(testmat2)
plot(pca2[,1], pca2[,2], col=rep(c("red", "green", "blue"), each=100),
     main = "PCA with Lynch distance")


###################################################
### code chunk number 16: polysattutorial.Rnw:390-393
###################################################
simgen2 <- deleteSamples(simgen, c("B59", "C30"))
simgen2 <- deleteLoci(simgen2, "loc2")
summary(simgen2)


###################################################
### code chunk number 17: polysattutorial.Rnw:401-405
###################################################
samToUse <- Samples(simgen2, populations=c("PopA", "PopB"), ploidies=4)
exclude <- c("A50", "A78", "B25", "B60", "B81")
samToUse <- samToUse[!samToUse %in% exclude]
samToUse


###################################################
### code chunk number 18: polysattutorial.Rnw:415-416
###################################################
summary(simgen2[samToUse, "loc1"])


###################################################
### code chunk number 19: polysattutorial.Rnw:423-428
###################################################
testmat3 <- meandistance.matrix(simgen2, samples = samToUse,
                                distmetric = Lynch.distance,
                                progress= FALSE)
pca3 <- cmdscale(testmat3)
plot(pca3[,1], pca3[,2], col=c("red", "blue")[PopInfo(simgen2)[samToUse]])


###################################################
### code chunk number 20: polysattutorial.Rnw:475-476
###################################################
simfreq <- read.table("vignettebuild/simfreq.txt")


###################################################
### code chunk number 21: polysattutorial.Rnw:479-484
###################################################
simfreq
simFst <- calcFst(simfreq)
simFst
simFst12 <- calcFst(simfreq, loci=c("loc1", "loc2"))
simFst12


###################################################
### code chunk number 22: polysattutorial.Rnw:500-501
###################################################
write.Structure(simgen, ploidy = 4, file="simgenStruct.txt")


###################################################
### code chunk number 23: polysattutorial.Rnw:551-552
###################################################
showClass("genambig")


###################################################
### code chunk number 24: polysattutorial.Rnw:560-563
###################################################
mysamples <- c("indA", "indB", "indC", "indD", "indE", "indF")
myloci <- c("loc1", "loc2", "loc3")
mydataset <- new("genambig", samples=mysamples, loci=myloci)


###################################################
### code chunk number 25: polysattutorial.Rnw:569-570
###################################################
mydataset


###################################################
### code chunk number 26: polysattutorial.Rnw:588-603
###################################################
Loci(mydataset)
Loci(mydataset) <- c("L1", "L2", "L3")
Loci(mydataset)
Samples(mydataset)
Samples(mydataset)[3] <- "indC1"
Samples(mydataset)
PopNames(mydataset) <- c("Yosemite", "Sequoia")
PopInfo(mydataset) <- c(1,1,1,2,2,2)
PopInfo(mydataset)
PopNum(mydataset, "Yosemite")
PopNum(mydataset, "Sequoia") <- 3
PopNames(mydataset)
PopInfo(mydataset)
Ploidies(mydataset) <- c(4,4,4,4,4,6)
Ploidies(mydataset)


###################################################
### code chunk number 27: polysattutorial.Rnw:606-621
###################################################
Ploidies(mydataset)["indC1",] <- 6
Ploidies(mydataset)
Usatnts(mydataset) <- c(2,2,2)
Usatnts(mydataset)
Description(mydataset) <- "Tutorial, part 2."
Description(mydataset)
Genotypes(mydataset, loci="L1") <- list(c(122, 124, 128), c(124,126),
                     c(120,126,128,130), c(122,124,130), c(128,130,132),
                     c(126,130))
Genotype(mydataset, "indB", "L3") <- c(150, 154, 160)
Genotypes(mydataset)
Genotype(mydataset, "indD", "L1")
Missing(mydataset)
Missing(mydataset) <- -1
Genotypes(mydataset)


###################################################
### code chunk number 28: polysattutorial.Rnw:627-629
###################################################
mydataset@Genotypes
mydataset@Genotypes[["indB","L1"]]


###################################################
### code chunk number 29: polysattutorial.Rnw:642-645
###################################################
isMissing(mydataset, "indA", "L2")
isMissing(mydataset, "indA", "L1")
isMissing(mydataset)


###################################################
### code chunk number 30: polysattutorial.Rnw:652-662
###################################################
moredata <- new("genambig", samples=c("indG", "indH"), loci=Loci(mydataset))
Usatnts(moredata) <- Usatnts(mydataset)
Description(moredata) <- Description(mydataset)
PopNames(moredata) <- "Kings Canyon"
PopInfo(moredata) <- c(1,1)
Ploidies(moredata) <- c(4,4)
Missing(moredata) <- Missing(mydataset)
Genotypes(moredata, loci="L1") <- list(c(126,130,136,138), c(124,126,128))
mydataset2 <- merge(mydataset, moredata)
mydataset2


###################################################
### code chunk number 31: polysattutorial.Rnw:713-726
###################################################
ploidyexample <- new("genambig")
Samples(ploidyexample)
Loci(ploidyexample)
Ploidies(ploidyexample)
ploidyexample <- reformatPloidies(ploidyexample, output="locus")
Ploidies(ploidyexample)
ploidyexample <- reformatPloidies(ploidyexample, output="sample")
Ploidies(ploidyexample)
ploidyexample <- reformatPloidies(ploidyexample, output="one")
Ploidies(ploidyexample)
Ploidies(ploidyexample) <- 4
ploidyexample <- reformatPloidies(ploidyexample, output="matrix")
Ploidies(ploidyexample)


###################################################
### code chunk number 32: polysattutorial.Rnw:734-735
###################################################
Ploidies(ploidyexample)["ind1", "loc1"]


###################################################
### code chunk number 33: polysattutorial.Rnw:743-746
###################################################
Ploidies(ploidyexample, "ind1", "loc1")
ploidyexample <- reformatPloidies(ploidyexample, output="one")
Ploidies(ploidyexample, "ind1", "loc1")


###################################################
### code chunk number 34: polysattutorial.Rnw:775-778
###################################################
simgenB <- genambig.to.genbinary(simgen)
Genotypes(simgenB, samples=paste("A", 1:20, sep=""), loci="loc1")
PopInfo(simgenB)[Samples(simgenB, ploidies=2)]


###################################################
### code chunk number 35: polysattutorial.Rnw:784-785
###################################################
write.table(Genotypes(simgenB), file="simBinaryData.txt")


###################################################
### code chunk number 36: polysattutorial.Rnw:835-839
###################################################
GDdata <- read.GenoDive("genodiveExample.txt")
Structdata <- read.Structure("structureExample.txt", ploidy = 8)
Spagdata <- read.SPAGeDi("spagediExample.txt")
PDdata <- read.POPDIST(c("POPDISTexample1.txt", "POPDISTexample2.txt"))


###################################################
### code chunk number 37: polysattutorial.Rnw:878-881
###################################################
GMdata <- read.GeneMapper(c("GeneMapperCBA15.txt",
                            "GeneMapperCBA23.txt",
                            "GeneMapperCBA28.txt"))


###################################################
### code chunk number 38: polysattutorial.Rnw:890-896
###################################################
# view the format
read.table("STRandExample.txt", sep="\t", header=TRUE)
# import the data
STRdata <- read.STRand("STRandExample.txt")
Samples(STRdata)
PopNames(STRdata)


###################################################
### code chunk number 39: polysattutorial.Rnw:905-907
###################################################
domdata <- read.table("dominantExample.txt", header=TRUE,
                      sep="\t", row.names=1)


###################################################
### code chunk number 40: polysattutorial.Rnw:915-920
###################################################
domdata
domdata <- as.matrix(domdata)
PAdata <- new("genbinary", samples=c("ind1", "ind2", "ind3"),
              loci=c("ABC1", "ABC2"))
Genotypes(PAdata) <- domdata


###################################################
### code chunk number 41: polysattutorial.Rnw:928-930
###################################################
PopInfo(PAdata) <- c(1,1,2)
PAdata <- genbinary.to.genambig(PAdata)


###################################################
### code chunk number 42: polysattutorial.Rnw:959-965
###################################################
myexcol <- array(c(rep(0:1, each=150), seq(0.1, 30, by=0.1)), dim=c(300,2),
                 dimnames = list(Samples(simgen), c("PopFlag", "Something")))
myexcol[1:10,]
write.Structure(simgen, ploidy=4, file="simgenStruct2.txt",
                writepopinfo = FALSE, extracols = myexcol,
                missingout = -1)


###################################################
### code chunk number 43: polysattutorial.Rnw:973-974
###################################################
write.GenoDive(simgen, file="simgenGD.txt")


###################################################
### code chunk number 44: polysattutorial.Rnw:986-987
###################################################
write.SPAGeDi(simgen, file="simgenSpag.txt")


###################################################
### code chunk number 45: polysattutorial.Rnw:1000-1002
###################################################
write.POPDIST(simgen, samples = Samples(simgen, ploidies=4),
              file = "simgenPOPDIST.txt")


###################################################
### code chunk number 46: polysattutorial.Rnw:1012-1013
###################################################
write.GeneMapper(simgen, file="simgenGM.txt")


###################################################
### code chunk number 47: polysattutorial.Rnw:1022-1025
###################################################
simgenPA <- genambig.to.genbinary(simgen)
write.table(Genotypes(simgenPA), file="simgenPA.txt", quote=FALSE,
            sep = ",")


###################################################
### code chunk number 48: polysattutorial.Rnw:1037-1039
###################################################
write.table(data.frame(Ploidies(simgen), row.names=Samples(simgen)),
            file="simgenPloidies.txt")


###################################################
### code chunk number 49: polysattutorial.Rnw:1090-1091
###################################################
load("vignettebuild/testmat4.RData")


###################################################
### code chunk number 50: polysattutorial.Rnw:1093-1096
###################################################
pca4 <- cmdscale(testmat4)
plot(pca4[,1], pca4[,2], col=c("red", "blue")[PopInfo(simgen)[samToUse]],
     main="Bruvo distance with meandistance.matrix2")


###################################################
### code chunk number 51: polysattutorial.Rnw:1104-1105
###################################################
hist(as.vector(testmat))


###################################################
### code chunk number 52: polysattutorial.Rnw:1107-1108
###################################################
hist(as.vector(testmat2))


###################################################
### code chunk number 53: polysattutorial.Rnw:1110-1111
###################################################
write.table(testmat2, file="simgenDistMat.txt")


###################################################
### code chunk number 54: polysattutorial.Rnw:1132-1142
###################################################
subsamples <- Samples(simgen, populations=1)
subsamples <- subsamples[!isMissing(simgen, subsamples, "loc1") &
                         !isMissing(simgen, subsamples, "loc2") &
                         !isMissing(simgen, subsamples, "loc3")]
Larray <- meandistance.matrix(simgen, samples=subsamples,
                              progress=FALSE,
                 distmetric=Lynch.distance, all.distances=TRUE)[[1]]
mdist1.2 <- meandist.from.array(Larray, loci=c("loc1","loc2"))
mdist2.3 <- meandist.from.array(Larray, loci=c("loc2","loc3"))
mdist1.3 <- meandist.from.array(Larray, loci=c("loc1","loc3"))


###################################################
### code chunk number 55: polysattutorial.Rnw:1165-1168
###################################################
clones <- assignClones(testmat, samples=paste("A", 1:100, sep=""),
                       threshold=0.2)
clones


###################################################
### code chunk number 56: polysattutorial.Rnw:1185-1188
###################################################
simal <- alleleDiversity(simgen)
simal$counts
simal$alleles[["PopA","loc1"]]


###################################################
### code chunk number 57: polysattutorial.Rnw:1227-1231
###################################################
simFst
simfreqSimple <- simpleFreq(simgen, samples = Samples(simgen, ploidies=4))
simFstSimple <- calcFst(simfreqSimple)
simFstSimple


###################################################
### code chunk number 58: polysattutorial.Rnw:1244-1245
###################################################
write.freq.SPAGeDi(simfreq, usatnts=Usatnts(simgen), file="SPAGfreq.txt")


###################################################
### code chunk number 59: polysattutorial.Rnw:1255-1256
###################################################
gpsimfreq <- freq.to.genpop(simfreq)


###################################################
### code chunk number 60: polysattutorial.Rnw:1286-1287
###################################################
load("vignettebuild/testmat5.RData")


###################################################
### code chunk number 61: polysattutorial.Rnw:1289-1293
###################################################
simdiv <- genotypeDiversity(simgen, d=testmat5, threshold=0.2, index=Shannon)
simdiv
simdiv2 <- genotypeDiversity(simgen, d=testmat5, threshold=0.2, index=Simpson)
simdiv2


###################################################
### code chunk number 62: polysattutorial.Rnw:1299-1303
###################################################
simdiv2var <- genotypeDiversity(simgen, d=testmat5, threshold=0.2,
                                index=Simpson.var)
simdiv2 - 2*simdiv2var^0.5
simdiv2 + 2*simdiv2var^0.5


###################################################
### code chunk number 63: polysattutorial.Rnw:1345-1347
###################################################
ATdata <- read.ATetra("ATetraExample.txt")
Tetdata <- read.Tetrasat("tetrasatExample.txt")


###################################################
### code chunk number 64: polysattutorial.Rnw:1358-1361
###################################################
write.ATetra(simgen, samples=Samples(simgen, ploidies=4), file="simgenAT.txt")
write.Tetrasat(simgen, samples=Samples(simgen, ploidies=4),
               file="simgenTet.txt")


###################################################
### code chunk number 65: polysattutorial.Rnw:1436-1440
###################################################
Present(simgenB) <- "P"
Absent(simgenB) <- 2
Missing(simgenB) <- 0
Genotypes(simgenB)[1:10, 1:6]


###################################################
### code chunk number 66: polysattutorial.Rnw:1446-1449
###################################################
genmat <- Genotypes(simgenB)
dimnames(genmat)[[2]] <- paste("M", 1:dim(genmat)[2], sep="")
genmat[1:10, 1:10]



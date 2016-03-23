## http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
## http://www.homolog.us/blogs/blog/2013/11/22/cshl-keynote-talk-lior-pachter/
## https://www.biostars.org/p/133488/


##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
pathToCountFiles=args[1]
target=args[2]
inputreads=as.numeric(args[3])
#print(pathToCountFiles)
print(target)

library("edgeR")

files <- dir(path=pathToCountFiles,pattern=paste("^", target, ".*quant.counts$", sep = ""))
print(files)
data <- readDGE(files)
head(data$counts)
rawCounts=rowSums(data$counts)  ## raw counts mapping to transcripts

data2=read.table("transcripts.lengthes", header=T)
length=data2$length

dataSummary=cbind(length,rawCounts)
#head(dataSummary)

normCount=dataSummary[,2]/(dataSummary[,1]/1000)  ## normalize the count by read length
normCountSum=sum(normCount)
calcTPM=(dataSummary[,2] * 10^6) / ((dataSummary[,1]/1000)*normCountSum)
calcTPMinputs=(dataSummary[,2] * 10^6) / ((dataSummary[,1]/1000)*inputreads)
dataSummary2=cbind(dataSummary,calcTPM,normCount,calcTPMinputs)
#head(dataSummary2)

#data3=read.table("gene_transcript.map")
#stopifnot(all(rownames(dataSummary2)==data3$V1))

write.table(dataSummary2, file=paste(target, 'dataSummary_comp', sep = "."), sep='\t', quote=F, row.names=T, col.names=NA)


################################
#pdf("rawCount-hist-plot.pdf")
#hist(rawCounts,xlim=c(0, 100), breaks = 10000000)
#dev.off()

#pdf("TPM-hist-plot.pdf")
#hist(calcTPM,xlim=c(0, 20), breaks = 100000)
#dev.off()

#pdf("isoformTPM-hist-plot.pdf")
#hist(as.numeric(isoformTPM),xlim=c(0, 100), breaks = 100)
#dev.off()

#pdf("isoformTPM-hist-plot-focus.pdf")
#hist(as.numeric(isoformTPM),xlim=c(0, 10), breaks = 100)
#dev.off()

#pdf("isoformTPM-hist-plot-focus2.pdf")
#hist(as.numeric(isoformTPM),xlim=c(0, 5), breaks = 1000)
#dev.off()

#reducedData=extCounts[totCounts>20,]
#write.table(reducedData, file='reducedData', sep='\t', quote=F, row.names=T, col.names=T)


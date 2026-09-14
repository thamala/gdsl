#Script for processing SNPs called with MUMmer.
#Tuomas Hämälä 20206

library(GenomicRanges)
library(fields)
library(reshape2)
library(ggplot2)

coords <- list.files(pattern="*.coords") #Coordinate files from each comparison (show-coords)
snps <- list.files(pattern="*.snps") #SNP from each comparison (show-snps)
target <- makeGRangesFromDataFrame(data.frame(Chr="chr3A",Start=(407654122-5e+5),End=(407655694+5e+5))) #Focal region

#Keep only regions aligned in each comparison
n <- length(coords)
for(i in 1:n){
	if(i==1|i==n|i%%2==0)  cat("File",i,"/",n,"\n")
	d <- read.table(coords[i])
	d <- makeGRangesFromDataFrame(data.frame(Chr=d$V8,Start=d$V1,End=d$V2))
	d <- subsetByOverlaps(d, target)
	if(i == 1) keep <- d
	else keep <- subsetByOverlaps(keep, d)
}

#Compile SNPs and defined REF alleles
n <- length(snps)
for(i in 1:n){
	if(i==1|i==n|i%%2==0)  cat("File",i,"/",n,"\n")
	d <- read.table(snps[i])
	d <- makeGRangesFromDataFrame(data.frame(Chr=d$V9,Start=d$V1,End=d$V1))
	temp <- as.data.frame(subsetByOverlaps(d, keep))
	temp <- data.frame(ps=temp$start,x=1)
	names(temp)[2] <- substr(snps[i],1,nchar(snps[i])-5)
	if(i == 1) out <- temp
	else out <- merge(out, temp, by="ps", all=T)
}
out[is.na(out)] <- 0
out <- out[order(out$ps),]

#Calculate allele frequencies in 50 bins, corresponding to 20 kb 
for(i in 2:ncol(out)){
	b <- stats.bin(out$ps, out[,i], N=50, prettyBins=T)
	d <- data.frame(ps=b$centers,x=b$stats[2,])
	names(d)[2] <- names(out)[i]
	if(i == 2) o <- d
	else o <- merge(o,d,by="ps")
}
out <- na.omit(o)

#PCA. Y-axis will be sorted according to PC1
pca <- prcomp(t(out[,-1]))
pc1 <- pca$x[,1]
pc1 <- pc1[order(pc1)]
out$ps <- as.factor(out$ps)
m <- melt(out)
ord <- factor(m$variable, levels=names(pc1))
ord_ps <- factor(m$ps, levels=out$ps)

#Plot SNP matrix
ggplot(m, aes(x=ord_ps,y=ord,fill=value,color=value))+
	geom_tile(color=NA)+
	labs(x="Chromosome 3A (407–408 mb)",y="Accession",fill="Genotype",title="3AG01900850")+
	scale_fill_gradient2(low="white",mid="#fee5d9",high="#de2d26",midpoint=0.3,limits=c(0,0.8),breaks=c(0,0.2,0.4,0.6,0.8),labels=c(0,0.2,0.4,0.6,0.8))+
	scale_color_gradient2(low="white",mid="#fee5d9",high="#de2d26",midpoint=0.3)+
	guides(color = "none")+
	theme(panel.background = element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		panel.border=element_rect(colour="black",fill=NA,linewidth=0.5),
		axis.line=element_blank(),
		axis.text=element_blank(),
		axis.ticks.length=unit(.15, "cm"),
		axis.ticks=element_blank(),
		plot.title=element_text(size=12, color="black"),
		axis.title.x=element_text(size=11, color="black"),
		axis.title.y=element_blank(),
		legend.title=element_text(size=12, color="black"),
		legend.text=element_text(size=11, color="black"),
		legend.key=element_blank())

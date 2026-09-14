#Scripts for analysing OrthoFinder results to define PAVs and CNVs in oat
#Tuomas Hämälä, April 2026

library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(viridis)

df <- read.table("N0.tsv",header=T,sep="\t",na.strings="") #Hierarchical orthogroups from OrthoFinder. Oat genomes were split by subgenome and named as "id.subgenome" (e.g., "ASLAK.A"). Note, file was modified to remove spaces after commas before loading into R (i.e., ", " replaced with ",").
all_genes <- read.table("OF_prot_id.txt") #File listing all genes used in OrthoFinder. Format: genome id, gene id.
df <- df[,-c(2:3)]

#Compile hog data
hog <- data.frame(hog=integer(),id=integer(),subgenome=integer(),genes=I(list()))
nrow <- nrow(df)
ncol <- ncol(df)
k <- 1
for(i in 1:nrow){
	for(j in 1:ncol){
		if(j == 1){
			for(l in 1:(ncol-1)) hog[k+l-1,1] <- df[i,1]
		}
		else{
			hog[k,2] <- unlist(strsplit(names(df)[j],".",fixed=T))[1]
			hog[k,3] <- unlist(strsplit(names(df)[j],".",fixed=T))[2]
			hog[[k,4]] <- unlist(strsplit(df[i,j],","))
			k <- k + 1
		}
	}
}


#Define PAVs and CNVs
pav <- data.frame(id=unique(hog$id),core=0,shell=0,cloud=0)
cnv <- data.frame(hog=integer(),id=integer(),n=integer())
h <- unique(hog$hog)
n <- length(h)
inds <- unique(hog$id)
k <- 1
for(i in 1:n){
	temp <- subset(hog, hog == h[i])
	p <- 0
	c <- 0
	counts <- c()
	for(j in 1:length(inds)){
		temp2 <- subset(temp, id == inds[j])
		if(all(is.na(temp2$genes))){
			p <- p + 1
			counts[j] <- 0
		} 
		else counts[j] <- sum(!is.na(unlist(temp2$genes)))
		c <- c + 1
		cnv[k,1] <- h[i]
		cnv[k,2] <- inds[j]
		cnv[k,3] <- sum(!is.na(unlist(temp2$genes)))
		k <- k + 1
		
	}
	p <- p / c
	if(p==0) pav$core <- pav$core + counts
	else if(p>0 & p<=0.85) pav$shell <- pav$shell + counts
	else pav$cloud <- pav$cloud + counts
}

#Add genes left out of the hogs (i.e., those unique for one assembly)
inds <- unique(all_genes$V1)
for(i in inds){
	temp <- subset(all_genes, V1 == i)
	temp2 <- subset(hog, id == i)
	pav$cloud[match(i,pav$id)] <- pav$cloud[match(i,pav$id)] + sum(!(temp$V2 %in% unlist(temp2$genes)))
	new_genes <- temp$V2[!(temp$V2 %in% unlist(temp2$genes))]
	for(j in new_genes){
		temp <- data.frame(hog=j,id=inds,n=0)
		temp$n[temp$id==i] <- 1
		cnv <- rbind(cnv,temp)
	}
}

#Plot PAV results
d <- melt(pav)
count <- data.frame(d %>% group_by(id) %>% summarise(sum(value)))
count <- count[order(count$sum.value.,decreasing=T),]
labels <- count$id
labels[labels == "TN1"] <- "TN1 (A. sterilis)"
labels[labels == "TN4"] <- "TN4 (A. sterilis)"
id_ord <- factor(d$id,levels=count$id, labels=labels)
var_ord <- factor(d$variable, levels=rev(c("core","shell","cloud")), labels=rev(c("Core","Shell","Cloud")))

ggplot(d,aes(x=id_ord,y=value,fill=var_ord))+
	geom_bar(position="stack",stat="identity")+
	scale_fill_manual(values=c("#FDE725FF","#00B395","#31688EFF"))+
	labs(x="Accession",y="Count")+
	theme(panel.background = element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		panel.border=element_blank(),
		axis.line=element_line(color="black",size=0.5),
		axis.text.x=element_text(size=10, hjust=1,vjust=0.5,angle=90,color="black"),
		axis.text.y=element_text(size=11,color="black"),
		axis.ticks.length=unit(.15, "cm"),
		axis.ticks=element_line(color="black",size=0.5),
		axis.title=element_text(size=12, color="black"),
		legend.title=element_blank(),
		legend.text=element_text(size=11, color="black"),
		legend.key=element_blank())


#Plot CNV results	
count <- data.frame(cnv %>% group_by(hog) %>% summarise(sum(n)))
names(count) <- c("hog","x")
count <- count[order(count$x,decreasing=F),]
r <- count[count$x==1,]
count[count$x==1,] <- r[sample(nrow(r)),]
ord_hog <- factor(cnv$hog, levels=count$hog)
labels <- unique(cnv$id)
labels[labels == "TN1"] <- "TN1 (A. sterilis)"
labels[labels == "TN4"] <- "TN4 (A. sterilis)"
colfunc <- colorRampPalette(c("white", "#31688EFF"))
colours <- c(colfunc(4),viridis(14,option="inferno", begin=0.4, end=0.92, direction=-1))
ggplot(data=cnv, aes(id, ord_hog, fill=n))+
	geom_tile(color="white")+
	scale_fill_gradientn(breaks=c(0,3,6,9,12,15),colors=colours, guide="legend")+
	labs(x="Accession",y="Hierarchical orthogroup", fill="Copy number")+
	scale_x_discrete(labels=labels)+
	theme(panel.background = element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		panel.border=element_blank(),
		axis.line=element_line(color="black",size=0.5),
		axis.text.x=element_text(size=10, hjust=1,vjust=0.5,angle=90,color="black"),
		axis.text.y=element_blank(),
		axis.ticks.length=unit(.15, "cm"),
		axis.ticks.x=element_line(color="black",size=0.5),
		axis.ticks.y=element_blank(),
		axis.title=element_text(size=12, color="black"),
		legend.title=element_text(size=12, color="black"),,
		legend.text=element_text(size=11, color="black"),
		legend.key=element_blank())

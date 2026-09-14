#Scripts for analysing the expression of GDSL genes
#Tuomas Hämälä 2026

library(vegan)
library(reshape2)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(viridis)

exp <- read.table("all_expression_TN1.txt", header=T) ##File listing expression data. Format: accession, gene, leaf, embryo, root, nod, panicle, caryopsis, ref gene in TN1

##RDA
#For RDA, data is compiled for each ref gene
inds <- unique(exp$id)
df <- NA
for(i in 1:length(inds)){
	temp <- subset(exp, id == inds[i])
	temp2 <- data.frame(temp[,c(10,3:8)])
	names(temp2) <- c("ref",paste0(inds[i],"_",colnames(temp)[3]),paste0(inds[i],"_",colnames(temp)[4]),paste0(inds[i],"_",colnames(temp)[5]),paste0(inds[i],"_",colnames(temp)[6]),paste0(inds[i],"_",colnames(temp)[7]),paste0(inds[i],"_",colnames(temp)[8]))
	if(i == 1) df <- temp2
	else df <- merge(df,temp2,by="ref")
}
d <- t(log10(df[,-1]+1))
colnames(d) <- df$ref
pred <- data.frame(geno=sapply(rownames(d), function(x){unlist(strsplit(as.character(x),"_"))[1]}),tissue=sapply(rownames(d), function(x){unlist(strsplit(as.character(x),"_"))[2]}))

#Fit partial RDA models
rda_full <- rda(d ~ geno+tissue,pred) #both genotype and tissue
rda_geno <- rda(d ~ geno+Condition(tissue),pred) #genotype only	
rda_tissue <- rda(d ~ tissue+Condition(geno),pred) #tissue only
RsquareAdj(rda_full)
RsquareAdj(rda_geno)
RsquareAdj(rda_tissue)

#Function for finding genotype-based outlier genes
find_outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)   
  o <- x[x < lims[1] | x > lims[2]]
  return(names(o))
}

#Plot RDA results from full model
scaling <- 2
rda <- as.data.frame(scores(rda_full, display="sites", scaling=scaling))
genes <- as.data.frame(scores(rda_full, display="species", scaling=scaling))
genes_geno <- scores(rda_geno, display="species")
outl <- c(find_outliers(genes_geno[,1],2),find_outliers(genes_geno[,2],2))
genes$outl <- ifelse(rownames(genes) %in% outl,1,0)
rda$tissue <- pred$tissue
xlab <- paste0("RDA1 (",round(summary(rda_full)$concont$importance[2]*100),"%)")
ylab <- paste0("RDA2 (",round(summary(rda_full)$concont$importance[5]*100),"%)")

ord <- factor(rda$tissue, levels=c("root","embryo","panicle","caryopsis","leaf","nod"),labels=c("Root","Embryo","Panicle","Caryopsis","Leaf","Internode"))

ggplot(rda, aes(RDA1, RDA2, color=ord))+
  geom_hline(yintercept=0,color="grey50",size=0.7,linetype="dashed")+
  geom_vline(xintercept=0,color="grey50",size=0.7,linetype="dashed")+
  geom_point(data=genes, aes(x=RDA1,y=RDA2,fill=as.factor(outl)), shape=21, size=2.5, color="NA", show.legend=F)+
  geom_point(size=5)+
  scale_color_manual(values=c("#F98400","#5BBCD6","#DD3226","#F2AD00","#00A08A","#748AA6"))+
  scale_fill_manual(values=c("grey55","grey20"))+
  labs(x=xlab, y=ylab, color="Tissue")+
  theme(panel.background = element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		panel.border=element_blank(),
		axis.line=element_line(color="black",size=0.5),
		axis.text=element_text(size=11,color="black"),
		axis.ticks.length=unit(.15, "cm"),
		axis.ticks=element_line(color="black",size=0.5),
		axis.title=element_text(size=12, color="black"),
		plot.title=element_text(size=14, color="black", hjust = 0.5),
		legend.text=element_text(size=11, color="black"),
		legend.title=element_text(size=12, color="black"),
		legend.key=element_blank(),
		aspect.ratio=1)


##SV-expression associations
exp <- read.table("exp_single_copy_TN1.txt", header=T) #Expression data from single-copy orthologs
ins <- read.table("INS.txt", header=T) #Insertions from SVGAP
del <- read.table("DEL.txt", header=T) #Deletions from SVGAP
gdsl <- read.table("TN1_gdsl.bed", header=T) #Locations of GDSL genes in TN1
temp <- rbind(ins,del)
temp$chr <- substr(temp$chr,5,nchar(temp$chr))
gr <- makeGRangesFromDataFrame(temp,keep.extra.columns=T)
grl <- makeGRangesFromDataFrame(gdsl,keep.extra.columns=T)
olaps <- findOverlaps(gr,grl,maxgap=10000) #Keep SVs within 10 kb of GDSL genes
q_matched <- gr[queryHits(olaps)]
mcols(q_matched) <- cbind(mcols(q_matched), mcols(grl)[subjectHits(olaps), ]) #Combine gene information with SV data
temp <- data.frame(q_matched)
names(temp) <- c(names(temp)[-ncol(temp)],"gene")
temp <- temp[,c(1:7,ncol(temp),8:(ncol(temp)-1))]

#Transform VCF genotypes into numeric genotypes
sv <- temp
n <- nrow(temp)
m <- ncol(temp)
for(i in 1:n){
	if(i==1|i==n|i%%100==0) cat("Line",i,"/",n,"\n")
	for(j in 9:m){
		g <- unlist(strsplit(temp[i,j],":"))[1]
		if(g == "0|0") ng <- 0
		else if(g == "1|1") ng <- 1
		else ng <- NA
		sv[i,j] <- ng
	}
}
mis <- apply(sv[,-c(1:8)],1,function(x){sum(is.na(x))/length(x)})
maf <- apply(sv[,-c(1:8)],1,function(x){sum(as.numeric(x),na.rm=T)/length(x)})
maf[maf>0.5] <- 1 - maf[maf>0.5]
filt_sv <- sv[mis<=0.2&maf>0.1,] #Filter for missing data and low MAF
filt_sv[,9:28] <- apply(filt_sv[,9:28],2,as.numeric)

#Calculate covariance matrix and perform PCA
temp <- filt_sv[,9:28]
p <- apply(temp,1,function(x)sum(x,na.rm=T)/length(x[!is.na(x)]))
n <- ncol(temp)
cov <- matrix(nrow=n,ncol=n)
for(i in 1:n){
	for(j in 1:i){
		cov[i,j] <- mean((temp[,i]-p)*(temp[,j]-p)/(p*(1-p)),na.rm=T)
		cov[j,i] <- cov[i,j]
	}	
}
pc <- prcomp(cov,scale=T)
pcs <- data.frame(id=names(temp),pc1=as.numeric(pc$x[,1]),pc2=as.numeric(pc$x[,2]),pc3=as.numeric(pc$x[,3])) #3 PCs used as cofactors

#Test for association between SVs and expression
d <- filt_sv
d$p <- NA
d$r2 <- 0
d$tis <- NA
n <- nrow(d)
for(i in 1:n){
	e <- subset(exp, TN1 == d$gene[i])
	if(nrow(e) == 0) next
	idx <- match(names(d)[9:28],e$id)
	idx2 <- match(names(d)[9:28],pcs$id)
	r2 <- 0
	for(j in 5:10){
		g <- as.numeric(d[i,9:28])
		p <- e[idx,j]
		p[is.na(g)] <- NA
		pc <- pcs[idx2,-1]
		pc[is.na(g),] <- NA
		m0 <- lm(log10(p+1)~pc$pc1+pc$pc2+pc$pc3)
		m1 <- lm(log10(p+1)~pc$pc1+pc$pc2+pc$pc3+g)
		p_val <- unlist(anova(m0,m1)[6])[2]
		if(deviance(m0) == 0 | deviance(m1) == 0) next
		r2 <- 1-deviance(m1)/deviance(m0)
		if(is.nan(r2)) next
		if(r2 > d$r2[i]){
			d$p[i] <- p_val
			d$r2[i] <- r2
			d$tis[i] <- names(e)[j]
		}
	}
}
d <- d[!is.na(d$p),]
d <- d[order(d$r2,decreasing=T),]
d$q <- qvalue(d$p)$qvalues

subset(d,q<0.05) #Significant associations

#Plot two genes with highest r2
e <- subset(exp, TN1=="AVESA.00700a.r1.3AG01907900.1")
idx <- match(names(d)[9:28],e$id)
g <- as.numeric(d[1,9:28])
p <- e[idx,10]
o <- data.frame(id=e$id[idx],g=g,p=p)
idx2 <- match(names(d)[9:28],pcs$id)
pc <- pcs[idx2,-1]
m0 <- lm(log(p+1)~pc$pc1+pc$pc2+pc$pc3)
m1 <- lm(log(p+1)~pc$pc1+pc$pc2+pc$pc3+g)
1-deviance(m1)/deviance(m0)
ggplot(o, aes(x=as.factor(g),y=log10(p+1),color=as.factor(g)))+
	geom_point(size=5, position=position_jitterdodge())+
	labs(x="Insertion genotype", y="Caryopsis expression", title="3AG01907900")+
	scale_y_continuous(limits=c(0,1.5))+
	theme(panel.background = element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		panel.border=element_blank(),
		axis.line=element_line(color="black",size=0.5),
		axis.text.x=element_text(size=11,color="black"),
		axis.text.y=element_text(size=11,color="black"),
		axis.ticks.length=unit(.15, "cm"),
		axis.ticks=element_line(color="black",size=0.5),
		axis.title=element_text(size=12, color="black"),
		plot.title=element_text(size=12,color="black"),
		legend.position="none")
		
e <- subset(exp, TN1=="AVESA.00700a.r1.6CG04833700.1")
idx <- match(names(d)[9:28],e$id)
g <- as.numeric(d[2,9:28])
p <- e[idx,7]
o <- data.frame(id=e$id[idx],g=g,p=p)
idx2 <- match(names(d)[9:28],pcs$id)
pc <- pcs[idx2,-1]
m0 <- lm(log(p+1)~pc$pc1+pc$pc2+pc$pc3)
m1 <- lm(log(p+1)~pc$pc1+pc$pc2+pc$pc3+o$g)
1-deviance(m1)/deviance(m0)
ggplot(o, aes(x=as.factor(g),y=log10(p+1),color=as.factor(g)))+
	geom_point(size=5, position=position_jitterdodge())+
	labs(x="Deletion genotype", y="Root expression", title="4CG03065080")+
	scale_y_continuous(limits=c(0,1))+
	theme(panel.background = element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		panel.border=element_blank(),
		axis.line=element_line(color="black",size=0.5),
		axis.text.x=element_text(size=11,color="black"),
		axis.text.y=element_text(size=11,color="black"),
		axis.ticks.length=unit(.15, "cm"),
		axis.ticks=element_line(color="black",size=0.5),
		axis.title=element_text(size=12, color="black"),
		plot.title=element_text(size=12,color="black"),
		legend.position="none")

##Average expression in each tissue
m <- melt(df)
m$id <- sapply(strsplit(as.character(m$variable), "_"), `[`, 1)
m$tissue <- sapply(strsplit(as.character(m$variable), "_"), `[`, 2)

mean <- as.data.frame(m %>% group_by(id,tissue) %>% summarise(value=mean(log10(value+1))))
ord <- factor(mean$tissue, levels=c("root","embryo","panicle","caryopsis","leaf","nod"),labels=c("Root","Embryo","Panicle","Caryopsis","Leaf","Internode"))

ggplot(mean, aes(ord,value,color=ord))+
	geom_point(size=6, position=position_jitterdodge(0.8))+
	labs(x="Tissue", y="Average expression level")+
	scale_y_continuous(limits=c(0.4,max(mean$value)))+
	scale_color_manual(values=c("#F98400","#5BBCD6","#DD3226","#F2AD00","#00A08A","#748AA6"))+
	theme(panel.background = element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		panel.border=element_blank(),
		axis.line=element_line(color="black",size=0.5),
		axis.text=element_text(size=11,color="black"),
		axis.ticks.length=unit(.15, "cm"),
		axis.ticks=element_line(color="black",size=0.5),
		axis.title=element_text(size=12, color="black"),
		plot.title=element_text(size=14, color="black", hjust = 0.5),
		legend.position="none")

##Tissue-specific expression
#Mean across all accessions
temp <- exp
temp[,3:8][temp[3:8] < 0.5] <- 0
temp[,3:8] <- log10(temp[,3:8]+1)
exp_mean <- data.frame(temp %>% group_by(ref) %>% summarise(leaf=mean(leaf,na.rm=T),embryo=mean(embryo,na.rm=T),root=mean(root,na.rm=T),nod=mean(nod,na.rm=T),panicle=mean(panicle,na.rm=T),caryopsis=mean(caryopsis,na.rm=T)))
exp_mean$tau <- apply(exp_mean[,-1], 1, function(x){sum(1-x/max(x))/(6-1)}) #Calculate tau-index
exp_mean$tissue <- apply(exp_mean[,2:7],1,function(x){names(exp_mean)[2:7][which.max(x)]}) #Which tissue has highest expression
o <- subset(exp_mean, tau >= 0.8) #Tissue-specific genes
tis_spe <- data.frame(table(o$tissue))
names(tis_spe) <- c("tissue","point")
tis_spe$low <- NA
tis_spe$high <- NA
boots <- rmultinom(1000,sum(tis_spe$point),tis_spe$point)
tis_spe[,3:4] <- t(apply(boots,1,quantile,probs=c(0.025,0.975))) #Define 95% CIs
ord_ts <- factor(tis_spe$tissue, levels=c("root","embryo","panicle","caryopsis","leaf","nod"),labels=c("Root","Embryo","Panicle","Caryopsis","Leaf","Internode"))

ggplot(tis_spe,aes(x=ord_ts,y=point,fill=ord_ts))+
	geom_bar(position="stack",stat="identity")+
	geom_errorbar(data=tis_spe, aes(ymin=low, ymax=high), position=position_dodge(width=0.6), stat="identity", width=0, size=0.7, color="black")+
	scale_fill_manual(values=brewer.pal(6, "Dark2"))+
	labs(x="Tissue",y="Tissue-specific genes")+
	theme(panel.background = element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		panel.border=element_blank(),
		axis.line=element_line(color="black",size=0.5),
		axis.text.x=element_text(size=11, hjust=1,vjust=1,angle=45,color="black"),
		axis.text.y=element_text(size=11,color="black"),
		axis.ticks.length=unit(.15, "cm"),
		axis.ticks=element_line(color="black",size=0.5),
		axis.title=element_text(size=12, color="black"),
		legend.position="none")

#Tissue-specificity in each accession
exp[,3:8][exp[3:8] < 0.5] <- 0
exp$tau <- apply(log10(exp[,3:8]+1), 1, function(x){sum(1-x/max(x))/(6-1)}) #Function for calculating tau-index
exp$tau[is.nan(exp$tau)] <- NA
for(i in 1:nrow(exp)){
	if(is.na(exp$tau[i]) == F) exp$max_tissue[i] <- names(exp)[3:8][which.max(exp[i,3:8])]
	else exp$max_tissue[i] <- NA
}
exp$max_tissue <- as.character(exp$max_tissue)
o <- subset(exp, tau >= 0.8)
tis_spe <- as.data.frame.matrix(table(o$max_tissue,o$id))
tis_spe$tissue <- row.names(tis_spe)
m <- melt(tis_spe)
ord_ts <- factor(m$tissue, levels=c("root","embryo","panicle","caryopsis","leaf","nod"),labels=c("Root","Embryo","Panicle","Caryopsis","Leaf","Internode"))

ggplot(m,aes(x=ord_ts,y=value,fill=variable))+
	geom_bar(position="dodge", stat="identity", color="black")+
	scale_fill_manual(values=viridis(20,option="turbo"))+
	labs(x="Tissue",y="Tissue-specific genes", fill="Accesssion")+
	theme(panel.background = element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		panel.border=element_blank(),
		axis.line=element_line(color="black",size=0.5),
		axis.text.x=element_text(size=11,color="black"),
		axis.text.y=element_text(size=11,color="black"),
		axis.ticks.length=unit(.15, "cm"),
		axis.ticks=element_line(color="black",size=0.5),
		axis.title=element_text(size=12, color="black"),
		legend.text=element_text(size=11,color="black"),
		legend.title=element_text(size=12, color="black"),
		legend.key=element_blank())

##Subgenome dominance
#Function for comparing observed and expected patterns
comp_triad <- function(x){
	ideal <- data.frame(bal=c(0.33,0.33,0.33),
		A_sup=c(0,0.5,0.5),
		C_sup=c(0.5,0,0.5),
		D_sup=c(0.5,0.5,0),
		A_dom=c(1,0,0),
		C_dom=c(0,1,0),
		D_dom=c(0,0,1))
	dist <- apply(ideal,2,function(y){sqrt(sum((x-y)^2))})
	return(names(dist)[which.min(dist)])
}

triads <- read.table("triads.txt", header=T) #Complete triads defined using hogs from OrthoFinder
inds <- unique(triads$id)
out <- data.frame(tissue=names(triads)[3:8],bal=0,A_sup=0,C_sup=0,D_sup=0,A_dom=0,C_dom=0,D_dom=0)
for(i in inds){
	temp <- subset(triads, id == i)
	hogs <- unique(temp$hog)
	for(h in hogs){
		temp2 <- subset(temp, hog == h)
		if(nrow(temp2) != 3) next
		temp2 <- temp2[order(temp2$subgenome),]
		for(j in 3:8){
			if(sum(temp2[,j]) < 0.5) next
			temp2[,j] <- log2(temp2[,j]+1)
			x <- comp_triad(temp2[,j]/sum(temp2[,j]))
			out[j-2,match(x,names(out))] <- out[j-2,match(x,names(out))] + 1
		}
	}
}

m <- reshape2::melt(out)
ord_sg <- factor(m$tissue, levels=c("root","embryo","panicle","caryopsis","leaf","nod"),labels=c("Root","Embryo","Panicle","Caryopsis","Leaf","Internode"))
ord_dom <- factor(m$variable, levels=c("A_dom","C_dom","D_dom","A_sup","C_sup","D_sup","bal"), labels=c("A dominant","C dominant","D dominant","A suppressed","C suppressed","D suppressed","Balanced"))

ggplot(m, aes(x=ord_sg, y=value, fill=ord_dom))+
	geom_bar(position="fill",stat="identity")+
	scale_fill_manual(values=paletteer_d("ggsci::legacy_tron"))+
	labs(x="Tissue",y="Proportion")+
	theme(panel.background = element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		panel.border=element_blank(),
		axis.line=element_line(color="black",size=0.5),
		axis.text.x=element_text(size=11, hjust=1,vjust=1,angle=45,color="black"),
		axis.text.y=element_text(size=11,color="black"),
		axis.ticks.length=unit(.15, "cm"),
		axis.ticks=element_line(color="black",size=0.5),
		axis.title=element_text(size=12, color="black"),
		legend.text=element_text(size=11, color="black"),
		legend.title=element_blank(),
		legend.key=element_blank())
	


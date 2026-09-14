#Scripts for estimating diversity metrics from multiple fasta format files
#Tuomas Hämälä 2026

library(phangorn)

#Function for extracting 0-fold and 4-fold sites
get_sites <- function(df, fold, focal=1){
	codon_table <- c("TTT","TTC","TTA","TTG","CTT","CTC","CTA","CTG","ATT","ATC","ATA","ATG","GTT","GTC","GTA","GTG","TCT","TCC","TCA","TCG","CCT","CCC","CCA","CCG","ACT","ACC","ACA","ACG","GCT","GCC","GCA","GCG","TAT","TAC","TAA","TAG","CAT","CAC","CAA","CAG","AAT","AAC","AAA","AAG","GAT","GAC","GAA","GAG","TGT","TGC","TGA","TGG","CGT","CGC","CGA","CGG","AGT","AGC","AGA","AGG","GGT","GGC","GGA","GGG")
	aa_table <- c("F","F","L","L","L","L","L","L","I","I","I","M","V","V","V","V","S","S","S","S","P","P","P","P","T","T","T","T","A","A","A","A","Y","Y","*","*","H","H","Q","Q","N","N","K","K","D","D","E","E","C","C","*","W","R","R","R","R","S","S","R","R","G","G","G","G")
	d <- data.frame(seq=toupper(df[,focal]),idx=1:nrow(df))
	d <- d[d$seq!="-",]
	if(nrow(d) %% 3 != 0){
		cat("ERROR: CDS length is not a multiple of 3!\n")
		return(NA)
	}
	sites <- integer()
	n <- nrow(d)
	for(i in seq(1,n,by=3)){
		codon <- d$seq[i:(i+2)]
		codon_idx <- match(paste(codon,collapse=""), codon_table)
		if(is.na(codon_idx)){
			cat("Warning: Codon",codon,"(sites",d$idx[i:(i+2)],")","not found\n")
			next
		}
		for(j in 1:3){
			aa <- NA
			k <- 1
			for(nuc in c("A","T","C","G")){
				temp <- codon
				temp[j] <- nuc
				codon_idx_temp <- match(paste(temp,collapse=""), codon_table)
				if(is.na(codon_idx_temp)) next
				aa[k] <- aa_table[codon_idx_temp]
				k <- k + 1
			}
			if(fold == 0 & length(unique(aa)) == 4) sites <- c(sites, d$idx[i+j-1])
			else if(fold == 4 & length(unique(aa)) == 1) sites <- c(sites, d$idx[i+j-1])
			else if(fold != 0 & fold != 4){
				cat("ERROR: fold can only be 0 or 4!\n")
				return(NA)
			}
		}
	}
	return(df[sites,])
}


#SFS without outgroups
est_sfs <- function(df,mis=0.2){
	sfs <- rep(0,ncol(df)+1)
	n <- nrow(df)
	for(i in 1:n){
		nucs <- as.character(df[i,])
		nucs_temp <- nucs[nucs!="-"]
		if(length(unique(nucs_temp)) > 2) next
		nmis <- sum(nucs=="-")
		if(nmis/length(nucs) > mis) next
		if(nmis >= 1){
			probs <- table(nucs[nucs!="-"])/sum(table(nucs[nucs!="-"]))
			nucs[nucs=="-"] <- sample(names(probs),nmis,prob=probs,replace=T)
		}
		o <- sample(nucs,1)
		x <- sum(nucs!=o)
		sfs[x+1] <- sfs[x+1] + 1
	}
	return(sfs)
}

#SFS with outgroups
est_sfs_k <- function(df,outgrp=1,mis=0.2){
	if(length(outgrp) > 1){
		match = apply(df[,outgrp], 1, function(x) all(x == x[1]) & all(!grepl("-",as.character(x))))
		df <- df[match,]
	}
	sfs <- rep(0,ncol(df)-length(outgrp)+1)
	n <- nrow(df)
	for(i in 1:n){
		o <- df[i,outgrp[1]]
		if(o=="-") next
		nucs <- as.character(df[i,-outgrp])
		nucs_temp <- nucs[nucs!="-"]
		if(length(unique(nucs_temp)) > 2) next
		nmis <- sum(nucs=="-")
		if(nmis/length(nucs) > mis) next
		if(nmis >= 1){
			probs <- table(nucs[nucs!="-"])/sum(table(nucs[nucs!="-"]))
			nucs[nucs=="-"] <- sample(names(probs),nmis,prob=probs,replace=T)
		}
		x <- sum(nucs!=o)
		sfs[x+1] <- sfs[x+1] + 1
	}
	return(sfs)
}


#Pairwise nucleotide diversity
est_pi <- function(sfs){
	n <- length(sfs)-1
	w <- seq(0,n)/n
	tP <- n/(n-1)*2*sum(sfs*w*(1-w))
	tP/sum(sfs)
}

#Watterson's theta
est_wat <- function(sfs){
	n <- length(sfs)-1
	S <- sum(sfs[-c(1,length(sfs))])
	tW <- S/sum(1/1:(n-1))
	tW/sum(sfs)
}

#Tajima's D
est_tajD <- function(sfs){
	n <- length(sfs)-1
	S <- sum(sfs[-c(1,length(sfs))])
	if(S == 0) return(NA)
	tW <- S/sum(1/1:(n-1))
	w <- seq(0,n)/n
	tP <- n/(n-1)*2*sum(sfs*w*(1-w))
	a1 <- sum(1/1:(n-1))
    a2 <- sum(1/(1:(n-1))^2)
    b1 <- (n+1)/(3.0*(n-1))
    b2 <- 2*((n*n)+n+3.0)/(9*n*(n-1))
    c1 <- b1-1/a1
    c2 <- b2-(n+2)/(a1*n)+(a2/(a1*a1))
    e1 <- c1/a1
    e2 <- c2/((a1*a1)+a2)
	(tP-tW)/sqrt(e1*S+e2*S*(S-1))
}

#Fay & Wu's H
est_fayH <- function(sfs) {
  n <- length(sfs)-1
  i <- 1:(n-1)
  S <- sum(sfs[-c(1,length(sfs))])
  if(S == 0) return(NA)
  tP <- sum((2*i*(n-i))/(n*(n-1))*sfs[-c(1,length(sfs))])
  tL <- sum((2*i^2)/(n*(n-1))*sfs[-c(1,length(sfs))])
  a2 <- sum(1/(1:(n-1))^2)
  var1 <- (n-2)/(6.0*(n-1))*S
  var2 <- ((18*n^2*(3*n+2)*a2-(88*n^3+9*n^2-13*n+6))/(9*n*(n-1)^2))*S^2
  H <- (tP-tL)/sqrt(var1+var2)
  if(is.nan(H)) return(NA)
  else return(H)
}

files <- list.files(pattern=".fas") #MFA files from PRANK
n <- length(files)
out <- data.frame(gene=integer(),nInd=integer(),pi=integer(),wat=integer(),tajD=integer(),fayH=integer(),piN=integer(),piS=integer())
j <- 1
for(i in 1:n){
	if(i==1|i==n|i%%10==0) cat("File",i,"/",n,"\n")
	df <-  data.frame(read.phyDat(files[i], format="fasta"))
	outgrp <- na.omit(match(c("TN1","TN4"),names(df))) #Outgroups
	nInd <- ncol(df)-length(outgrp)
	if(nInd < 5) next
	nsyn <- get_sites(df,0,match("TN1",names(df))) #Extract 0-fold sites
	syn <- get_sites(df,4,match("TN1",names(df))) #Extract 4-fold sites
	n_sfs <- est_sfs(nsyn[,-outgrp]) #0-fold sfs
	s_sfs <- est_sfs(syn[,-outgrp]) #4-fold sfs
	sfs <- est_sfs(df[,-outgrp]) #sfs without polarisation
	sfs_k <- est_sfs_k(df,outgrp) #polarised sfs
	out[j,1] <- substr(files[i],1,nchar(files[i])-14)
	out[j,2] <- nInd
	out[j,3] <- est_pi(sfs)
	out[j,4] <- est_wat(sfs)
	out[j,5] <- est_tajD(sfs)
	out[j,6] <- est_fayH(sfs_k)
	out[j,7] <- est_pi(n_sfs)
	out[j,8] <- est_pi(s_sfs)
	j <- j + 1
}


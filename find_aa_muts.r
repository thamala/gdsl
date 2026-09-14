#Script for finding amino acid substitutions from MFA files with at a specified derived allele frequency. Output is information required by PPVED.
#Tuomas Hämälä 2026

library(phangorn)

#Function requires a data frame, focal sequence that will be used to define the amino acid sequence, outgroups, derived allele frequency, and minimum number of samples 
find_muts <- function(df, focal=1, outgrp=1, min_freq=0, min_n=5){
	if(ncol(df)-length(outgrp) < min_n) return(NA)
	df <- df[df[,focal]!="-",]
	rownames(df) <- 1:nrow(df)
	save <- df
	if(length(outgrp) > 1){
		match = apply(df[,outgrp], 1, function(x) all(x == x[1]) & all(!grepl("-",as.character(x))))
		df <- df[match,]
	}
	n <- nrow(df)
	ok <- 0
	for(i in 1:n){
		o <- df[i,outgrp[1]]
		if(o=="-") next
		aa <- as.character(df[i,-outgrp])
		aa <- aa[aa!="-"]
		if(all(aa == o) == F){
			tab <- data.frame(table(aa))
			tab$Prop <- tab$Freq / sum(tab$Freq)
			tab <- tab[tab$aa != o,]
			if(ok == 0){
				out <- data.frame(ps=rownames(df)[i],wild=o,mut=tab$aa,freq=tab$Prop)
				ok <- 1
			}
			else out <- rbind(out, data.frame(ps=rownames(df)[i],wild=o,mut=tab$aa,freq=tab$Prop))
		}
	}
	out <- subset(out, freq >= min_freq)
	if(nrow(out) == 0) return(NA)
	rownames(out) <- 1:nrow(out)
	if(ok == 0) return(NA)
	seq <- paste(save[,focal],collapse="")
	return(list(seq,out))
}

#Read MFA files and print results
files <- list.files(pattern=".mfa")
for(f in files){
	df <-  data.frame(read.phyDat(f, format="fasta", type="AA"))
	outgrp <- na.omit(match(c("TN1","TN4"),names(df)))
	o <- find_muts(df,focal=1,outgrp=outgrp,min_freq=0.9)
	if(is.na(o[1]) == F) print(list(f,o))
}

## if working with notepad++, set the language to R if not automatically done so.

# set your working directory
setwd("/Users/daria/Downloads")


#################################  
# missingness per individual plot
#################################

# read in *.imiss file
imiss=read.table("test.imiss",h=T)  
plot(imiss$IID, imiss$F_MISS, main="Scatterplot sample missingness",
   xlab="sampleID", ylab="F_MISS", pch=19)  


   
#########################################  
# heterozygosity rate per individual plot
#########################################

# read in *.het file
het=read.table("test.het",h=T)
# calculate mean heterozygosity
het$meanHet = (het$N.NM. - het$O.HOM.)/het$N.NM.
# write as .*txt file in your working directory fyi 
write.table(het,file="GWA_raw_meanHet.txt",col.names=TRUE,row.names=FALSE,quote=F,sep="\t")

# scatterplot heterozygosity rate for all individuals
plot(het$IID, het$meanHet, main="Scatterplot heterozygosity",
   xlab="sampleID", ylab="meanHet", pch=19)
  

  
##################  
# miss-vs-het plot
##################

# read in *.het file
het=read.table("test.het",h=T)
# calculate mean heterozygosity
het$meanHet = (het$N.NM. - het$O.HOM.)/het$N.NM.

# read in *.imiss file
imiss=read.table("test.imiss",h=T)
imiss$logF_MISS = log10(imiss[,6])

# geneplotter is part of Bioconductor. If the Bioconductor package was not yet installed: 
source("https://bioconductor.org/biocLite.R")
biocLite("geneplotter")

library("geneplotter")
colors  <- densCols(imiss$logF_MISS,het$meanHet)
#pdf("miss-vs-het.pdf") # use this if you want the plot to be outputted as a pdf (or change accordingly if you want it to be a jpeg or other type of figure file)
plot(imiss$logF_MISS,het$meanHet, col=colors, xlim=c(-3,0),ylim=c(0,0.5),pch=20, xlab="Proportion of missing genotypes", ylab="Heterozygosity rate",axes=F)
axis(2,at=c(0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5),tick=T)
axis(1,at=c(-3,-2,-1,0),labels=c(0.001,0.01,0.1,1))
abline(h=mean(het$meanHet)-(2*sd(het$meanHet)),col="RED",lty=2)
abline(h=mean(het$meanHet)+(2*sd(het$meanHet)),col="RED",lty=2)
abline(v=-0.8860566, col="RED", lty=2)
# when working interactively, don't forget following line:
#dev.off()




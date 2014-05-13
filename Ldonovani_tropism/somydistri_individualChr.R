options <- commandArgs(trailingOnly = T)
FILEIN <- options[1]
print(FILEIN)

somydata0 <- read.table(FILEIN,  header=TRUE)
somydata1<-t(somydata0)
chromosomes <- colnames(somydata1)
color<-c("darkgreen", "red", "cyan", "blue")
libcodes <- rownames(somydata1)

 for (i in 1:ncol(somydata1)){
        chr<-chromosomes[i]
        print(chr)
        pngfile <- paste("LinJ.",chr, "_somy.png", sep="")
        print(pngfile)
        png(pngfile)
        title=paste("Chromosome LinJ.", chr, "\nDistribution of Somy across \n Sri Lankan isolates/strains")
        barplot(somydata1[,i], col=color, main=title, ylim=c(0,4), names.arg = libcodes, las=2)
	grid(nx=NA, ny=NULL, col="gray")
        dev.off()
 }



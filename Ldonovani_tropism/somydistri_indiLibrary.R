options <- commandArgs(trailingOnly = T)
FILEIN <- options[1]
print(FILEIN)

somydata1 <- read.table(FILEIN,  header=TRUE)
libs <- colnames(somydata1)
color<-c("darkgreen", "red", "cyan", "blue", "orange", "purple", "yellow", "brown", "limegreen", "black", "maroon", "gray", "#b0c4de", "#fa8072", "#ff00ff", "#ffc0cb", "#dcdcdc", "#008080", "#00fa9a" )
chrs <- rownames(somydata1)

 for (i in 1:ncol(somydata1)){
        lib<-libs[i]
        print(lib)
        pngfile <- paste(lib, "_somy.png", sep="")
        print(pngfile)
        png(pngfile)
        title=paste("Librarycode", lib, "\nDistribution of Somy across chromosomes")
        barplot(somydata1[,i], col=color, main=title, ylim=c(0,4), names.arg = chrs, las=2)
	grid(nx=NA, ny=NULL, col="gray")
        dev.off()
 }



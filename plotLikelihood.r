args = paste(commandArgs(trailingOnly=TRUE))
likefile = args[[1]]
outfile = args[[2]]
l = scan(likefile,quiet=TRUE)
png(outfile)
plot(l,type='l',ylab="log likelihood",xlab="iterations")
g <- dev.off()

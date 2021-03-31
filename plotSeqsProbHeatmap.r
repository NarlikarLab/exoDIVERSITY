library("RColorBrewer")
library("plotfunctions")
args=commandArgs(trailingOnly=TRUE)
infofile = args[1]
inpfile = args[2]
outfile = args[3]

seqprobs = read.table(inpfile)
info = read.table(infofile)
modes = as.integer(rownames(table(info[,2])))
modeOrder = seq(min(modes),max(modes))

lseqs = cbind(info[,2],seqprobs[,2:dim(seqprobs)[2]])
newdf = lseqs[which(lseqs[,1]==modeOrder[1]),]

rmo = max(modes):min(modes)
for (i in 2:length(modeOrder)){
    newdf = rbind(newdf,lseqs[which(lseqs[,1]==modeOrder[i]),])
}

df = newdf[,2:dim(newdf)[2]]
ord = order(rmo)
df = df[,ord]

shape = dim(df)
lenOfSeq = shape[2]
totalSeqs = shape[1]
tdf = t(df)
w=4
h=10
modes = table(newdf[,1])
pal = colorRampPalette(c("white","gray10"))(256)

png(outfile,w=w,h=h,units="in",res=600)
par(xpd=TRUE,cex.lab=0.7,cex.axis=0.7,tcl=-0.05,mgp=c(0.5,0.2,0),las=1)
image(1:dim(tdf)[1],1:dim(tdf)[2],tdf,col=pal,xlab="",ylab="",xaxt="n",yaxt="n",zlim=c(0,1))
axis(1,labels=modeOrder,at=seq(1,length(modeOrder)),col=NA,col.ticks=1,cex.axis=0.7)
ylabels = c()
ylabelposns = c()
j = (length(modeOrder)-1)
for (i in 1:(length(modeOrder)-1)){
    if (i==1){
        yval = strtoi(modes[[as.character(modeOrder[i])]])
        mid = yval/2
        lines(c(0.5,lenOfSeq+0.5),c(yval,yval),lwd=0.8)
        #text(x= -0.6,y=mid,labels=j,cex=0.6,col="black")
        ylabels[i] = j
        ylabelposns[i] = mid
    }
    else{
        mid = yval
        yval = yval + strtoi(modes[[as.character(modeOrder[i])]])
        mid = (mid+yval)/2
        lines(c(0.5,lenOfSeq+0.5),c(yval,yval),lwd=0.8)
        #text(x= -0.6,y=mid,labels=j,cex=0.6,col="black")
        ylabels[i] = j
        ylabelposns[i] = mid
    }
    j=j-1
}
i=i+1
mid = yval
yval = yval + strtoi(modes[[as.character(modeOrder[i])]])
mid = (mid+yval)/2
ylabels[i] = j
ylables = modeOrder
ylabelposns[i] = mid
#text(x=-0.6,y=mid,labels=j,cex=0.6,col="black")
axis(2,labels=ylabels,at=ylabelposns,col=NA,col.ticks=1,cex.axis=0.7)
for (i in 2:length(modeOrder)){
    lines(c((i-0.5),(i-0.5)),c(0,totalSeqs),lwd=0.5,lty=3,col="black")
}

par(cex=0.5)
gradientLegend(valRange=c(0,1),pos = c(1,-(0.04*totalSeqs),lenOfSeq,-(0.048*totalSeqs)),coords=T,side=1,color=pal,cex.label=0.7,seq=1)
g = dev.off()

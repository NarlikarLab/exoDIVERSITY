suppressMessages(require("RColorBrewer"))
suppressMessages(require("plotfunctions"))
args = commandArgs(trailingOnly=TRUE)
color = args[1]
rfile = args[2]
outfile = args[3]
modeorder = args[4:length(args)]


df = read.table(rfile)

modes = unlist(df[,1])
seqsInMode = table(modes)

newdf = df[which(df[,1]==modeorder[1]),]
for (i in 2:length(modeorder)){
     newdf = rbind(newdf,df[which(df[,1]==modeorder[i]),])
}
df = newdf[,2:dim(newdf)[2]]
### Binarize based on whole data median
#med = median(unlist(df))
#df[df>med] = 1
#df[df<=med] = 0

shape = dim(df)
totalSeqs = shape[1]
lenOfSeq = shape[2]

q = quantile(unlist(df),0.98)
df[df>=q]=q

tdf = t(df)
#w = 2.5
#h = 7
#pdf(outfile,height=h,width=w)
width = 4
height = 10
png(outfile,w=width,h=height,units="in",res=600)
par(xpd=TRUE,cex.lab=0.8,cex.axis=0.8,tcl= -0.05,mgp=c(1.5,0.05,0),las=1)
pal = colorRampPalette(c("ivory2",color))(q+1)
image(1:dim(tdf)[1],1:dim(tdf)[2],tdf,col=pal,xlab="",ylab="",axes=FALSE)
o1 = modeorder
modenames = modeorder
i=0
yval=0
j = length(modeorder)
if (length(o1)>1){
for (i in 1:(length(o1)-1)){
    if (i==1){
        yval = strtoi(seqsInMode[o1[i]])
        mid = yval/2
        lines(c(1,lenOfSeq),c(yval,yval),lwd=1)
        text(x=-5,y=mid,labels=paste(modenames[i],sep=""),cex=0.8,col="black")
    }
    else{
        mid= yval
        yval = yval+strtoi(seqsInMode[o1[[i]]])
        mid = (mid+yval)/2
        lines(c(1,lenOfSeq),c(yval,yval),lwd=1)
        text(x=-5,y=mid,labels=paste(modenames[i],sep=""),cex=0.8,col="black")
    }
    j= j-1
}
}

i=i+1
mid = yval
yval = yval+strtoi(seqsInMode[o1[i]])
mid = (mid+yval)/2
text(x=-5,y=mid,labels=paste(modenames[i],sep=""),cex=0.8,col="black")

par(cex=1)
gradientLegend(valRange=c(0,q),pos=c(10,-(0.06*totalSeqs),(lenOfSeq-10),-(0.08*totalSeqs)),coords=T,side=1,color=pal,cex.label=1,n.seg=1)
g = dev.off()

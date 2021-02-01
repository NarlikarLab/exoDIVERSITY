draw <- function(infile,outfile)
{
	df <- read.table(infile)

        newdf = df[which(df[,1]==modeorder[1]),]
        for (i in 2:length(modeorder)){
            newdf = rbind(newdf,df[which(df[,1]==modeorder[i]),])
        }
        modes = table(newdf[,1])
        df = newdf[,2:dim(newdf)[2]]
        shape = dim(df)
        lenOfSeq = shape[2]
        totalSeqs = shape[1]
        
	if(any(apply(df,2,function(x) any(x == 4))))
	{
		colorVector <- c("green","blue","orange","red","black")
		letters <- c("A","C","G","T","N")
	}
	else
	{
		colorVector <- c("green","blue","orange","red")
		letters <- c("A","C","G","T")
	}
	tdf = t(df)
        #w = 2.5
	#h = 7
	#pdf(outfile,height=h,width=w)
        height = 10
        width = 4
        png(outfile,w=width,h=height,units="in",res=600)
	par(xpd=TRUE,cex.lab=0.8,cex.axis=0.8,tcl = -0.05,mgp=c(0.5,0.05,0),las=1)
	
	image(1:dim(tdf)[1],1:dim(tdf)[2],tdf,col=colorVector,xlab="",ylab="",axes=FALSE)
        
        order1 = modeorder
        i=0
        yvalue=0
        modenames = order1
        j=length(order1)
        if (length(order1)>1){
	for (i in 1:(length(order1)-1))
	{
		if (i == 1)
		{
		 yvalue<-strtoi(modes[[order1[i]]])
                 mid = yvalue/2
		 lines(c(1,lenOfSeq),c(yvalue,yvalue),lwd=1)
                 text(x=-4.5,y=mid,labels=paste(modenames[i],sep=""),cex=0.8,col="black")
		}
		else
		{
                 mid = yvalue
		 yvalue <- yvalue + strtoi(modes[[order1[i]]])
                 mid = (mid+yvalue)/2
		 lines(c(1,lenOfSeq),c(yvalue,yvalue),lwd=1)
                 text(x=-4.5,y=mid,labels=paste(modenames[i],sep=""),cex=0.8,col="black")
		}
                j=j-1
	}
        }
        
        i=i+1
        mid = yvalue
	yvalue = yvalue + strtoi(modes[[order1[i] ]])
        
        mid = (mid+yvalue)/2
	text(x=-4.5,y=mid,labels=paste(modenames[i],sep=""),cex=0.8,col="black")
        
	legend(0,30,inset=c(0,-0.15),title="",legend=letters,fill=colorVector,horiz=TRUE,cex=1,bty="n")
	garbage <- dev.off()
}

array <- commandArgs(trailingOnly=TRUE)	
inputfile <- array[1]
outputfile <- array[2]
modeorder <- array[3:length(array)]
draw(inputfile,outputfile)

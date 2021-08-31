library(ggplot2)
meth.data <- read.table("MGMTe1methallsamples_wo.txt",sep="\t")
meth.data <- t(meth.data)
colnames(meth.data) <- meth.data[1,]
meth.data <- meth.data[-1,]
r.names <- meth.data[,1]
meth.data <- meth.data[,-1]
meth.data <- as.data.frame(apply(meth.data,2,as.numeric))
row.names(meth.data) <- r.names
e1.data <- meth.data

meth.data <- read.table("MGMTupmethallsamples_august_wo.txt",sep="\t")
meth.data <- t(meth.data)
colnames(meth.data) <- meth.data[1,]
meth.data <- meth.data[-1,]
r.names <- meth.data[,1]
r.names <- gsub("-","",r.names)
meth.data <- meth.data[,-1]
meth.data <- as.data.frame(apply(meth.data,2,as.numeric))
row.names(meth.data) <- r.names
up.data <- meth.data

meth.data <- read.table("MGMTi1allsamples_wo.txt",sep="\t")
meth.data <- t(meth.data)
colnames(meth.data) <- meth.data[1,]
meth.data <- meth.data[-1,]
r.names <- meth.data[,1]
meth.data <- meth.data[,-1]
meth.data <- as.data.frame(apply(meth.data,2,as.numeric))
row.names(meth.data) <- r.names
i1.data <- meth.data

int.samples <- intersect(colnames(e1.data),intersect(colnames(i1.data),colnames(up.data)))
e1.data <- e1.data[,int.samples]
up.data <- up.data[,int.samples]
i1.data <- i1.data[,int.samples]
mean.e1 <- apply(e1.data,2,mean)
mean.up <- apply(up.data,2,mean)
mean.i1 <- apply(i1.data,2,mean)
is.high <- mean.e1>0.3
is.low <- mean.e1<0.09
label <- rep("intermediately methylated",ncol(e1.data))
label[is.high] <- "highly methylated"
label[is.low] <- "lowly methylated"
to.plot <- data.frame(UpMeth=mean.up,Exon1=mean.e1,Intron1=mean.i1,Class=label)
#to.plot <- to.plot[-37,]
to.plot <- reshape2::melt(to.plot,id="Class")
colnames(to.plot) <- c("Class","Region","Methylation")
my_theme <- theme_bw()+theme(panel.grid=element_blank(),text=element_text(size=18,color="black"),
                  axis.ticks=element_line(color="black"),plot.title = element_text(size=5,color="black",hjust = 1,angle=90),
                  axis.text = element_text(size=15,color="black"))

plot <- ggplot(to.plot,aes(x=Region,y=100*Methylation,color=Class,group=Class))+geom_point()+geom_smooth()+
	my_theme+scale_color_manual(values=rev(c("#4575b4","#300040","#d73027")))
ggsave("overall_description.pdf",plot)

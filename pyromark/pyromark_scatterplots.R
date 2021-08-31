library(ggplot2)
cpg.anno <- read.csv("CpG_anno_e1.csv",header=T)
row.names(cpg.anno) <- cpg.anno$CpG
cpg.anno <- cpg.anno[,-1]
cols.rows <- list(Esteller=c(forward="#469990",reverse="#ffe119"),
	Bady=c("Bady"="#800000"),
	Pyromark=c(Pyromark="#4363d8"),
	Felsberg=c(forward="#469990",reverse="#ffe119"),
	EPIC=c(EPIC="#808000"),
	A450k=c(A450k="#000075"))
my_theme <- theme_bw()+theme(panel.grid=element_blank(),text=element_text(size=18,color="black"),
                  axis.ticks=element_line(color="black"),plot.title = element_text(size=18,color="black",hjust = .5),
                  axis.text = element_text(size=15,color="black"))

te.meth <- read.table("MGMTe1methallsamples_august_wo.txt",sep="\t")
te.meth <- t(te.meth)
colnames(te.meth) <- te.meth[1,]
te.meth <- te.meth[-1,]
r.names <- te.meth[,1]
te.meth <- te.meth[,-1]
te.meth <- as.data.frame(apply(te.meth,2,as.numeric))
row.names(te.meth) <- r.names
row.names(cpg.anno) <- row.names(te.meth)

pyro.res <- read.csv("pyromark_methylation_august.csv")
row.names(pyro.res) <- pyro.res[,1]
pyro.res <- pyro.res[,-1]
pyro.res <- t(pyro.res)
s.amples <- intersect(colnames(pyro.res),colnames(te.meth))
te.meth <- te.meth[paste0("CpG.",10:13),s.amples]
pyro.res <- pyro.res[,s.amples]

avg.miseq <- apply(te.meth,2,mean,na.rm=T)
to.plot <- data.frame(Pyromark=pyro.res[5,],MiSeq=avg.miseq)
cori <- round(cor(na.omit(to.plot$Pyromark),na.omit(to.plot$MiSeq)),3)
cor.p <- cor.test(na.omit(to.plot)$Pyromark,na.omit(to.plot)$MiSeq)$p.value
plot <- ggplot(to.plot,aes(x=Pyromark,y=MiSeq))+geom_point()+geom_smooth(method="lm")+my_theme+xlim(-0.01,1.01)+ylim(-0.01,1.01)+geom_text(label=paste0("Pearson correlation: ",cori),x=0.1,y=0.5)
ggsave("average.pdf",plot)

for(i in 1:nrow(te.meth)){
	to.plot <- data.frame(Pyromark=pyro.res[i,],MiSeq=unlist(te.meth[i,,drop=T]))
	cori <- round(cor(na.omit(to.plot)$Pyromark,na.omit(to.plot)$MiSeq),3)
	cor.p <- c(cor.p,cor.test(na.omit(to.plot)$Pyromark,na.omit(to.plot)$MiSeq)$p.value)
	plot <- ggplot(to.plot,aes(x=Pyromark,y=MiSeq))+geom_point()+geom_smooth(method="lm")+my_theme+xlim(-0.01,1.01)+ylim(-0.01,1.01)+geom_text(label=paste0("Pearson correlation: ",cori),x=0.1,y=0.5)
	ggsave(paste0("CpG",9+i,".pdf"),plot)
}
names(cor.p) <- c("Average",row.names(te.meth))

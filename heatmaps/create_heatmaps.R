# exon1
library(pheatmap)
cpg.anno <- read.csv("CpG_anno_e1.csv",header=T)
row.names(cpg.anno) <- cpg.anno$CpG
cpg.anno <- cpg.anno[,-1]
cols.rows <- list(Esteller=c(forward="#469990",reverse="#ffe119"),
	Bady=c("Bady"="#800000"),
	Pyromark=c(Pyromark="#4363d8"),
	Felsberg=c(forward="#469990",reverse="#ffe119"),
	EPIC=c(EPIC="#808000"),
	A450k=c(A450k="#000075"))

s.anno <- read.csv('sample_ids_new.csv')
te.meth <- read.table("MGMTe1methallsamples.txt",sep="\t")
te.meth <- t(te.meth)
colnames(te.meth) <- te.meth[1,]
te.meth <- te.meth[-1,]
r.names <- te.meth[,1]
te.meth <- te.meth[,-1]
te.meth <- as.data.frame(apply(te.meth,2,as.numeric))
row.names(te.meth) <- r.names
matchi <- match(colnames(te.meth),s.anno$ID)
colnames(te.meth) <- s.anno$Spalte1[matchi]
te.meth <- te.meth[,!(colnames(te.meth)%in%c('#81','#82','#83'))]
#cpg.anno$Context <- gsub("[[:punct:]].","",row.names(te.meth))
#cpg.anno$Context <- gsub("[0-9]","",cpg.anno$Context)
row.names(cpg.anno) <- row.names(te.meth)

png("MGMTe1_heatmap.png")
pheatmap(t(te.meth), annotation_col = cpg.anno[,-c(4:6)], annotation_colors=cols.rows,
	cluster_cols=F,fontsize_col=10,fontsize_row=6)
dev.off()
pdf("MGMTe1_heatmap.pdf")
pheatmap(t(te.meth), annotation_col = cpg.anno[,-c(4:6)], annotation_colors=cols.rows,
	cluster_cols=F,fontsize_col=10,fontsize_row=6)
dev.off()

te.meth <- te.meth[,sort(colnames(te.meth))]
clust <- hclust(dist(t(te.meth)))
ord <- clust$order

# intron1
library(pheatmap)
cpg.anno <- read.csv("CpG_anno_i1.csv",header=T)
row.names(cpg.anno) <- cpg.anno$CpG
cpg.anno <- cpg.anno[,-1]
cols.rows <- list(Esteller=c(forward="#469990",reverse="#ffe119"),
	Bady=c("Bady"="#800000"),
	Pyromark=c(Pyromark="#4363d8"),
	Felsberg=c(forward="#469990",reverse="#ffe119"),
	EPIC=c(EPIC="#808000"),
	A450k=c(A450k="#000075"))

i.meth <- read.table("MGMTi1allsamples.txt",sep="\t")
i.meth <- t(i.meth)
colnames(i.meth) <- i.meth[1,]
i.meth <- i.meth[-1,]
r.names <- i.meth[,1]
i.meth <- i.meth[,-1]
i.meth <- as.data.frame(apply(i.meth,2,as.numeric))
row.names(i.meth) <- r.names
matchi <- match(colnames(i.meth),s.anno$ID)
colnames(i.meth) <- s.anno$Spalte1[matchi]
i.meth <- i.meth[,!(colnames(i.meth)%in%c('#81','#82','#83'))]
#cpg.anno$Context <- gsub("[[:punct:]].","",row.names(te.meth))
#cpg.anno$Context <- gsub("[0-9]","",cpg.anno$Context)
row.names(cpg.anno) <- row.names(i.meth)
i.meth <- i.meth[,sort(colnames(i.meth))]
i.meth <- i.meth[,ord]

png("MGMTi1_heatmap.png")
pheatmap(t(i.meth), annotation_col = cpg.anno[,-3], annotation_colors=cols.rows,
	cluster_cols=F,cluster_rows=F,fontsize_col=10,fontsize_row=6)
dev.off()
pdf("MGMTi1_heatmap.pdf")
pheatmap(t(i.meth), annotation_col = cpg.anno[,-3], annotation_colors=cols.rows,
	cluster_cols=F,cluster_rows=F,fontsize_col=10,fontsize_row=6)
dev.off()

# upmeth
library(pheatmap)
cpg.anno <- read.csv("CpG_anno_up.csv",header=T)
row.names(cpg.anno) <- cpg.anno$CpG
cpg.anno <- cpg.anno[,-1]
cols.rows <- list(Esteller=c(forward="#469990",reverse="#ffe119"),
	Bady=c("Bady"="#800000"),
	Pyromark=c(Pyromark="#4363d8"),
	Felsberg=c(forward="#469990",reverse="#ffe119"),
	EPIC=c(EPIC="#808000"),
	A450k=c(A450k="#000075"))

up.meth <- read.table("MGMTupmethallsamplesneu.txt",sep="\t")
up.meth <- t(up.meth)
colnames(up.meth) <- up.meth[1,]
up.meth <- up.meth[-1,]
r.names <- up.meth[,1]
up.meth <- up.meth[,-1]
up.meth <- as.data.frame(apply(up.meth,2,as.numeric))
row.names(up.meth) <- r.names
matchi <- match(colnames(up.meth),s.anno$ID)
colnames(up.meth) <- s.anno$Spalte1[matchi]
up.meth <- up.meth[,!(colnames(up.meth)%in%c('#81','#82','#83'))]

#cpg.anno$Context <- gsub("[[:punct:]].","",row.names(te.meth))
#cpg.anno$Context <- gsub("[0-9]","",cpg.anno$Context)
row.names(cpg.anno) <- row.names(up.meth)
up.meth <- up.meth[,sort(colnames(up.meth))]
up.meth <- up.meth[,ord]

png("MGMTup_heatmap.png")
pheatmap(t(up.meth), annotation_col = cpg.anno, annotation_colors=cols.rows,
         cluster_cols=F,cluster_rows=F,fontsize_col=10,fontsize_row=5)
dev.off()
pdf("MGMTup_heatmap.pdf")
pheatmap(t(up.meth), annotation_col = cpg.anno, annotation_colors=cols.rows,
         cluster_cols=F,cluster_rows=F,fontsize_col=10,fontsize_row=5)
dev.off()


library(survival)
library(RnBeads)
library(survminer)
all.p.vals <- c()
meth.data <- read.table("MGMTupmethallsamples_august_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("sample_annotation_reduced.csv")
conf.factors <- c("Sex","Age","PFS")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
p.vals <- apply(meth.data,1,function(cpg){
	used.dat <- data.frame(Methylation=cpg,ph)
	mod <- tryCatch(coxph(form,used.dat),error=function(e)return(e))
	if(inherits(mod,"error")){
		NA
	}else{
		summary(mod)$coefficient["Methylation","Pr(>|z|)"]
	}
})
all.p.vals <- c(all.p.vals,p.vals)

library(survival)
library(RnBeads)
library(survminer)
meth.data <- read.table("MGMTe1methallsamples_august_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("sample_annotation_reduced.csv")
conf.factors <- c("Sex","Age","PFS")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
p.vals <- apply(meth.data,1,function(cpg){
	used.dat <- data.frame(Methylation=cpg,ph)
	mod <- tryCatch(coxph(form,used.dat),error=function(e)return(e))
	if(inherits(mod,"error")){
		NA
	}else{
		summary(mod)$coefficient["Methylation","Pr(>|z|)"]
	}
})
all.p.vals <- c(all.p.vals,p.vals)

library(survival)
library(RnBeads)
library(survminer)
meth.data <- read.table("MGMTi1allsamples_august_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("sample_annotation_reduced.csv")
conf.factors <- c("Sex","Age","PFS")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
p.vals <- apply(meth.data,1,function(cpg){
	used.dat <- data.frame(Methylation=cpg,ph)
	mod <- tryCatch(coxph(form,used.dat),error=function(e)return(e))
	if(inherits(mod,"error")){
		NA
	}else{
		summary(mod)$coefficient["Methylation","Pr(>|z|)"]
	}
})
all.p.vals <- c(all.p.vals,p.vals)

library(survival)
library(RnBeads)
library(survminer)
meth.data <- read.table("MGMTupmethallsamples_august_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("sample_annotation_reduced.csv")
conf.factors <- c("Sex","Age","PFS","path")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
dat <- data.frame(ph,Methylation=ifelse(meth.data["CpG.21",]>0.2,"high","low"))
mod <- survfit(Surv(PFS,status)~Methylation,data=dat)
ggsurvplot(mod,data=dat,pval=T)

#dat <- read.csv("/DEEP_fhgfs/projects/mscherer/data/other/MGMT_Homburg/Alea/final_run/survival_NGS_CpG20.csv")
#mod <- survfit(Surv(time,status=="Dead")~Methylation>0.2,data=dat)
#ggsurvplot(mod,data=dat,pval=T)

all.p.vals.wo <- all.p.vals
library(survival)
library(RnBeads)
library(survminer)
all.p.vals <- c()
meth.data <- read.table("MGMTupmethallsamplesneu_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("sample_annotation_reduced_with.csv")
conf.factors <- c("Sex","Age","PFS","path")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
p.vals <- apply(meth.data,1,function(cpg){
	used.dat <- data.frame(Methylation=cpg,ph)
	mod <- tryCatch(coxph(form,used.dat),error=function(e)return(e))
	if(inherits(mod,"error")){
		NA
	}else{
		summary(mod)$coefficient["Methylation","Pr(>|z|)"]
	}
})
all.p.vals <- c(all.p.vals,p.vals)

library(survival)
library(RnBeads)
library(survminer)
meth.data <- read.table("MGMTe1methallsamples_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("sample_annotation_reduced_with.csv")
conf.factors <- c("Sex","Age","PFS","path")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
p.vals <- apply(meth.data,1,function(cpg){
	used.dat <- data.frame(Methylation=cpg,ph)
	mod <- tryCatch(coxph(form,used.dat),error=function(e)return(e))
	if(inherits(mod,"error")){
		NA
	}else{
		summary(mod)$coefficient["Methylation","Pr(>|z|)"]
	}
})
all.p.vals <- c(all.p.vals,p.vals)

library(survival)
library(RnBeads)
library(survminer)
meth.data <- read.table("MGMTi1allsamples_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("sample_annotation_reduced_with.csv")
conf.factors <- c("Sex","Age","PFS","path")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
p.vals <- apply(meth.data,1,function(cpg){
	used.dat <- data.frame(Methylation=cpg,ph)
	mod <- tryCatch(coxph(form,used.dat),error=function(e)return(e))
	if(inherits(mod,"error")){
		NA
	}else{
		summary(mod)$coefficient["Methylation","Pr(>|z|)"]
	}
})
all.p.vals <- c(all.p.vals,p.vals)

my_theme <- theme_bw()+theme(panel.grid=element_blank(),text=element_text(size=18,color="black"),
                  axis.ticks=element_line(color="black"),plot.title = element_text(size=18,color="black",hjust = .5),
                  axis.text = element_text(size=15,color="black"),
		axis.text.y = element_text(size=12,color="black"))
to.plot <- data.frame(Without=all.p.vals.wo,With=all.p.vals)
to.plot$CpG <- factor(c(paste0("CpG.",1:22),paste("CpG.",1:27)),levels=rev(c(paste0("CpG.",1:22),paste("CpG.",1:27))))
to.plot <- reshape2::melt(to.plot,id="CpG")
colnames(to.plot) <- c("CpG","Outliers","PValue")
plot <- ggplot(to.plot,aes(x=Outliers,y=CpG,fill=-log10(PValue)))+geom_tile()+my_theme+scale_fill_gradient(low='white',high='firebrick2')
ggsave("all_p_vals_heatmap.pdf",plot)

plot <- ggplot(to.plot,aes(x=Outliers,y=CpG,fill=ifelse(PValue>0.05,"white","red")))+geom_tile()+my_theme+theme(legend.position="none")
ggsave("all_p_vals_binary.pdf",plot)

# Sex
library(survival)
library(RnBeads)
library(survminer)
all.p.vals <- c()
meth.data <- read.table("MGMTupmethallsamplesneu_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("sample_annotation_reduced.csv")
conf.factors <- c("Sex","Age","PFS","path")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
p.vals <- apply(meth.data,1,function(cpg){
	used.dat <- data.frame(Methylation=cpg,ph)
	mod <- tryCatch(coxph(form,used.dat),error=function(e)return(e))
	if(inherits(mod,"error")){
		NA
	}else{
		summary(mod)$coefficient["Sexw","Pr(>|z|)"]
	}
})
all.p.vals <- c(all.p.vals,p.vals)

library(survival)
library(RnBeads)
library(survminer)
meth.data <- read.table("MGMTe1methallsamples_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("sample_annotation_reduced.csv")
conf.factors <- c("Sex","Age","PFS","path")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
p.vals <- apply(meth.data,1,function(cpg){
	used.dat <- data.frame(Methylation=cpg,ph)
	mod <- tryCatch(coxph(form,used.dat),error=function(e)return(e))
	if(inherits(mod,"error")){
		NA
	}else{
		summary(mod)$coefficient["Sexw","Pr(>|z|)"]
	}
})
all.p.vals <- c(all.p.vals,p.vals)

library(survival)
library(RnBeads)
library(survminer)
meth.data <- read.table("MGMTi1allsamples_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("sample_annotation_reduced.csv")
conf.factors <- c("Sex","Age","PFS","path")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
p.vals <- apply(meth.data,1,function(cpg){
	used.dat <- data.frame(Methylation=cpg,ph)
	mod <- tryCatch(coxph(form,used.dat),error=function(e)return(e))
	if(inherits(mod,"error")){
		NA
	}else{
		summary(mod)$coefficient["Sexw","Pr(>|z|)"]
	}
})
all.p.vals <- c(all.p.vals,p.vals)

all.p.vals.wo <- all.p.vals
library(survival)
library(RnBeads)
library(survminer)
all.p.vals <- c()
meth.data <- read.table("MGMTupmethallsamplesneu_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("ample_annotation_reduced_with.csv")
conf.factors <- c("Sex","Age","PFS","path")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
p.vals <- apply(meth.data,1,function(cpg){
	used.dat <- data.frame(Methylation=cpg,ph)
	mod <- tryCatch(coxph(form,used.dat),error=function(e)return(e))
	if(inherits(mod,"error")){
		NA
	}else{
		summary(mod)$coefficient["Sexw","Pr(>|z|)"]
	}
})
all.p.vals <- c(all.p.vals,p.vals)

library(survival)
library(RnBeads)
library(survminer)
meth.data <- read.table("MGMTe1methallsamples_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("sample_annotation_reduced_with.csv")
conf.factors <- c("Sex","Age","PFS","path")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
p.vals <- apply(meth.data,1,function(cpg){
	used.dat <- data.frame(Methylation=cpg,ph)
	mod <- tryCatch(coxph(form,used.dat),error=function(e)return(e))
	if(inherits(mod,"error")){
		NA
	}else{
		summary(mod)$coefficient["Sexw","Pr(>|z|)"]
	}
})
all.p.vals <- c(all.p.vals,p.vals)

library(survival)
library(RnBeads)
library(survminer)
meth.data <- read.table("MGMTi1allsamples_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("sample_annotation_reduced_with.csv")
conf.factors <- c("Sex","Age","PFS","path")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
p.vals <- apply(meth.data,1,function(cpg){
	used.dat <- data.frame(Methylation=cpg,ph)
	mod <- tryCatch(coxph(form,used.dat),error=function(e)return(e))
	if(inherits(mod,"error")){
		NA
	}else{
		summary(mod)$coefficient["Sexw","Pr(>|z|)"]
	}
})
all.p.vals <- c(all.p.vals,p.vals)

my_theme <- theme_bw()+theme(panel.grid=element_blank(),text=element_text(size=18,color="black"),
                  axis.ticks=element_line(color="black"),plot.title = element_text(size=18,color="black",hjust = .5),
                  axis.text = element_text(size=15,color="black"),
		axis.text.y = element_text(size=12,color="black"))
to.plot <- data.frame(Without=all.p.vals.wo,With=all.p.vals)
to.plot$CpG <- factor(c(paste0("CpG.",1:22),paste("CpG.",1:27)),levels=rev(c(paste0("CpG.",1:22),paste("CpG.",1:27))))
to.plot <- reshape2::melt(to.plot,id="CpG")
colnames(to.plot) <- c("CpG","Outliers","PValue")
plot <- ggplot(to.plot,aes(x=Outliers,y=CpG,fill=-log10(PValue)))+geom_tile()+my_theme+scale_fill_gradient(low='white',high='firebrick2')
ggsave("all_p_vals_heatmap_sex.pdf",plot)

plot <- ggplot(to.plot,aes(x=Outliers,y=CpG,fill=ifelse(PValue>0.05,"white","red")))+geom_tile()+my_theme+theme(legend.position="none")
ggsave("all_p_vals_binary_sex.pdf",plot)

# Age
library(survival)
library(RnBeads)
library(survminer)
all.p.vals <- c()
meth.data <- read.table("MGMTupmethallsamplesneu_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("sample_annotation_reduced.csv")
conf.factors <- c("Sex","Age","PFS","path")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
p.vals <- apply(meth.data,1,function(cpg){
	used.dat <- data.frame(Methylation=cpg,ph)
	mod <- tryCatch(coxph(form,used.dat),error=function(e)return(e))
	if(inherits(mod,"error")){
		NA
	}else{
		summary(mod)$coefficient["Age","Pr(>|z|)"]
	}
})
all.p.vals <- c(all.p.vals,p.vals)

library(survival)
library(RnBeads)
library(survminer)
meth.data <- read.table("MGMTe1methallsamples_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("sample_annotation_reduced.csv")
conf.factors <- c("Sex","Age","PFS","path")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
p.vals <- apply(meth.data,1,function(cpg){
	used.dat <- data.frame(Methylation=cpg,ph)
	mod <- tryCatch(coxph(form,used.dat),error=function(e)return(e))
	if(inherits(mod,"error")){
		NA
	}else{
		summary(mod)$coefficient["Age","Pr(>|z|)"]
	}
})
all.p.vals <- c(all.p.vals,p.vals)

library(survival)
library(RnBeads)
library(survminer)
meth.data <- read.table("MGMTi1allsamples_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("sample_annotation_reduced.csv")
conf.factors <- c("Sex","Age","PFS","path")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
p.vals <- apply(meth.data,1,function(cpg){
	used.dat <- data.frame(Methylation=cpg,ph)
	mod <- tryCatch(coxph(form,used.dat),error=function(e)return(e))
	if(inherits(mod,"error")){
		NA
	}else{
		summary(mod)$coefficient["Age","Pr(>|z|)"]
	}
})
all.p.vals <- c(all.p.vals,p.vals)

all.p.vals.wo <- all.p.vals
library(survival)
library(RnBeads)
library(survminer)
all.p.vals <- c()
meth.data <- read.table("MGMTupmethallsamplesneu_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("sample_annotation_reduced_with.csv")
conf.factors <- c("Sex","Age","PFS","path")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
p.vals <- apply(meth.data,1,function(cpg){
	used.dat <- data.frame(Methylation=cpg,ph)
	mod <- tryCatch(coxph(form,used.dat),error=function(e)return(e))
	if(inherits(mod,"error")){
		NA
	}else{
		summary(mod)$coefficient["Age","Pr(>|z|)"]
	}
})
all.p.vals <- c(all.p.vals,p.vals)

library(survival)
library(RnBeads)
library(survminer)
meth.data <- read.table("MGMTe1methallsamples_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("sample_annotation_reduced_with.csv")
conf.factors <- c("Sex","Age","PFS","path")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
p.vals <- apply(meth.data,1,function(cpg){
	used.dat <- data.frame(Methylation=cpg,ph)
	mod <- tryCatch(coxph(form,used.dat),error=function(e)return(e))
	if(inherits(mod,"error")){
		NA
	}else{
		summary(mod)$coefficient["Age","Pr(>|z|)"]
	}
})
all.p.vals <- c(all.p.vals,p.vals)

library(survival)
library(RnBeads)
library(survminer)
meth.data <- read.table("MGMTi1allsamples_wo.txt",sep='\t')
row.names(meth.data) <- meth.data[,1]
meth.data <- meth.data[,-1]
colnames(meth.data) <- as.character(unlist(meth.data[1,]))
meth.data <- meth.data[-1,]
meth.data <- t(meth.data)
meth.data <- apply(meth.data,c(1,2),as.numeric)
ph <- read.csv("sample_annotation_reduced_with.csv")
conf.factors <- c("Sex","Age","PFS","path")
row.names(ph) <- ph[,1]
ph <- ph[,-1]
ph$status <- rep(1,nrow(ph))
form <- as.formula(Surv(ph$PFS,status)~Methylation+Sex+Age,path)
s.int <- intersect(colnames(meth.data),row.names(ph))
ph <- ph[s.int,]
meth.data <- meth.data[,s.int]
ph$PFS <- as.numeric(as.character(ph$PFS))
p.vals <- apply(meth.data,1,function(cpg){
	used.dat <- data.frame(Methylation=cpg,ph)
	mod <- tryCatch(coxph(form,used.dat),error=function(e)return(e))
	if(inherits(mod,"error")){
		NA
	}else{
		summary(mod)$coefficient["Age","Pr(>|z|)"]
	}
})
all.p.vals <- c(all.p.vals,p.vals)

my_theme <- theme_bw()+theme(panel.grid=element_blank(),text=element_text(size=18,color="black"),
                  axis.ticks=element_line(color="black"),plot.title = element_text(size=18,color="black",hjust = .5),
                  axis.text = element_text(size=15,color="black"),
		axis.text.y = element_text(size=12,color="black"))
to.plot <- data.frame(Without=all.p.vals.wo,With=all.p.vals)
to.plot$CpG <- factor(c(paste0("CpG.",1:22),paste("CpG.",1:27)),levels=rev(c(paste0("CpG.",1:22),paste("CpG.",1:27))))
to.plot <- reshape2::melt(to.plot,id="CpG")
colnames(to.plot) <- c("CpG","Outliers","PValue")
plot <- ggplot(to.plot,aes(x=Outliers,y=CpG,fill=-log10(PValue)))+geom_tile()+my_theme+scale_fill_gradient(low='white',high='firebrick2')
ggsave("all_p_vals_heatmap_age.pdf",plot)

plot <- ggplot(to.plot,aes(x=Outliers,y=CpG,fill=ifelse(PValue>0.05,"white","red")))+geom_tile()+my_theme+theme(legend.position="none")
ggsave("all_p_vals_binary_age.pdf",plot)


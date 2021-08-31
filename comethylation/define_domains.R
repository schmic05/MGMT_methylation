########################################################################
# define_domains.R
# This file defines co-methylation domains for the MGMT exon1, intron 1,
# and upstream region accordingly

library(eDMA)
meth.data <- read.table("MGMTupmethallsamples_august_wo.txt",sep="\t")
meth.data <- t(meth.data)
colnames(meth.data) <- meth.data[1,]
meth.data <- meth.data[-1,]
r.names <- meth.data[,1]
r.names <- gsub("-","",r.names)
meth.data <- meth.data[,-1]
meth.data <- as.data.frame(apply(meth.data,2,as.numeric))
row.names(meth.data) <- r.names
mean.cpgs <- apply(meth.data,1,mean)
var.cpgs <- apply(meth.data,1,var)
dist.mean <- c()
dist.var <- c()
for(i in 1:length(mean.cpgs)-1){
	dist.mean <- c(dist.mean,abs(mean.cpgs[i]-mean.cpgs[i+1]))
	dist.var <- c(dist.var,abs(var.cpgs[i]-var.cpgs[i+1]))
}
thres.holds <- (dist.mean>quantile(dist.mean,.9))|(dist.var>quantile(dist.var,.9))
thres.holds[length(thres.holds)+1] <- TRUE

all.domains <- list()
dom.class <- list()
for(i in 0:(length(which(thres.holds))-1)){
	if(i==0){
		my.domain <- meth.data[0:which(thres.holds)[1],,drop=F]
	}else{
		my.domain <- meth.data[(which(thres.holds)[i]+1):(which(thres.holds)[i+1]),,drop=F]
	}
	if(nrow(my.domain)<2){
		all.domains[i] <- NA
		next
	}
	mean.meth <- apply(my.domain,2,mean)
	if(mean(mean.meth)>0.5){
		out.put <- mean.meth<(mean(mean.meth)-4*sqrt(sd(mean.meth)/length(mean.meth)))
	}else{
		out.put <- mean.meth>(mean(mean.meth)+4*sqrt(sd(mean.meth)/length(mean.meth)))
	}
	all.models <- PowerSet(nrow(my.domain))[-1]
	pred.dat <- data.frame(Outcome=out.put,t(my.domain))
	all.res <- list()
	dom.class[[i+1]] <- c(Positive=sum(out.put),Negative=sum(!out.put))
	for(j in 1:length(all.models)){
		sset <- all.models[[j]]
		form <- as.formula(paste("Outcome~",paste(row.names(my.domain)[as.vector(sset)+1],collapse=":")))
		log.model <- glm(form,pred.dat,family="binomial")
		mean.positive <- mean(apply(my.domain[as.vector(sset)+1,pred.dat$Outcome],1,mean))
		names(mean.positive) <- paste0("MeanPositive",names(mean.positive))
		sd.positive <- sd(apply(my.domain[as.vector(sset)+1,pred.dat$Outcome],1,mean),na.rm=T)
		names(sd.positive) <- paste0("SDPositive",names(sd.positive))
		mean.negative <- mean(apply(my.domain[as.vector(sset)+1,!pred.dat$Outcome],1,mean))
		names(mean.negative) <- paste0("MeanNegative",names(mean.negative))
		sd.negative <- sd(apply(my.domain[as.vector(sset)+1,!pred.dat$Outcome],1,mean),na.rm=T)
		names(sd.negative) <- paste0("SDNegative",names(sd.negative))
		all.res[[j]] <- c(Name=paste0(row.names(my.domain)[as.vector(sset)+1],collapse="&"),
				P.value=summary(log.model)$coefficient[2,"Pr(>|z|)"],
				mean.positive,
				sd.positive,
				mean.negative,
				sd.negative)
	}
	all.domains[[i+1]] <- all.res
}
names(all.domains) <- paste0("Domain",1:length(all.domains))
for(dom in 1:length(all.domains)){
	#for(pp in 1:length(all.domains[[dom]])){
		write.table(dom.class[[dom]],file=paste0("comethylation/upmeth_meth_",names(all.domains)[dom],".csv"),append=T)
		dat <- t(as.data.frame(all.domains[[dom]]))
		if(nrow(dat)>1){
			colnames(dat) <- c("Name",
			"P.value",
			"Mean.Mean.Positive",
			"SD.Mean.Positive",
			"Mean.Mean.Negative",
			"SD.Mean.Negativ")
			write.table(dat,file=paste0("comethylation/upmeth_meth_",names(all.domains)[dom],".csv"),append=T,row.names=F)
		}
	#}
}

library(eDMA)
meth.data <- read.table("MGMTe1methallsamples_wo.txt",sep="\t")
meth.data <- t(meth.data)
colnames(meth.data) <- meth.data[1,]
meth.data <- meth.data[-1,]
r.names <- meth.data[,1]
meth.data <- meth.data[,-1]
meth.data <- as.data.frame(apply(meth.data,2,as.numeric))
row.names(meth.data) <- r.names
mean.cpgs <- apply(meth.data,1,mean)
var.cpgs <- apply(meth.data,1,var)
dist.mean <- c()
dist.var <- c()
for(i in 1:length(mean.cpgs)-1){
	dist.mean <- c(dist.mean,abs(mean.cpgs[i]-mean.cpgs[i+1]))
	dist.var <- c(dist.var,abs(var.cpgs[i]-var.cpgs[i+1]))
}
thres.holds <- (dist.mean>quantile(dist.mean,.9))|(dist.var>quantile(dist.var,.9))
thres.holds[length(thres.holds)+1] <- TRUE

all.domains <- list()
dom.class <- list()
for(i in 0:(length(which(thres.holds))-1)){
	if(i==0){
		my.domain <- meth.data[0:which(thres.holds)[1],,drop=F]
	}else{
		my.domain <- meth.data[(which(thres.holds)[i]+1):(which(thres.holds)[i+1]),,drop=F]
	}
	if(nrow(my.domain)<2){
		all.domains[i] <- NA
		next
	}
	mean.meth <- apply(my.domain,2,mean)
	if(mean(mean.meth)>0.5){
		out.put <- mean.meth<(mean(mean.meth)-2*sqrt(sd(mean.meth)/length(mean.meth)))
	}else{
		out.put <- mean.meth>(mean(mean.meth)+2*sqrt(sd(mean.meth)/length(mean.meth)))
	}
	all.models <- PowerSet(nrow(my.domain))[-1]
	pred.dat <- data.frame(Outcome=out.put,t(my.domain))
	all.res <- list()
	dom.class[[i+1]] <- c(Positive=sum(out.put),Negative=sum(!out.put))
	for(j in 1:length(all.models)){
		sset <- all.models[[j]]
		form <- as.formula(paste("Outcome~",paste(row.names(my.domain)[as.vector(sset)+1],collapse=":")))
		log.model <- glm(form,pred.dat,family="binomial")
		mean.positive <- mean(apply(my.domain[as.vector(sset)+1,pred.dat$Outcome],1,mean))
		names(mean.positive) <- paste0("MeanPositive",names(mean.positive))
		sd.positive <- sd(apply(my.domain[as.vector(sset)+1,pred.dat$Outcome],1,mean),na.rm=T)
		names(sd.positive) <- paste0("SDPositive",names(sd.positive))
		mean.negative <- mean(apply(my.domain[as.vector(sset)+1,!pred.dat$Outcome],1,mean))
		names(mean.negative) <- paste0("MeanNegative",names(mean.negative))
		sd.negative <- sd(apply(my.domain[as.vector(sset)+1,!pred.dat$Outcome],1,mean),na.rm=T)
		names(sd.negative) <- paste0("SDNegative",names(sd.negative))
		all.res[[j]] <- c(Name=paste0(row.names(my.domain)[as.vector(sset)+1],collapse="&"),
				P.value=summary(log.model)$coefficient[2,"Pr(>|z|)"],
				mean.positive,
				sd.positive,
				mean.negative,
				sd.negative)
	}
	all.domains[[i+1]] <- all.res
}
names(all.domains) <- paste0("Domain",1:length(all.domains))
for(dom in 1:length(all.domains)){
	#for(pp in 1:length(all.domains[[dom]])){
		write.table(dom.class[[dom]],file=paste0("comethylation/exon1_meth_",names(all.domains)[dom],".csv"),append=T)
		dat <- t(as.data.frame(all.domains[[dom]]))
		if(nrow(dat)>1){
			colnames(dat) <- c("Name",
			"P.value",
			"Mean.Mean.Positive",
			"SD.Mean.Positive",
			"Mean.Mean.Negative",
			"SD.Mean.Negativ")
			write.table(dat,file=paste0("comethylation/exon1_meth_",names(all.domains)[dom],".csv"),append=T,row.names=F)
		}
	#}
}

library(eDMA)
meth.data <- read.table("MGMTi1allsamples_wo.txt",sep="\t")
meth.data <- t(meth.data)
colnames(meth.data) <- meth.data[1,]
meth.data <- meth.data[-1,]
r.names <- meth.data[,1]
meth.data <- meth.data[,-1]
meth.data <- as.data.frame(apply(meth.data,2,as.numeric))
row.names(meth.data) <- r.names
mean.cpgs <- apply(meth.data,1,mean)
var.cpgs <- apply(meth.data,1,var)
dist.mean <- c()
dist.var <- c()
for(i in 1:length(mean.cpgs)-1){
	dist.mean <- c(dist.mean,abs(mean.cpgs[i]-mean.cpgs[i+1]))
	dist.var <- c(dist.var,abs(var.cpgs[i]-var.cpgs[i+1]))
}
thres.holds <- (dist.mean>quantile(dist.mean,.9))|(dist.var>quantile(dist.var,.9))
thres.holds[length(thres.holds)+1] <- TRUE

all.domains <- list()
dom.class <- list()
for(i in 0:(length(which(thres.holds))-1)){
	if(i==0){
		my.domain <- meth.data[0:which(thres.holds)[1],,drop=F]
	}else{
		my.domain <- meth.data[(which(thres.holds)[i]+1):(which(thres.holds)[i+1]),,drop=F]
	}
	if(nrow(my.domain)<2){
		all.domains[i] <- NA
		next
	}
	mean.meth <- apply(my.domain,2,mean)
	if(mean(mean.meth)>0.5){
		out.put <- mean.meth<(mean(mean.meth)-2*sqrt(sd(mean.meth)/length(mean.meth)))
	}else{
		out.put <- mean.meth>(mean(mean.meth)+2*sqrt(sd(mean.meth)/length(mean.meth)))
	}
	all.models <- PowerSet(nrow(my.domain))[-1]
	pred.dat <- data.frame(Outcome=out.put,t(my.domain))
	all.res <- list()
	dom.class[[i+1]] <- c(Positive=sum(out.put),Negative=sum(!out.put))
	for(j in 1:length(all.models)){
		sset <- all.models[[j]]
		form <- as.formula(paste("Outcome~",paste(row.names(my.domain)[as.vector(sset)+1],collapse=":")))
		log.model <- glm(form,pred.dat,family="binomial")
		mean.positive <- mean(apply(my.domain[as.vector(sset)+1,pred.dat$Outcome],1,mean))
		names(mean.positive) <- paste0("MeanPositive",names(mean.positive))
		sd.positive <- sd(apply(my.domain[as.vector(sset)+1,pred.dat$Outcome],1,mean),na.rm=T)
		names(sd.positive) <- paste0("SDPositive",names(sd.positive))
		mean.negative <- mean(apply(my.domain[as.vector(sset)+1,!pred.dat$Outcome],1,mean))
		names(mean.negative) <- paste0("MeanNegative",names(mean.negative))
		sd.negative <- sd(apply(my.domain[as.vector(sset)+1,!pred.dat$Outcome],1,mean),na.rm=T)
		names(sd.negative) <- paste0("SDNegative",names(sd.negative))
		all.res[[j]] <- c(Name=paste0(row.names(my.domain)[as.vector(sset)+1],collapse="&"),
				P.value=summary(log.model)$coefficient[2,"Pr(>|z|)"],
				mean.positive,
				sd.positive,
				mean.negative,
				sd.negative)
	}
	all.domains[[i+1]] <- all.res
}
names(all.domains) <- paste0("Domain",1:length(all.domains))
for(dom in 1:length(all.domains)){
	#for(pp in 1:length(all.domains[[dom]])){
		write.table(dom.class[[dom]],file=paste0("comethylation/intron1_meth_",names(all.domains)[dom],".csv"),append=T)
		dat <- t(as.data.frame(all.domains[[dom]]))
		if(nrow(dat)>1){
			colnames(dat) <- c("Name",
			"P.value",
			"Mean.Mean.Positive",
			"SD.Mean.Positive",
			"Mean.Mean.Negative",
			"SD.Mean.Negativ")
			write.table(dat,file=paste0("comethylation/intron1_meth_",names(all.domains)[dom],".csv"),append=T,row.names=F)
		}
	#}
}


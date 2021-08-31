library(RnBeads)
library(survival)
rnb.set <- load.rnb.set("report2019-12-05/rnbSet_unnormalized/")
rem.samples <- rep(FALSE,length(samples(rnb.set)))
rem.samples[seq(2,length(samples(rnb.set)),by=2)] <- TRUE
rnb.set <- remove.samples(rnb.set,rem.samples)
gr <- readRDS("/DEEP_fhgfs/projects/mscherer/data/other/MGMT_Homburg/Alea/final_run/up_CpG_positions.RDS")
anno.rnb <- makeGRangesFromDataFrame(annotation(rnb.set))
op <- findOverlaps(gr[21],anno.rnb)
sel.meth <- meth(rnb.set)[subjectHits(op),]
ph <- pheno(rnb.set)
to.export <- data.frame(cg01341123=sel.meth,ph[,c("days_to_death","vital_status","age_at_index","gender","race")])
to.export$survival <- as.numeric(as.character(to.export$days_to_death))/7
to.export$status <- ifelse(to.export$vital_status%in%"Dead",1,0)
colnames(to.export)[4:6] <- c("Age","Sex","Ethnicity")
write.csv(to.export,"survival_data_TCGA_cg01341123.csv")

form <- as.formula(Surv(survival,status)~cg01341123+Sex+Age+Ethnicity,to.export)
mod <- coxph(form,to.export)

# Bady CpG cg12981137, our CpG cg01341123, others cg23998405 cg02022136 cg25946389
library(RnBeads)
library(survival)
rnb.set <- load.rnb.set("report2019-12-05/rnbSet_unnormalized/")
rem.samples <- rep(FALSE,length(samples(rnb.set)))
rem.samples[seq(2,length(samples(rnb.set)),by=2)] <- TRUE
rnb.set <- remove.samples(rnb.set,rem.samples)
c.bio <- read.table("idh_mutants.tsv",sep="\t",header=T)
rem.samples <- pheno(rnb.set)$submitter_id%in%c.bio$Patient.ID
rnb.set <- remove.samples(rnb.set,rem.samples)
ph <- pheno(rnb.set)
sel.meth <- meth(rnb.set,row.names=T)["cg25946389",]
to.export <- data.frame(cg25946389=sel.meth,ph[,c("days_to_death","vital_status","age_at_index","gender")])
to.export$survival <- as.numeric(as.character(to.export$days_to_death))/7
to.export$status <- ifelse(to.export$vital_status%in%"Dead",1,0)
colnames(to.export)[4:6] <- c("Age","Sex")
write.csv(to.export,"survival_data_TCGA_450k_Bady_cg25946389_woIDH.csv")

form <- as.formula(Surv(survival,status)~cg25946389+Sex+Age,to.export)
mod <- coxph(form,to.export)

# TCGA 27k cohort
library(RnBeads)
library(survival)
rnb.set <- load.rnb.set("rnb_report_20200815/cluster_run/import_RnBSet/")
c.bio <- read.table("idh_mutants.tsv",sep="\t",header=T)
rem.samples <- pheno(rnb.set)$case_submitter_id%in%c.bio$Patient.ID
rnb.set <- remove.samples(rnb.set,rem.samples)
ph <- pheno(rnb.set)
ph$age_at_index <- as.numeric(as.character(ph$age_at_index))
sel.meth <- meth(rnb.set,row.names=T)["cg12981137",]
to.export <- data.frame(cg12981137=sel.meth,ph[,c("days_to_death","vital_status","age_at_index","gender","race")])
to.export$survival <- as.numeric(as.character(to.export$days_to_death))/7
to.export$status <- ifelse(to.export$vital_status%in%"Dead",1,0)
colnames(to.export)[4:6] <- c("Age","Sex","Ethnicity")
write.csv(to.export,"survival_data_TCGA_27k_Bady_cg12981137_woIDH.csv")

form <- as.formula(Surv(survival,status)~cg12981137+Sex+Age+Ethnicity,to.export)
mod <- coxph(form,to.export)


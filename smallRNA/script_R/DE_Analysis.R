set.seed(1234)
HOME <- "/Users/sja517/Documents/projects/Ellen/"
DATA <- paste(HOME,"dataset/",sep="")
SCRIPT<-paste(HOME,"script_R/",sep="")
OUTPUT_T<-paste(HOME,"output_T/",sep="")
setwd(HOME)
library(DESeq)

## load small RNA-seq data
mirna.data<-read.table(paste(DATA,"counts_mtx.txt",sep=""), header=T)
mirna.data<-mirna.data[!duplicated(mirna.data$Geneid),]
rownames(mirna.data)<-mirna.data$Geneid; 
mirna.data$Geneid<-NULL
colData<-read.csv("dataset/colData.csv")
# ---------------------------------------------------------------------------------|
#    G0 vs G2 day3   (no replicates)                                               |    
# ---------------------------------------------------------------------------------|
# SWcolData<-read.csv("dataset/SWcolData.csv")

SW1<-mirna.data[,c(1,2)]
SWcolData1<-colData[c(1,2),]
rownames(SWcolData1)<-colnames(SW1)
SWcolData1$index<-NULL
cond<-factor(SWcolData1$condition)
SW1 <-SW1[rowSums(SW1 > 0) != 0, ]
# SW1 <- SW1+1  ## avoid zero fold change
cds<-newCountDataSet(SW1, cond)
cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds, method='blind',sharingMode='fit-only', fitType="local")
res<-nbinomTest(cds,"G0","G2")
res.ordered<-res[order(res$pval),]
res.ordered.pval05<-res.ordered[which(res.ordered$pval <0.05),]
write.csv(res.ordered.pval05,file=paste(OUTPUT_T,"G0_vs_G2_day3.pval05.csv",sep=""))
rm(res, res.ordered, res.ordered.pval05, SW1, cds, cond, SWcolData1)

# ---------------------------------------------------------------------------------|
#    G0 vs G2 day3   (replicates)                                               |    
# ---------------------------------------------------------------------------------|
SW2<-mirna.data[,c(3,4,13,14)]
SWcolData2<-colData[c(3,4,13,14),]
rownames(SWcolData2)<-colnames(SW2)
SWcolData2$index<-NULL
SW2<-SW2[rowSums(SW2 > 0) != 0, ]
cond<-factor(SWcolData2$condition);
cds<-newCountDataSet(SW2, cond)
cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds, method="per-condition", fitType="local", sharingMode="fit-only")
res<-nbinomTest(cds,"G0","G2")
res.ordered<-res[order(res$pval),]
res.ordered.pval05<-res.ordered[which(res.ordered$pval <0.05),]
write.csv(res.ordered.pval05,file=paste(OUTPUT_T,"G0_vs_G2_day7.pval05.per-condition.csv",sep=""))
rm(res, res.ordered, res.ordered.pval05, SW2, cds, cond, SWcolData2)

# ---------------------------------------------------------------------------------|
#    G0 vs G2 day3   (replicates)                                               |    
# ---------------------------------------------------------------------------------|
SW2<-mirna.data[,c(3,4,13,14)]
SWcolData2<-colData[c(3,4,13,14),]
rownames(SWcolData2)<-colnames(SW2)
SWcolData2$index<-NULL
SW2<-SW2[rowSums(SW2 > 0) != 0, ]
cond<-factor(SWcolData2$condition);
cds<-newCountDataSet(SW2, cond)
cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds, method="pooled", fitType="local", sharingMode="fit-only")
res<-nbinomTest(cds,"G0","G2")
res.ordered<-res[order(res$pval),]
res.ordered.pval05<-res.ordered[which(res.ordered$pval <0.05),]
write.csv(res.ordered.pval05,file=paste(OUTPUT_T,"G0_vs_G2_day7.pval05.pooled.csv",sep=""))
rm(res, res.ordered, res.ordered.pval05, SW2, cds, cond, SWcolData2)

# ---------------------------------------------------------------------------------|
#    B6 vs MRL normal day0   (no replicates)                                       |    
# ---------------------------------------------------------------------------------|
B6vsMRL_day0<-mirna.data[,c(5,11)]
B6vsMRL_day0_colData<-colData[c(5,11),]
rownames(B6vsMRL_day0_colData)<-colnames(B6vsMRL_day0)
B6vsMRL_day0_colData$index<-NULL
B6vsMRL_day0_cond<-factor(B6vsMRL_day0_colData$tissuetype)
B6vsMRL_day0 <-B6vsMRL_day0[rowSums(B6vsMRL_day0 > 0) != 0, ]
# B6vsMRL_day0 <- B6vsMRL_day0+1  ## avoid zero fold change
B6vsMRL_day0_cds<-newCountDataSet(B6vsMRL_day0, B6vsMRL_day0_cond)
B6vsMRL_day0_cds<-estimateSizeFactors(B6vsMRL_day0_cds)
B6vsMRL_day0_cds<-estimateDispersions(B6vsMRL_day0_cds, method='blind',fitType="local", sharingMode='fit-only')
B6vsMRL_day0_res<-nbinomTest(B6vsMRL_day0_cds,"B6","MRL")
B6vsMRL_day0_res.ordered<-B6vsMRL_day0_res[order(B6vsMRL_day0_res$pval),]
B6vsMRL_day0_res.ordered.pval05<-B6vsMRL_day0_res.ordered[which(B6vsMRL_day0_res.ordered$pval <0.05),]
write.csv(B6vsMRL_day0_res.ordered.pval05,file=paste(OUTPUT_T,"B6_vs_MRL_normal_day0.pval05.csv",sep=""))
rm(B6vsMRL_day0,B6vsMRL_day0_cds,B6vsMRL_day0_colData,B6vsMRL_day0_cond,B6vsMRL_day0_res,B6vsMRL_day0_res.ordered,B6vsMRL_day0_res.ordered.pval05)

# ---------------------------------------------------------------------------------|
#    B6 vs MRL normal day7   (no replicates)                                       |    
# ---------------------------------------------------------------------------------|
B6vsMRL_day7<-mirna.data[,c(6,12)]
B6vsMRL_day7_colData<-colData[c(6,12),]
rownames(B6vsMRL_day7_colData)<-colnames(B6vsMRL_day7)
B6vsMRL_day7_colData$index<-NULL
B6vsMRL_day7_cond<-factor(B6vsMRL_day7_colData$tissuetype)
B6vsMRL_day7 <-B6vsMRL_day7[rowSums(B6vsMRL_day7 > 0) != 0, ]
# B6vsMRL_day7 <- B6vsMRL_day7+1  ## avoid zero fold change
B6vsMRL_day7_cds<-newCountDataSet(B6vsMRL_day7, B6vsMRL_day7_cond)
B6vsMRL_day7_cds<-estimateSizeFactors(B6vsMRL_day7_cds)
B6vsMRL_day7_cds<-estimateDispersions(B6vsMRL_day7_cds, method='blind',fitType="local", sharingMode='fit-only')
B6vsMRL_day7_res<-nbinomTest(B6vsMRL_day7_cds,"B6","MRL")
B6vsMRL_day7_res.ordered<-B6vsMRL_day7_res[order(B6vsMRL_day7_res$pval),]
B6vsMRL_day7_res.ordered.pval05<-B6vsMRL_day7_res.ordered[which(B6vsMRL_day7_res.ordered$pval <0.05),]
write.csv(B6vsMRL_day7_res.ordered.pval05,file=paste(OUTPUT_T,"B6_vs_MRL_normal_day7.pval05.csv",sep=""))
rm(B6vsMRL_day7,B6vsMRL_day7_cds,B6vsMRL_day7_colData,B6vsMRL_day7_cond,B6vsMRL_day7_res,B6vsMRL_day7_res.ordered,B6vsMRL_day7_res.ordered.pval05)

# ---------------------------------------------------------------------------------|
#    B6 vs MRL HFD day0   (no replicates)                                          |    
# ---------------------------------------------------------------------------------|
B6vsMRL_day0<-mirna.data[,c(8,7)]
B6vsMRL_day0_colData<-colData[c(8,7),]
rownames(B6vsMRL_day0_colData)<-colnames(B6vsMRL_day0)
B6vsMRL_day0_colData$index<-NULL
B6vsMRL_day0_cond<-factor(B6vsMRL_day0_colData$tissuetype)
B6vsMRL_day0 <-B6vsMRL_day0[rowSums(B6vsMRL_day0 > 0) != 0, ]
# B6vsMRL_day0 <- B6vsMRL_day0+1  ## avoid zero fold change
B6vsMRL_day0_cds<-newCountDataSet(B6vsMRL_day0, B6vsMRL_day0_cond)
B6vsMRL_day0_cds<-estimateSizeFactors(B6vsMRL_day0_cds)
B6vsMRL_day0_cds<-estimateDispersions(B6vsMRL_day0_cds, method='blind',fitType="local", sharingMode='fit-only')
B6vsMRL_day0_res<-nbinomTest(B6vsMRL_day0_cds,"B6","MRL")
B6vsMRL_day0_res.ordered<-B6vsMRL_day0_res[order(B6vsMRL_day0_res$pval),]
B6vsMRL_day0_res.ordered.pval05<-B6vsMRL_day0_res.ordered[which(B6vsMRL_day0_res.ordered$pval <0.05),]
write.csv(B6vsMRL_day0_res.ordered.pval05,file=paste(OUTPUT_T,"B6_vs_MRL_HFD_day0.pval05.csv",sep=""))
rm(B6vsMRL_day0,B6vsMRL_day0_cds,B6vsMRL_day0_colData,B6vsMRL_day0_cond,B6vsMRL_day0_res,B6vsMRL_day0_res.ordered,B6vsMRL_day0_res.ordered.pval05)

# ---------------------------------------------------------------------------------|
#    B6 vs MRL HFD day7   (no replicates)                                          |    
# ---------------------------------------------------------------------------------|
B6vsMRL_day7<-mirna.data[,c(9,10)]
B6vsMRL_day7_colData<-colData[c(9,10),]
rownames(B6vsMRL_day7_colData)<-colnames(B6vsMRL_day7)
B6vsMRL_day7_colData$index<-NULL
B6vsMRL_day7_cond<-factor(B6vsMRL_day7_colData$tissuetype)
B6vsMRL_day7 <-B6vsMRL_day7[rowSums(B6vsMRL_day7 > 0) != 0, ]
# B6vsMRL_day7 <- B6vsMRL_day7+1  ## avoid zero fold change
B6vsMRL_day7_cds<-newCountDataSet(B6vsMRL_day7, B6vsMRL_day7_cond)
B6vsMRL_day7_cds<-estimateSizeFactors(B6vsMRL_day7_cds)
B6vsMRL_day7_cds<-estimateDispersions(B6vsMRL_day7_cds, method='blind',fitType="local", sharingMode='fit-only')
B6vsMRL_day7_res<-nbinomTest(B6vsMRL_day7_cds,"B6","MRL")
B6vsMRL_day7_res.ordered<-B6vsMRL_day7_res[order(B6vsMRL_day7_res$pval),]
B6vsMRL_day7_res.ordered.pval05<-B6vsMRL_day7_res.ordered[which(B6vsMRL_day7_res.ordered$pval <0.05),]
write.csv(B6vsMRL_day7_res.ordered.pval05,file=paste(OUTPUT_T,"B6_vs_MRL_HFD_day7.pval05.csv",sep=""))
rm(B6vsMRL_day7,B6vsMRL_day7_cds,B6vsMRL_day7_colData,B6vsMRL_day7_cond,B6vsMRL_day7_res,B6vsMRL_day7_res.ordered,B6vsMRL_day7_res.ordered.pval05)

# ---------------------------------------------------------------------------------|
#    B6 HFD (8) vs B6 normal day0 (5)   (no replicates)                            |    
# ---------------------------------------------------------------------------------|
B6_HFDvsNormal_day0<-mirna.data[,c(8,5)]
B6_HFDvsNormal_day0_colData<-colData[c(8,5),]
rownames(B6_HFDvsNormal_day0_colData)<-colnames(B6_HFDvsNormal_day0)
B6_HFDvsNormal_day0_colData$index<-NULL
B6_HFDvsNormal_day0_cond<-factor(B6_HFDvsNormal_day0_colData$condition)
B6_HFDvsNormal_day0 <-B6_HFDvsNormal_day0[rowSums(B6_HFDvsNormal_day0 > 0) != 0, ]
# B6vsMRL_day0 <- B6vsMRL_day0+1  ## avoid zero fold change
B6_HFDvsNormal_day0_cds<-newCountDataSet(B6_HFDvsNormal_day0, B6_HFDvsNormal_day0_cond)
B6_HFDvsNormal_day0_cds<-estimateSizeFactors(B6_HFDvsNormal_day0_cds)
B6_HFDvsNormal_day0_cds<-estimateDispersions(B6_HFDvsNormal_day0_cds, method='blind',fitType="local", sharingMode='fit-only')
B6_HFDvsNormal_day0_res<-nbinomTest(B6_HFDvsNormal_day0_cds,"HFD","normal")
B6_HFDvsNormal_day0_res.ordered<-B6_HFDvsNormal_day0_res[order(B6_HFDvsNormal_day0_res$pval),]
B6_HFDvsNormal_day0_res.ordered.pval05<-B6_HFDvsNormal_day0_res.ordered[which(B6_HFDvsNormal_day0_res.ordered$pval <0.05),]
write.csv(B6_HFDvsNormal_day0_res.ordered.pval05,file=paste(OUTPUT_T,"B6_HFDvsNormal_day0.pval05.csv",sep=""))
rm(B6_HFDvsNormal_day0,B6_HFDvsNormal_day0_cds,B6_HFDvsNormal_day0_colData,B6_HFDvsNormal_day0_cond,B6_HFDvsNormal_day0_res,B6_HFDvsNormal_day0_res.ordered,B6_HFDvsNormal_day0_res.ordered.pval05)

# ---------------------------------------------------------------------------------|
#    B6 HFD (7) vs B6 normal day7 (11)   (no replicates)                            |    
# ---------------------------------------------------------------------------------|
B6_HFDvsNormal_day7<-mirna.data[,c(7,11)]
B6_HFDvsNormal_day7_colData<-colData[c(7,11),]
rownames(B6_HFDvsNormal_day7_colData)<-colnames(B6_HFDvsNormal_day7)
B6_HFDvsNormal_day7_colData$index<-NULL
B6_HFDvsNormal_day7_cond<-factor(B6_HFDvsNormal_day7_colData$condition)
B6_HFDvsNormal_day7 <-B6_HFDvsNormal_day7[rowSums(B6_HFDvsNormal_day7 > 0) != 0, ]
# B6vsMRL_day7 <- B6vsMRL_day7+1  ## avoid zero fold change
B6_HFDvsNormal_day7_cds<-newCountDataSet(B6_HFDvsNormal_day7, B6_HFDvsNormal_day7_cond)
B6_HFDvsNormal_day7_cds<-estimateSizeFactors(B6_HFDvsNormal_day7_cds)
B6_HFDvsNormal_day7_cds<-estimateDispersions(B6_HFDvsNormal_day7_cds, method='blind',fitType="local", sharingMode='fit-only')
B6_HFDvsNormal_day7_res<-nbinomTest(B6_HFDvsNormal_day7_cds,"HFD","normal")
B6_HFDvsNormal_day7_res.ordered<-B6_HFDvsNormal_day7_res[order(B6_HFDvsNormal_day7_res$pval),]
B6_HFDvsNormal_day7_res.ordered.pval05<-B6_HFDvsNormal_day7_res.ordered[which(B6_HFDvsNormal_day7_res.ordered$pval <0.05),]
write.csv(B6_HFDvsNormal_day7_res.ordered.pval05,file=paste(OUTPUT_T,"B6_HFDvsNormal_day7.pval05.csv",sep=""))
rm(B6_HFDvsNormal_day7,B6_HFDvsNormal_day7_cds,B6_HFDvsNormal_day7_colData,B6_HFDvsNormal_day7_cond,B6_HFDvsNormal_day7_res,B6_HFDvsNormal_day7_res.ordered,B6_HFDvsNormal_day7_res.ordered.pval05)

# ---------------------------------------------------------------------------------|
#    MRL HFD (7) vs MRL normal day0 (11)   (no replicates)                            |    
# ---------------------------------------------------------------------------------|
MRL_HFDvsNormal_day0<-mirna.data[,c(7,11)]
MRL_HFDvsNormal_day0_colData<-colData[c(7,11),]
rownames(MRL_HFDvsNormal_day0_colData)<-colnames(MRL_HFDvsNormal_day0)
MRL_HFDvsNormal_day0_colData$index<-NULL
MRL_HFDvsNormal_day0_cond<-factor(MRL_HFDvsNormal_day0_colData$condition)
MRL_HFDvsNormal_day0 <-MRL_HFDvsNormal_day0[rowSums(MRL_HFDvsNormal_day0 > 0) != 0, ]
# B6vsMRL_day0 <- B6vsMRL_day0+1  ## avoid zero fold change
MRL_HFDvsNormal_day0_cds<-newCountDataSet(MRL_HFDvsNormal_day0, MRL_HFDvsNormal_day0_cond)
MRL_HFDvsNormal_day0_cds<-estimateSizeFactors(MRL_HFDvsNormal_day0_cds)
MRL_HFDvsNormal_day0_cds<-estimateDispersions(MRL_HFDvsNormal_day0_cds, method='blind',fitType="local", sharingMode='fit-only')
MRL_HFDvsNormal_day0_res<-nbinomTest(MRL_HFDvsNormal_day0_cds,"HFD","normal")
MRL_HFDvsNormal_day0_res.ordered<-MRL_HFDvsNormal_day0_res[order(MRL_HFDvsNormal_day0_res$pval),]
MRL_HFDvsNormal_day0_res.ordered.pval05<-MRL_HFDvsNormal_day0_res.ordered[which(MRL_HFDvsNormal_day0_res.ordered$pval <0.05),]
write.csv(MRL_HFDvsNormal_day0_res.ordered.pval05,file=paste(OUTPUT_T,"MRL_HFDvsNormal_day0.pval05.csv",sep=""))
rm(MRL_HFDvsNormal_day0,MRL_HFDvsNormal_day0_cds,MRL_HFDvsNormal_day0_colData,MRL_HFDvsNormal_day0_cond,MRL_HFDvsNormal_day0_res,MRL_HFDvsNormal_day0_res.ordered,MRL_HFDvsNormal_day0_res.ordered.pval05)

# ---------------------------------------------------------------------------------|
#    MRL HFD (9) vs MRL normal day7 (12)   (no replicates)                            |    
# ---------------------------------------------------------------------------------|
MRL_HFDvsNormal_day7<-mirna.data[,c(9,12)]
MRL_HFDvsNormal_day7_colData<-colData[c(9,12),]
rownames(MRL_HFDvsNormal_day7_colData)<-colnames(MRL_HFDvsNormal_day7)
MRL_HFDvsNormal_day7_colData$index<-NULL
MRL_HFDvsNormal_day7_cond<-factor(MRL_HFDvsNormal_day7_colData$condition)
MRL_HFDvsNormal_day7 <-MRL_HFDvsNormal_day7[rowSums(MRL_HFDvsNormal_day7 > 0) != 0, ]
# B6vsMRL_day7 <- B6vsMRL_day7+1  ## avoid zero fold change
MRL_HFDvsNormal_day7_cds<-newCountDataSet(MRL_HFDvsNormal_day7, MRL_HFDvsNormal_day7_cond)
MRL_HFDvsNormal_day7_cds<-estimateSizeFactors(MRL_HFDvsNormal_day7_cds)
MRL_HFDvsNormal_day7_cds<-estimateDispersions(MRL_HFDvsNormal_day7_cds, method='blind',fitType="local", sharingMode='fit-only')
MRL_HFDvsNormal_day7_res<-nbinomTest(MRL_HFDvsNormal_day7_cds,"HFD","normal")
MRL_HFDvsNormal_day7_res.ordered<-MRL_HFDvsNormal_day7_res[order(MRL_HFDvsNormal_day7_res$pval),]
MRL_HFDvsNormal_day7_res.ordered.pval05<-MRL_HFDvsNormal_day7_res.ordered[which(MRL_HFDvsNormal_day7_res.ordered$pval <0.05),]
write.csv(MRL_HFDvsNormal_day7_res.ordered.pval05,file=paste(OUTPUT_T,"MRL_HFDvsNormal_day7.pval05.csv",sep=""))
rm(MRL_HFDvsNormal_day7,MRL_HFDvsNormal_day7_cds,MRL_HFDvsNormal_day7_colData,MRL_HFDvsNormal_day7_cond,MRL_HFDvsNormal_day7_res,MRL_HFDvsNormal_day7_res.ordered,MRL_HFDvsNormal_day7_res.ordered.pval05)

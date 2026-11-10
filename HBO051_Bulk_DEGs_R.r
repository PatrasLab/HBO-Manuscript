
library( "DESeq2" )
library("dplyr") 
library("tibble")

####Pre-processing####

#import counts
new_hbo_data <- read.csv(file = "gene_count.csv", header = T, row.names = 1)

#import metadata
hbo_coldata<-read.csv("HBO_metadata.csv", header=T)

#make dataframes
countdata <- as.data.frame(new_hbo_data)
metadata <- as.data.frame(hbo_coldata)

#store gene names for later
genes <- countdata$gene_name
countdata$gene_name <- NULL

#Checking names are the same
names(countdata)
metadata$ID

#DEseq
ddsM_new_HBO <- DESeqDataSetFromMatrix(countData = countdata,
                                        colData = metadata,
                                        design = ~ Treatment_Full)

ddsDEM_new_HBO <- DESeq(ddsM_new_HBO, minReplicatesForReplace=Inf)

#Normalize Counts
HBO_new_normCounts <- counts(ddsDEM_new_HBO, normalized = T)

#replace IDs with readable gene names
rownames(ddsDEM_new_HBO) <- genes
rownames(HBO_new_normCounts) <- genes

write.csv(HBO_new_normCounts, "HBO_new_Normalized.counts.csv")

####DEG analysis####

#E. coli vs Urine
res_EC_4h_vs_Urine_4h <- results(ddsDEM_new_HBO, contrast=c("Treatment_Full", "H51_EC_4","H51_Urine_4"), alpha = 0.05, cooksCutoff=FALSE, lfcThreshold = 1)
write.csv(as.data.frame(res_EC_4h_vs_Urine_4h), file="res_new_HBO051_EC_4h_vs_Urine_4h.csv")

sigsub_EC_4h_vs_Urine_4h<- subset(res_EC_4h_vs_Urine_4h, padj < 0.05)
write.csv(as.data.frame(sigsub_EC_4h_vs_Urine_4h), file="sigsub_new_HBO051_EC_4h_vs_Urine_4h.csv")
res_EC_4h_vs_Urine_4h_sig <-read.csv("res_new_HBO051_EC_4h_vs_Urine_4h.csv", header = T)


#4h HP3 vs Urine
res_HP3_4h_vs_Urine_4h <- results(ddsDEM_new_HBO, contrast=c("Treatment_Full", "H51_HP3_4","H51_Urine_4"), alpha = 0.05, cooksCutoff=FALSE, lfcThreshold = 1)
write.csv(as.data.frame(res_new_HP3_4h_vs_Urine_4h), file="res_new_HBO051_HP3_4h_vs_Urine_4h.csv")

sigsub_HP3_4h_vs_Urine_4h<- subset(res_HP3_4h_vs_Urine_4h, padj < 0.05)
write.csv(as.data.frame(sigsub_HP3_4h_vs_Urine_4h), file="sigsub_new_HBO051_HP3_4h_vs_Urine_4h.csv")
res_HP3_4h_vs_Urine_4h_sig <-read.csv("res_new_HBO051_HP3_4h_vs_Urine_4h.csv", header = T)


#24h E. coli vs Urine
res_EC_24h_vs_Urine_24h <- results(ddsDEM_new_HBO, contrast=c("Treatment_Full", "H51_EC_24","H51_Urine_24"), alpha = 0.05, cooksCutoff=FALSE, lfcThreshold = 1)
write.csv(as.data.frame(res_EC_24h_vs_Urine_24h), file="res_new_HBO051_EC_24h_vs_Urine_24h.csv")

sigsub_EC_24h_vs_Urine_24h<- subset(res_EC_24h_vs_Urine_24h, padj < 0.05)
write.csv(as.data.frame(sigsub_EC_24h_vs_Urine_24h), file="sigsub_new_HBO051_EC_24h_vs_Urine_24h.csv")
res_EC_24h_vs_Urine_24h_sig <-read.csv("res_new_HBO051_EC_24h_vs_Urine_24h.csv", header = T)


#24h HP3 vs Urine
res_HP3_24h_vs_Urine_24h <- results(ddsDEM_new_HBO, contrast=c("Treatment_Full", "H51_HP3_24","H51_Urine_24"), alpha = 0.05, cooksCutoff=FALSE, lfcThreshold = 1)
write.csv(as.data.frame(res_HP3_24h_vs_Urine_24h), file="res_new_HBO051_HP3_24h_vs_Urine_24h.csv")

sigsub_HP3_24h_vs_Urine_24h<- subset(res_HP3_24h_vs_Urine_24h, padj < 0.05)
write.csv(as.data.frame(sigsub_HP3_24h_vs_Urine_24h), file="sigsub_new_HBO051_HP3_24h_vs_Urine_24h.csv")
res_HP3_24h_vs_Urine_24h_sig <-read.csv("res_new_HBO051_HP3_24h_vs_Urine_24h.csv", header = T)



##EC+HP3 vs EC Only

#4h EC+HP3 vs 4h EC Only
res_ECHP3_4h_vs_EC_4h <- results(ddsDEM_new_HBO, contrast=c("Treatment_Full", "H51_ECHP3_4","H51_EC_4"), alpha = 0.05, cooksCutoff=FALSE, lfcThreshold = 1)
write.csv(as.data.frame(res_ECHP3_4h_vs_EC_4h), file="res_new_HBO051_ECHP3_4h_vs_HBO051_EC_4h.csv")

sigsub_res_ECHP3_4h_vs_EC_4h<- subset(res_ECHP3_4h_vs_EC_4h, padj < 0.05)
write.csv(as.data.frame(sigsub_res_new_ECHP3_4h_vs_EC_4h), file="sigsub_res_new_HBO051_ECHP3_4h_vs_HBO051_EC_4h.csv")
res_new_ECHP3_4h_vs_EC_4h_sig <-read.csv("res_new_HBO051_ECHP3_4h_vs_HBO051_EC_4h.csv", header = T)


#24h EC+HP3 vs 24h EC Only
res_ECHP3_24h_vs_EC_24h <- results(ddsDEM_new_HBO, contrast=c("Treatment_Full", "H51_ECHP3_24","H51_EC_24"), alpha = 0.05, cooksCutoff=FALSE, lfcThreshold = 1)
write.csv(as.data.frame(res_ECHP3_24h_vs_EC_24h), file="res_new_HBO051_ECHP3_24h_vs_HBO051_EC_24h.csv")

sigsub_res_ECHP3_24h_vs_EC_24h<- subset(res_ECHP3_24h_vs_EC_24h, padj < 0.05)
write.csv(as.data.frame(sigsub_res_ECHP3_24h_vs_EC_24h), file="sigsub_res_new_HBO051_ECHP3_24h_vs_HBO051_EC_24h.csv")
res_ECHP3_24h_vs_EC_24h_sig <-read.csv("res_new_HBO051_ECHP3_24h_vs_HBO051_EC_24h.csv", header = T)


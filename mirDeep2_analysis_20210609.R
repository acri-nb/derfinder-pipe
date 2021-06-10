
#Dependencies
set.seed(234)
library("rtracklayer")
library("GenomicRanges")
library("systemsbio")
library("derfinder")
library("ggplot2")
library("polyester")
library(BSgenome.Hsapiens.UCSC.hg19)
library("caret")
library("VennDiagram")
library("RColorBrewer")
library("org.Hs.eg.db")
library("clusterProfiler")
library("reshape2")
library("data.table")
library(srnadiff)
library("gprofiler2")


#plot style
Eric_theme<-theme(panel.grid.major.x=element_blank(),
                  panel.grid.minor.x=element_blank(),
                  axis.title = element_text(size = 24),
                  axis.text = element_text(size = 18),
                  panel.background = element_rect(fill = "white"),
                  panel.border = element_blank(),
                  axis.line = element_line(colour = "black"),
                  panel.grid.major.y = element_blank(),
                  legend.title = element_text(size = 18))


#Assign wdir, filenames and GTF annotation
setwd("D:/ALS_project")
files <- list.files(pattern = "\\.bam$") 
files <- files[1:14]

gffdata<-import.gff("D:/ALS_project/Analysis20200109/gencode.v19.annotation.gtf")
gffdata$score<-NULL
seqlevelsStyle(gffdata) <- "UCSC"

# We used a single feature-type annotation to simplify benchmarking, in this case gene
gffdata<- gffdata[gffdata$type == "gene",]

#Keep a list of contigs and genes to exclude from downstream analyses (useful later, X, Y and M)
gff_sub <- gffdata[seqnames(gffdata) %in% c("chrX", "chrY", "chrM"),]
supp_genes <- unique(gff_sub$gene_id)


#Standardize DGE analysis
DGE_EdgeR<-function(mat){
  library("edgeR")
  
  design<- matrix(ncol=2,nrow=14)
  rownames(design)<-colnames(mat)
  colnames(design)<-c("group","dummy")
  design[,1]<-c(rep(0,3),rep(1,8), rep(0,3))
  
  pheno <- new("AnnotatedDataFrame", data = as.data.frame(design))
  factor1<-as.factor(pheno$group)
  dge <- DGEList(mat)
  dge <- calcNormFactors(dge, method = "TMM")
  pseudo_TMM <- log2(scale(dge$counts,center=FALSE,scale=dge$samples$norm.factors)+1)
  
  
  design_matrix <- model.matrix(~ factor1)
  
  
  #Estimate dispertions, a necessary step to assess gene-level scatter in our data.
  dge <- estimateGLMCommonDisp(dge, design_matrix)
  dge <- estimateGLMTrendedDisp(dge, design_matrix)
  dge <- estimateGLMTagwiseDisp(dge, design_matrix)
  
  # par(mfrow=c(1,1))
  # plotBCV(dge, main = paste0("BCV plot"))
  
  fit<-glmFit(dge,design_matrix)
  res <- glmLRT(fit, coef = 2)
  
  pvals<- data.frame(raw.pvalue = res$table$PValue,
                     FDR=p.adjust(res$table$PValue, "BH"))
  # hist(pvals$raw.pvalue, xlim=c(0,1))
  
  res<-topTags(res,nrow(dge$counts))
  alpha=0.05
  
  #Quick check to see how many genes are flagges as DE and if there are any imbalances.
  return(res$table)
}


#Prepare DGE results by tool

#FeatureCounts
setwd("D:/ALS_project/Analysis20200630")
FC_counts<- read.table("FC_Bowtie_gencode19_2020.txt", header = T, row.names = 1)
FC_counts<- FC_counts[,-c(1,2,3,4,5)]
FC_counts<- FC_counts[!(rownames(FC_counts) %in% supp_genes),]
FC_counts<- FC_counts[rowMeans(FC_counts) > 20,]
FC_DGE<- DGE_EdgeR(FC_counts)

#HTSeq
HT_counts<- read.table("HT_gencode19_2020.txt", header = T, row.names = 1)
HT_counts<- HT_counts[!(rownames(HT_counts) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")),]
HT_counts<- HT_counts[!(rownames(HT_counts) %in% supp_genes),]
HT_counts<- HT_counts[rowMeans(HT_counts) > 20,]
HT_DGE<- DGE_EdgeR(HT_counts)

#Derfinder
Der_counts<- read.table("Derfinder_nocutoff_nonorm_20200707.txt", header = T, sep= "\t")
rownames(Der_counts)<- Der_counts$gene_id
Der_counts<- Der_counts[26:39]
Der_counts<- Der_counts[rowMeans(Der_counts) > 20,]
Der_counts<- Der_counts[,names(FC_counts)]
Der_DGE<- DGE_EdgeR(Der_counts)


#srnadiff

setwd("D:/ALS_project/Analysis20200109")
SampleInfo<- matrix(nrow = 14, ncol = 3)
SampleInfo[,1]<-files
SampleInfo[,2]<-files
SampleInfo[,3]<-c(rep(0,3),rep(1,8), rep(0,3))
colnames(SampleInfo)<-c("FileName", "SampleName", "Condition")
SampleInfo <- as.data.frame(SampleInfo)
annotReg   <- readAnnotation("gencode.v19.annotation.gtf", feature="gene", source=NULL)
setwd("D:/ALS_project")
srnaExp    <- srnadiffExp(files, SampleInfo, annotReg)
parameters(srnaExp) <- list(minDepth=0, minSize=10, minGap = 0, minOverlap = 1)

srnaExp    <- srnadiff(srnaExp,segMethod = "annotation")


CM<- srnaExp@countMatrix
CM<- CM[order(rownames(CM)),]
trs<- srnaExp@annotReg[names(srnaExp@annotReg) %in% rownames(CM),]
trs<-cbind(names(trs),trs$gene_id)
colnames(trs) <- c("Internal_ids", "ENSG_ids")
trs <- data.frame(trs)
trs<- trs[order(trs$Internal_ids),]
rownames(CM)<-trs$ENSG_ids
CM<- CM[!(rownames(CM) %in% supp_genes),]
srna_counts<- CM[rowMeans(CM) > 20,]
srna_counts <- srna_counts[,names(FC_counts)]
CM <- CM[,names(FC_counts)]
srna_DGE<- DGE_EdgeR(srna_counts)


##shortstack data##

ss<- read.table(file.choose(),header = T, sep = "\t", row.names = 1)

ssplit1<- colsplit(rownames(ss), ":", names=c("chr","start"))
ssplit2<- colsplit(ssplit1$start, "-", names=c("start","end"))
ss <- cbind(ssplit1$chr, ssplit2, ss[,3:16])
col.order <- c("ssplit1$chr", "start", "end", "ACRI40", "ACRI42", "ACRI62", "ALS16", "ALS1", "ALS2", "ALS3_cutadapt", "ALS4", "ALS5", "ALS6", "ALS7", "Control_orange", "Control_purple", "Control_red")
ss<- ss[,col.order]

#find overlapping regions with annotation
ss_gr<- makeGRangesFromDataFrame(ss, keep.extra.columns = T, seqnames.field = "ssplit1$chr")
ss_gr<- ss_gr[!(seqnames(ss_gr) %in% c("chrX", "chrY", "chrM")),]
ss_ol<- subsetByOverlaps.keepAllMeta(gffdata, ss_gr)
mt_dfr<-cbind(as.character(seqnames(ss_ol)),start(ss_ol),end(ss_ol))
colnames(mt_dfr)<-c("contig","start","end")
mt_dfr<-cbind(mt_dfr,values(ss_ol))
print(head(mt_dfr))

#Get sum of counts per gene / meta-feature
summing<-function(x) sum(na.omit(sapply(strsplit(x,";"),as.integer))) 
for (i in names(ss)[3:17]){
  sums<-sapply(mt_dfr[[i]], summing)
  #replace with new values
  mt_dfr[[i]]<-sums
}
mt_dfr$contig<-as.character(mt_dfr$contig)
mt_dfr$start<-as.character(mt_dfr$start)
mt_dfr$end<-as.character(mt_dfr$end)

#DGE analysis for ShortStack
sstack<- mt_dfr
rownames(sstack) <- sstack$gene_id
sstack_counts <- sstack[,25:38]
sstack_counts <- as.matrix(sstack_counts)
sstack_counts <- sstack_counts[rowMeans(sstack_counts) > 20,]
sstack_counts <- as.data.frame(sstack_counts)
ss_DGE <- DGE_EdgeR(sstack_counts)






# ############### Comparison with mirDeep2 ######################
# Read in files... expression results from mirDeep2 were used with mature clipped miRNA from human and chimp. 
mirdeep <- read.table(file.choose(), header = T, stringsAsFactors = F, sep = "\t")
mirdeep <- aggregate(.~mirdeep$ENSG, mirdeep,FUN=max)
rownames(mirdeep) <- mirdeep$ENSG

#Remove unecessary columns
mirdeep <- mirdeep[,-c(1,2)]
mirdeep <- as.matrix(mirdeep)
mirdeep1 <- apply(mirdeep, 2, as.numeric)
rownames(mirdeep1) <- rownames(mirdeep)
mirdeep1<- mirdeep1[rowSums(mirdeep1)>0,]
colnames(mirdeep1)<-c("Sample1", "Sample2", "Sample3") # these were samples ALS1, ALS2 and ALS4


# Remove trailing version suffix from ENSMBL annotation and get only mirs common to mirDeep2 results
names<- colsplit(rownames(FC_counts), "[.]", names=c("col1","col2"))
rownames(FC_counts)<- names$col1
FC_mir <- FC_counts[rownames(FC_counts) %in% rownames(mirdeep1),]
names<- colsplit(rownames(HT_counts), "[.]", names=c("col1","col2"))
rownames(HT_counts)<- names$col1
HT_mir <- HT_counts[rownames(HT_counts) %in% rownames(mirdeep1),]
names<- colsplit(rownames(Der_counts), "[.]", names=c("col1","col2"))
rownames(Der_counts)<- names$col1
Der_mir <- Der_counts[rownames(Der_counts) %in% rownames(mirdeep1),]
names<- colsplit(rownames(sstack_counts), "[.]", names=c("col1","col2"))
rownames(sstack_counts) <- names$col1
ss_mir <- sstack_counts[rownames(sstack_counts) %in% rownames(mirdeep1),]

names<- colsplit(rownames(CM), "[.]", names=c("col1","col2"))
rownames(CM)<- names$col1
srna_mir <- CM[rownames(CM) %in% rownames(mirdeep1),]
srna_mir <- as.data.frame(srna_mir)
srna_mir[,2:4]<-lapply(srna_mir[,2:4], as.integer)

#Quick check of correlation between mirdeep results and other tools
#Only compare to mirdeep for 3 samples, this was to aleviate time, as mirdeep2 is run as command-line tool and takes some time
FC_mir <- cbind(rownames(FC_mir),FC_mir[,c(5,6,8)])
HT_mir <- cbind(rownames(HT_mir),HT_mir[,c(5,6,8)])
Der_mir <- cbind(rownames(Der_mir),Der_mir[,c(5,6,8)])
srna_mir <- cbind(rownames(srna_mir),srna_mir[,c(5,6,8)])
ss_mir <- cbind(rownames(ss_mir),ss_mir[,c(5,6,8)])

names(FC_mir)<-c("ENSG" ,"Sample1", "Sample2", "Sample3")
names(HT_mir)<-c("ENSG" ,"Sample1", "Sample2", "Sample3")
names(Der_mir)<-c("ENSG" ,"Sample1", "Sample2", "Sample3")
names(srna_mir)<-c("ENSG" ,"Sample1", "Sample2", "Sample3")
names(ss_mir)<-c("ENSG" ,"Sample1", "Sample2", "Sample3")

#join rownames with data from each tool
combi <- matrix(ncol = 4, nrow = 454)
rownames(combi) <- rownames(mirdeep1)    
#need full matrices for comparisons with mirs
FC_full <- matrix(ncol=4, nrow=454)
rownames(FC_full) <- rownames(mirdeep1)
HT_full <- matrix(ncol=4, nrow=454)
rownames(HT_full) <- rownames(mirdeep1)
Der_full <- matrix(ncol=4, nrow=454)
rownames(Der_full) <- rownames(mirdeep1)
srna_full <- matrix(ncol=4, nrow=454)
rownames(srna_full) <- rownames(mirdeep1)
ss_full <- matrix(ncol=4, nrow=454)
rownames(ss_full) <- rownames(mirdeep1)
for (i in seq(1:4)){
  FC_full[,i] <- FC_mir[,i][match(rownames(FC_full), rownames(FC_mir))]
  HT_full[,i] <- HT_mir[,i][match(rownames(HT_full), rownames(HT_mir))]
  Der_full[,i] <- Der_mir[,i][match(rownames(Der_full), rownames(Der_mir))]
  srna_full[,i] <- srna_mir[,i][match(rownames(srna_full), rownames(srna_mir))]
  ss_full[,i] <- ss_mir[,i][match(rownames(ss_full), rownames(ss_mir))]
}  

#format the dataframes for plotting
FC_full <- FC_full[,-1]
HT_full <- HT_full[,-1]
Der_full <- Der_full[,-1]
srna_full <- srna_full[,-1]
ss_full <- ss_full[,-1]
class(FC_full) <- "numeric"
class(HT_full) <- "numeric"
class(Der_full) <- "numeric"
class(srna_full) <- "numeric"
class(ss_full) <- "numeric"
FC_full <- as.data.frame(FC_full)
HT_full <- as.data.frame(HT_full)
Der_full <- as.data.frame(Der_full)
srna_full <- as.data.frame(srna_full)
ss_full <- as.data.frame(ss_full)
names(FC_full) <- c("Sample1", "Sample2", "Sample3")
names(HT_full) <- c("Sample1", "Sample2", "Sample3")
names(Der_full) <- c("Sample1", "Sample2", "Sample3")
names(srna_full) <- c("Sample1", "Sample2", "Sample3")
names(ss_full) <- c("Sample1", "Sample2", "Sample3")
FC_full$ID <- rownames(FC_full)
HT_full$ID <- rownames(HT_full)
Der_full$ID <- rownames(Der_full)
srna_full$ID <- rownames(srna_full)
ss_full$ID <- rownames(ss_full)


#Get data in long format for ggplot
FC_mir<- melt(FC_full)
HT_mir<- melt(HT_full)
Der_mir<- melt(Der_full)
srna_mir<- melt(srna_full)
ss_mir <- melt(ss_full)
mirdeep1<- melt(mirdeep1)


FC_mir$tool<- rep("featureCounts", 1362)
HT_mir$tool<- rep("HTSeq", 1362)
Der_mir$tool<- rep("derfinder", 1362)
srna_mir$tool<- rep("srnadiff", 1362)
ss_mir$tool <- rep("shortstack", 1362)
mirdeep1$tool<- rep("mirdeep2", 1362)


names(FC_mir)<-c("ID", "Sample", "CountsFC", "tool")
names(HT_mir)<-c("ID", "Sample", "CountsHT", "tool")
names(Der_mir)<-c("ID", "Sample", "CountsDER", "tool")
names(srna_mir)<-c("ID", "Sample", "CountsSRNA", "tool")
names(ss_mir)<-c("ID", "Sample", "CountsSS", "tool")
names(mirdeep1)<-c("ID", "Sample", "CountsMIR", "tool")

#plot and correlations
gd<- cbind(mirdeep1,FC_mir$CountsFC)
g<-ggplot(gd, aes(x=log(CountsMIR), y=log(FC_mir$CountsFC), colour = Sample))+geom_point(alpha = .5)+xlab("mirDeep2 (log-counts)")+ylab("featureCounts (log-counts)")+Eric_theme
#save 500x500


cor.test(mirdeep1$CountsMIR, FC_mir$CountsFC, method = "spearman")









######## COPARISON WITH DIRECT ALIGNMENT TO ANNOTATION#########
## checking against old mir data from Saucier et al. 

mirdeep <- read.table(file.choose(), header = T, stringsAsFactors = F, sep = "\t")
mirdeep <- aggregate(.~mirdeep$ENSG, mirdeep,FUN=max)
rownames(mirdeep) <- mirdeep$ENSG
mirdeep <- mirdeep[,-c(1,2)]
mirdeep <- as.matrix(mirdeep)
mirdeep1 <- apply(mirdeep, 2, as.numeric)
rownames(mirdeep1) <- rownames(mirdeep)
mirdeep1<- mirdeep1[rowSums(mirdeep1)>0,]
colnames(mirdeep1)<-c("Sample1", "Sample2", "Sample3")

# Remove trailing version suffix from ENSMBL annotation and get only mirs common to mirDeep2
names<- colsplit(rownames(FC_counts), "[.]", names=c("col1","col2"))
rownames(FC_counts)<- names$col1
FC_mir <- FC_counts[rownames(FC_counts) %in% rownames(mirdeep1),]
names<- colsplit(rownames(HT_counts), "[.]", names=c("col1","col2"))
rownames(HT_counts)<- names$col1
HT_mir <- HT_counts[rownames(HT_counts) %in% rownames(mirdeep1),]
names<- colsplit(rownames(Der_counts), "[.]", names=c("col1","col2"))
rownames(Der_counts)<- names$col1
Der_mir <- Der_counts[rownames(Der_counts) %in% rownames(mirdeep1),]
names<- colsplit(rownames(sstack_counts), "[.]", names=c("col1","col2"))
rownames(sstack_counts) <- names$col1
ss_mir <- sstack_counts[rownames(sstack_counts) %in% rownames(mirdeep1),]

names<- colsplit(rownames(CM), "[.]", names=c("col1","col2"))
rownames(CM)<- names$col1
srna_mir <- CM[rownames(CM) %in% rownames(mirdeep1),]
srna_mir <- as.data.frame(srna_mir)
srna_mir[,2:4]<-lapply(srna_mir[,2:4], as.integer)


#convert IDs
gc1<- gconvert(query=rownames(FC_mir), organism="hsapiens", target="MIRBASE", mthreshold=1, filter_na=T)
gc2<- gconvert(query=rownames(HT_mir), organism="hsapiens", target="MIRBASE", mthreshold=1, filter_na=T)
gc3<- gconvert(query=rownames(Der_mir), organism="hsapiens", target="MIRBASE", mthreshold=1, filter_na=T)
gc4<- gconvert(query=rownames(srna_mir), organism="hsapiens", target="MIRBASE", mthreshold=1, filter_na=T)
gc5<- gconvert(query=rownames(ss_mir), organism="hsapiens", target="MIRBASE", mthreshold=1, filter_na=T)

FC_mir$name <- gc1$target
HT_mir$name <- gc2$target
Der_mir$name <- gc3$target
srna_mir$name <- gc4$target
ss_mir$name <-gc5$target

newnames<- gsub(pattern = "-1$", replacement = "", x=FC_mir$name)
newnames<- gsub(pattern = "-2$", replacement = "", x=newnames)
newnames<- gsub(pattern = "-3$", replacement = "", x=newnames)
FC_mir$name <- newnames
FCM<- aggregate(as.matrix(FC_mir[,1:14])~FC_mir$name, data=FC_mir,FUN=sum, na.rm=T)
rownames(FCM)<- FCM$`FC_mir$name`
FCM <- FCM[,-1]

newnames<- gsub(pattern = "-1$", replacement = "", x=HT_mir$name)
newnames<- gsub(pattern = "-2$", replacement = "", x=newnames)
newnames<- gsub(pattern = "-3$", replacement = "", x=newnames)
HT_mir$name <- newnames
HTM<- aggregate(as.matrix(HT_mir[,1:14])~HT_mir$name, data=HT_mir,FUN=sum, na.rm=T)
rownames(HTM)<- HTM$`HT_mir$name`
HTM <- HTM[,-1]

newnames<- gsub(pattern = "-1$", replacement = "", x=Der_mir$name)
newnames<- gsub(pattern = "-2$", replacement = "", x=newnames)
newnames<- gsub(pattern = "-3$", replacement = "", x=newnames)
Der_mir$name <- newnames
DerM<- aggregate(as.matrix(Der_mir[,1:14])~Der_mir$name, data=Der_mir,FUN=sum, na.rm=T)
rownames(DerM)<- DerM$`Der_mir$name`
DerM <- DerM[,-1]

newnames<- gsub(pattern = "-1$", replacement = "", x=srna_mir$name)
newnames<- gsub(pattern = "-2$", replacement = "", x=newnames)
newnames<- gsub(pattern = "-3$", replacement = "", x=newnames)
srna_mir$name <- newnames
srnaM<- aggregate(as.matrix(srna_mir[,1:14])~srna_mir$name, data=srna_mir,FUN=sum, na.rm=T)
rownames(srnaM)<- srnaM$`srna_mir$name`
srnaM <- srnaM[,-1]

newnames<- gsub(pattern = "-1$", replacement = "", x=ss_mir$name)
newnames<- gsub(pattern = "-2$", replacement = "", x=newnames)
newnames<- gsub(pattern = "-3$", replacement = "", x=newnames)
ss_mir$name <- newnames
ssM<- aggregate(as.matrix(ss_mir[,1:14])~ss_mir$name, data=ss_mir,FUN=sum, na.rm=T)
rownames(ssM)<- ssM$`ss_mir$name`
ssM <- ssM[,-1]



##Run differential expr analysis
FC_DME<- DGE_EdgeR(FCM)
HT_DME<- DGE_EdgeR(HTM)
Der_DME<- DGE_EdgeR(DerM)
srna_DME<- DGE_EdgeR(srnaM)
ss_DME<- DGE_EdgeR(ssM)


##read-in results for all 14 samples with mirbase alignment 

mirbase <- read.table(file.choose(), header = T, sep = "\t")
newnames<- gsub(pattern = "-1$", replacement = "", x=mirbase$X)
newnames<- gsub(pattern = "-2$", replacement = "", x=newnames)
mirbase$X <- newnames
mirbase1 <- aggregate(cbind(mirbase$Control_ACRI_40,mirbase$Control_ACRI_42, mirbase$Control_ACRI_62, mirbase$ALS_1_visit1, mirbase$ALS_16_visit1, mirbase$ALS_2_visit1, mirbase$ALS_3_visit1, mirbase$ALS_4_visit1, mirbase$ALS_5_visit1, mirbase$ALS_6_visit1, mirbase$ALS_7_visit1, mirbase$Control_plasma_orange, mirbase$Control_plasma_purple, mirbase$Control_plasma_red)~mirbase$X, data=mirbase,FUN=sum, na.rm=T)
rownames(mirbase1) <- mirbase1$`mirbase$X`
mirbase1 <- mirbase1[,-1]
names(mirbase1)<- names(FC_mir)[1:14]
mirbase_DME <- DGE_EdgeR(mirbase1)


FC_DME <- FC_DME[order(rownames(FC_DME)),]
HT_DME <- HT_DME[order(rownames(HT_DME)),]
Der_DME <- Der_DME[order(rownames(Der_DME)),]
srna_DME <- srna_DME[order(rownames(srna_DME)),]
ss_DMA <- ss_DME[order(rownames(ss_DME)),]
mirbase_DME <- mirbase_DME[order(rownames(mirbase_DME)),]




# restrict analysis to detected mirs
all_DE_genes<- unique(c(rownames(FC_DME),rownames(HT_DME),rownames(ss_DME), rownames(srna_DME), rownames(Der_DME) ))
mirbase2<- matrix(ncol=14, nrow=250)
rownames(mirbase2) <- all_DE_genes
FC_full <- matrix(ncol=14, nrow=250)
#FC_full <- matrix(ncol=14, nrow=4671)
rownames(FC_full) <- all_DE_genes
HT_full <- matrix(ncol=14, nrow=250)
#HT_full <- matrix(ncol=14, nrow=4671)
rownames(HT_full) <- all_DE_genes
Der_full <- matrix(ncol=14, nrow=250)
#Der_full <- matrix(ncol=14, nrow=4671)
rownames(Der_full) <- all_DE_genes
srna_full <- matrix(ncol=14, nrow=250)
#srna_full <- matrix(ncol=14, nrow=4671)
rownames(srna_full) <- all_DE_genes
ss_full <- matrix(ncol=14, nrow=250)
#ss_full <- matrix(ncol=14, nrow=4671)
rownames(ss_full) <- all_DE_genes
for (i in seq(1:14)){
  FC_full[,i] <- FCM[,i][match(rownames(FC_full), rownames(FCM))]
  HT_full[,i] <- HTM[,i][match(rownames(HT_full), rownames(HTM))]
  Der_full[,i] <- DerM[,i][match(rownames(Der_full), rownames(DerM))]
  srna_full[,i] <- srnaM[,i][match(rownames(srna_full), rownames(srnaM))]
  ss_full[,i] <- ssM[,i][match(rownames(ss_full), rownames(ssM))]
  mirbase2[,i] <- mirbase1[,i][match(rownames(mirbase2), rownames(mirbase1))]
}  


#Set NA values to 0 count
FC_full[is.na(FC_full)] <-0
HT_full[is.na(HT_full)] <-0
Der_full[is.na(Der_full)] <-0
srna_full[is.na(srna_full)] <-0
ss_full[is.na(ss_full)] <-0
mirbase2[is.na(mirbase2)] <-0



##Correlation analysis for count matrices

for (i in seq(1:14)){
  corr_table <- cbind(FC_full[,i],HT_full[,i], Der_full[,i], srna_full[,i], ss_full[,i], mirbase2[,i])
  rownames(corr_table)<- rownames(srna_full)
  #corr_table<- log(corr_table)
  corr_table<- as.data.frame(corr_table)
  names(corr_table) <- c("featureCounts", "HTSeq", "derfinder", "srnadiff", "shortstack", "mirbase2")
  print(cor(corr_table, use="complete.obs"))
  
}




##Correlation analysis for DGE results. 
mirbase2<- matrix(ncol=5, nrow=250)
rownames(mirbase2) <- all_DE_genes
FC_full <- matrix(ncol=5, nrow=250)
#FC_full <- matrix(ncol=14, nrow=4671)
rownames(FC_full) <- all_DE_genes
HT_full <- matrix(ncol=5, nrow=250)
#HT_full <- matrix(ncol=14, nrow=4671)
rownames(HT_full) <- all_DE_genes
Der_full <- matrix(ncol=5, nrow=250)
#Der_full <- matrix(ncol=14, nrow=4671)
rownames(Der_full) <- all_DE_genes
srna_full <- matrix(ncol=5, nrow=250)
#srna_full <- matrix(ncol=14, nrow=4671)
rownames(srna_full) <- all_DE_genes
ss_full <- matrix(ncol=5, nrow=250)
#ss_full <- matrix(ncol=14, nrow=4671)
rownames(ss_full) <- all_DE_genes
for (i in seq(1:5)){
  FC_full[,i] <- FC_DME[,i][match(rownames(FC_full), rownames(FC_DME))]
  HT_full[,i] <- HT_DME[,i][match(rownames(HT_full), rownames(HT_DME))]
  Der_full[,i] <- Der_DME[,i][match(rownames(Der_full), rownames(Der_DME))]
  srna_full[,i] <- srna_DME[,i][match(rownames(srna_full), rownames(srna_DME))]
  ss_full[,i] <- ss_DME[,i][match(rownames(ss_full), rownames(ss_DME))]
  mirbase2[,i] <- mirbase_DME[,i][match(rownames(mirbase2), rownames(mirbase_DME))]
}  

for (i in seq(1:5)){
  corr_table <- cbind(FC_full[,i],HT_full[,i], Der_full[,i], srna_full[,i], ss_full[,i], mirbase2[,i])
  rownames(corr_table)<- rownames(srna_full)
  #corr_table<- log(corr_table)
  corr_table<- as.data.frame(corr_table)
  names(corr_table) <- c("featureCounts", "HTSeq", "derfinder", "srnadiff", "shortstack", "mirbase2")
  print(cor(corr_table, use="complete.obs"))
  
}
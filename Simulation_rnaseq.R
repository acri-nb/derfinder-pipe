
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
library(GGally)

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





setwd("D:/ALS_project")
files <- list.files(pattern = "\\.bam$") 
files <- files[15:26]

gffdata<-import.gff(file.choose())
gffdata1<- gffdata[seqnames(gffdata) == "chr13",]

# ##### Simulation results
# 
# gffdata1<- gffdata[seqnames(gffdata) == "chr13",]
# chr13<- Hsapiens$chr13
# ch13_a<-DNAStringSet(chr13)
# names(ch13_a) = "chr13"
# 
# 
# #subset to chromosome 13 to reduce computationnal load. Simulate a complete known annotation
# Tx_list<- seq_gtf("exons_c13_sim.gff", ch13_a, feature = "transcript", exononly = TRUE, idfield = "transcript_id", attrsep = "; ")
# 
# #Write the transcripts as fasta to use as reference with polyester
# writeXStringSet(Tx_list, "Transcripts_hg19_c13.fa", format = "fasta")
# writeXStringSet(ch13_a, "Genome_hg19_c13.fa", format = "fasta")
# 
# #Create the design matrix for polyester
# Sim_mat<- matrix(ncol=2, nrow = 4173)
# Sim_mat[,1]<-rep(1,4173)
# Sim_mat[,2]<-c(rep(0.5,834),rep(1,2505), rep(2,834))
# 
# simulate_experiment("Transcripts_hg19_c13.fa",
#                     num_reps = c(6,6),
#                     reads_per_transcript = rep(200, length(Tx_list)),
#                     fold_changes = Sim_mat,
#                     paired = FALSE,
#                     readlen = 22,
#                     distr ="normal",
#                     error_model = "uniform",
#                     error_rate = 0.01)
# 
# Der_res <- Get_Classic_Annotated_Matrix(chrs = "chr13",cutoff =1,bams = list.files(pattern="\\.sorted.bam$"),Len = 22,gffdata)


### Differential expression analysis of simulated reads ###

DGE_EdgeR<-function(mat){
  library("edgeR")
  
  design<- matrix(ncol=2,nrow=12)
  rownames(design)<-colnames(mat)
  colnames(design)<-c("group","dummy")
  design[,1]<-c(rep(0,6),rep(1,6))
  
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

# Create a design table with file name and condition (necessary for srnadiff)

SampleInfo<- matrix(nrow = 12, ncol = 3)
SampleInfo[,1]<-files
SampleInfo[,2]<-files
SampleInfo[,3]<-c(rep(0,6),rep(1,6))
colnames(SampleInfo)<-c("FileName", "SampleName", "Condition")
SampleInfo <- as.data.frame(SampleInfo)
annotReg   <- readAnnotation("D:/ALS_project/exons_c13_sim.gff", feature="transcript", source=NULL)
srnaExp    <- srnadiffExp(files, SampleInfo, annotReg)
srnaExp    <- srnadiff(srnaExp,segMethod = "annotation")
CM<- srnaExp@countMatrix
CM<- CM[order(rownames(CM)),]
trs<- srnaExp@annotReg[names(srnaExp@annotReg) %in% rownames(CM),]
trs<-cbind(names(trs),trs$transcript_id)
colnames(trs) <- c("Internal_ids", "ENST_ids")
trs <- data.frame(trs)
trs<- trs[order(trs$Internal_ids),]
rownames(CM)<-trs$ENST_ids
res_srna<-DGE_EdgeR(CM)

#upload previous derfinder results, filter low counts and get DGE results

raw_der <-read.table("D:/ALS_project/Full_sim_20191218.txt", header = T, stringsAsFactors = F)
raw_der <- raw_der[!(is.na(raw_der$transcr_names)),]
raw_counts <- as.matrix(raw_der[,9:20])
rownames(raw_counts)<-raw_der$transcr_names
raw_counts_der <- raw_counts[rowMeans(raw_counts)>20,]
Der_DE<-DGE_EdgeR(raw_counts_der)

#upload previous featureCounts results, filter low counts and get DGE results

FC_res<-read.table("D:/ALS_project/FC_sim_Bowtie_20191218.txt", header = T, stringsAsFactors = F)
raw_counts <-as.matrix(FC_res[7:18])
rownames(raw_counts)<-FC_res$Geneid
raw_counts_FC <- raw_counts[rowMeans(raw_counts)>20,]
raw_counts_FC<-raw_counts_FC[rownames(raw_counts_FC) %in% gffdata$transcript_id,]
FC_DE<-DGE_EdgeR(raw_counts_FC)

#upload previous HTSeq results, filter low counts and get DGE results

HT_res<-read.table("D:/ALS_project/HT_sim_20191218.txt", header = F, stringsAsFactors = F)
names(HT_res)<-c("Geneid",names(FC_res[7:18]))
raw_counts <-as.matrix(HT_res[2:13])
rownames(raw_counts)<-HT_res$Geneid
raw_counts_HT <- raw_counts[rowMeans(raw_counts)>20,]
raw_counts_HT<-raw_counts_HT[rownames(raw_counts_HT) %in% gffdata$transcript_id,]
HT_DE<-DGE_EdgeR(raw_counts_HT)

#Filter out low-expression genes for srnadiff

CM1 <- CM[rowMeans(CM)>20,]

venn.diagram(
  x = list(rownames(raw_counts_FC), rownames(raw_counts_HT), rownames(raw_counts_der), rownames(CM1)),
  category.names = c("featureCounts" , "HTSeq" , "derfinder", "srnadiff"),
  filename = 'venn_diagramm_sim.png',
  output=TRUE
)

# get only the intersecting expressed RNAs for comparisons

int1<- intersect(rownames(CM1), rownames(raw_counts_der))
raw_counts_der <- raw_counts_der[rownames(raw_counts_der) %in% int1,]
raw_counts_der <- raw_counts_der[order(rownames(raw_counts_der)),]
raw_counts_FC <- raw_counts_FC[rownames(raw_counts_FC) %in% int1,]
raw_counts_FC <- raw_counts_FC[order(rownames(raw_counts_FC)),]
raw_counts_HT <- raw_counts_HT[rownames(raw_counts_HT) %in% int1,]
raw_counts_HT <- raw_counts_HT[order(rownames(raw_counts_HT)),]
CM1 <- CM1[rownames(CM1) %in% int1,]
CM1 <- CM1[order(rownames(CM1)),]

#Get correlation matrices / tables
for (i in seq(1:12)){
  corr_table <- cbind(raw_counts_FC[,i],raw_counts_HT[,i], raw_counts_der[,i], CM1[,i])
  rownames(corr_table)<- rownames(CM1)
  corr_table<- log(corr_table)
  corr_table<- as.data.frame(corr_table)
  names(corr_table) <- c("featureCounts", "HTSeq", "derfinder", "srnadiff")
  print(cor(corr_table, method = "spearman"))
  
}


# ### ANALYSIS OF PERFORMANCE METRICS FOR SIMULATED READS RESULTS ###
Truth <-read.table("D:/ALS_project/sim_tx_info.txt", sep = "\t", stringsAsFactors = F, header = TRUE)
#FC_DE<-read.table(file.choose(), sep = "\t", stringsAsFactors = F, header = TRUE)
#HT_DE<-read.table(file.choose(), sep = "\t", stringsAsFactors = F, header = TRUE)
# Der_DE<-read.table(file.choose(), sep = "\t", stringsAsFactors = F, header = TRUE)
# 


#Set dummy factor based on significance threshold
FC_DE$is_de<-cut(FC_DE$FDR,breaks = c(-Inf,0.05,Inf),labels =c("TRUE","FALSE"))
HT_DE$is_de<-cut(HT_DE$FDR,breaks = c(-Inf,0.05,Inf),labels =c("TRUE","FALSE"))
Der_DE$is_de<-cut(Der_DE$FDR,breaks = c(-Inf,0.05,Inf),labels =c("TRUE","FALSE"))
res_srna$is_de <- cut(res_srna$FDR,breaks = c(-Inf,0.05,Inf),labels =c("TRUE","FALSE"))



#Get only intersecting transcript names for all tools 
gene_names<- intersect(intersect(rownames(Der_DE), rownames(HT_DE)), intersect(rownames(FC_DE), rownames(res_srna)))

FC_DE<-FC_DE[rownames(FC_DE) %in% gene_names,]
HT_DE<-HT_DE[rownames(HT_DE) %in% gene_names,]
Der_DE<-Der_DE[rownames(Der_DE) %in% gene_names,]
res_srna1<- res_srna[rownames(res_srna) %in% gene_names,]


#Order rownames
FC_DE<-FC_DE[order(rownames(FC_DE)),]
HT_DE<-HT_DE[order(rownames(HT_DE)),]
Der_DE<-Der_DE[order(rownames(Der_DE)),]
res_srna1<-res_srna1[order(rownames(res_srna1)),]

Truth<-Truth[Truth$transcriptid %in% gene_names,]
Truth<-Truth[order(Truth$transcriptid),]


#combine all data into a single dataframe for plotting

All_data<-cbind(Truth$transcriptid,Truth$DEstatus.2,HT_DE$is_de,FC_DE$is_de, Der_DE$is_de, res_srna1$is_de)
replace<-ifelse(All_data[,2] == "TRUE","1","2")
colnames(All_data)<-c("id","truth","HTseq","FeatCount","DerFinder", "srnadiff")
All_data[,2]<-replace

#Calculate sensitivity and specificity

sensitivity(table(All_data[,2],All_data[,3]))
specificity(table(All_data[,2],All_data[,3]))

sensitivity(table(All_data[,2],All_data[,4]))
specificity(table(All_data[,2],All_data[,4]))

sensitivity(table(All_data[,2],All_data[,5]))
specificity(table(All_data[,2],All_data[,5]))

sensitivity(table(All_data[,2],All_data[,6]))
specificity(table(All_data[,2],All_data[,6]))

FC_DE <- FC_DE[order(rownames(FC_DE)),]
HT_DE <- HT_DE[order(rownames(HT_DE)),]
Der_DE <- Der_DE[order(rownames(Der_DE)),]
res_srna1 <- res_srna1[order(rownames(res_srna1)),]


# plot correlation scatter plots, taking NA values into account. 
grphic<- FC_DE
grphic$HT <- rep(NA, length(grphic$logFC))
grphic$Der <- rep(NA, length(grphic$logFC))
grphic$srna <- rep(NA, length(grphic$logFC))

rn<- intersect(rownames(HT_DE), rownames(FC_DE))
rpl <- HT_DE$logFC
grphic$HT[rownames(grphic) %in% rn] <- rpl
rn<- intersect(rownames(FC_DE), rownames(Der_DE))
rpl <- Der_DE$logFC
grphic$Der[rownames(grphic) %in% rn] <- rpl
rn<- intersect(rownames(FC_DE), rownames(res_srna1))
rpl <- res_srna1$logFC
grphic$srna[rownames(grphic) %in% rn] <- rpl
grphic<- grphic[,c(1,7,8,9)]
names(grphic) <- c("featureCounts", "HTSeq", "derfinder", "srnadiff")


ggpairs(grphic)


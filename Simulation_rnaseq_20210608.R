
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
gffdata1<- gffdata1[gffdata1$type == "transcript",]

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


##shortstack##
ss<- read.table(file.choose(),header = T, sep = "\t", row.names = 1)

ssplit1<- colsplit(rownames(ss), ":", names=c("chr","start"))
ssplit2<- colsplit(ssplit1$start, "-", names=c("start","end"))
ss <- cbind(ssplit1$chr, ssplit2, ss[,3:14])
#col.order <- c("ssplit1$chr", "start", "end", "sample_01", "sample_02", "sample_03", "sample_04", "sample_05", "sample_06", "sample_07", "sample_08", "sample_09", "sample_10", "sample_11", "sample_12")
#ss<- ss[,col.order]

#find overlapping regions with annotation
ss_gr<- makeGRangesFromDataFrame(ss, keep.extra.columns = T, seqnames.field = "ssplit1$chr")
ss_ol<- subsetByOverlaps.keepAllMeta(gffdata1, ss_gr)
mt_dfr<-cbind(as.character(seqnames(ss_ol)),start(ss_ol),end(ss_ol))
colnames(mt_dfr)<-c("contig","start","end")
mt_dfr<-cbind(mt_dfr,values(ss_ol))
print(head(mt_dfr))

#Get sum of counts per gene
summing<-function(x) sum(na.omit(sapply(strsplit(x,";"),as.integer))) 
for (i in names(ss)[3:15]){
  sums<-sapply(mt_dfr[[i]], summing)
  #replace with new values
  mt_dfr[[i]]<-sums
}
mt_dfr$contig<-as.character(mt_dfr$contig)
mt_dfr$start<-as.character(mt_dfr$start)
mt_dfr$end<-as.character(mt_dfr$end)

#DGE analysis
sstack<- mt_dfr
rownames(sstack) <- sstack$transcript_id
sstack_counts <- sstack[,26:37]
sstack_counts <- as.matrix(sstack_counts)
sstack_counts <- sstack_counts[rowMeans(sstack_counts) > 20,]
sstack_counts <- as.data.frame(sstack_counts)
ss_DGE <- DGE_EdgeR(sstack_counts)





venn.diagram(
  x = list(rownames(raw_counts_FC), rownames(raw_counts_HT), rownames(raw_counts_der), rownames(CM1), rownames(sstack_counts)),
  category.names = c("featureCounts" , "HTSeq" , "derfinder", "srnadiff", "shortstack"),
  filename = 'venn_sim_w_ss.png',
  output=TRUE
)

# get only the intersecting expressed RNAs for comparisons


#combine all DGE results into one DF
all_DE_genes<- unique(c(rownames(raw_counts_FC),rownames(raw_counts_HT),rownames(sstack_counts), rownames(CM1), rownames(raw_counts_der) ))
combi <- matrix(ncol = 5, nrow = 4093)
rownames(combi) <- all_DE_genes        
colnames(combi) <- c("FC", "HT", "Der", "srna", "sstack")
combi[,1] <- FC_DE$logFC[match(rownames(combi), rownames(FC_DE))]
combi[,2] <- HT_DE$logFC[match(rownames(combi), rownames(HT_DE))]
combi[,3] <- Der_DE$logFC[match(rownames(combi), rownames(Der_DE))]
combi[,4] <- res_srna$logFC[match(rownames(combi), rownames(res_srna))]
combi[,5] <- ss_DGE$logFC[match(rownames(combi), rownames(ss_DGE))]

grphic <- as.data.frame(combi)
combi[is.na(combi)] <- 0
names(grphic) <- c("featureCounts", "HTSeq", "derfinder", "srnadiff", "shortstack")
ggpairs(grphic) #save 800x800





FC_full <- matrix(ncol=12, nrow=4093)
rownames(FC_full) <- all_DE_genes
HT_full <- matrix(ncol=12, nrow=4093)
rownames(HT_full) <- all_DE_genes
Der_full <- matrix(ncol=12, nrow=4093)
rownames(Der_full) <- all_DE_genes
srna_full <- matrix(ncol=12, nrow=4093)
rownames(srna_full) <- all_DE_genes
ss_full <- matrix(ncol=12, nrow=4093)
rownames(ss_full) <- all_DE_genes
for (i in seq(1:12)){
  FC_full[,i] <- raw_counts_FC[,i][match(rownames(FC_full), rownames(raw_counts_FC))]
  HT_full[,i] <- raw_counts_HT[,i][match(rownames(HT_full), rownames(raw_counts_HT))]
  Der_full[,i] <- raw_counts_der[,i][match(rownames(Der_full), rownames(raw_counts_der))]
  srna_full[,i] <- CM1[,i][match(rownames(srna_full), rownames(CM1))]
  ss_full[,i] <- sstack_counts[,i][match(rownames(ss_full), rownames(sstack_counts))]
}  

FC_full[is.na(FC_full)] <-0
HT_full[is.na(HT_full)] <-0
Der_full[is.na(Der_full)] <-0
srna_full[is.na(srna_full)] <-0
ss_full[is.na(ss_full)] <-0


#Get correlation matrices / tables
for (i in seq(1:12)){
  corr_table <- cbind(FC_full[,i],HT_full[,i], Der_full[,i], srna_full[,i], ss_full[,1])
  rownames(corr_table)<- rownames(Der_full)
  corr_table<- log(corr_table)
  corr_table<- as.data.frame(corr_table)
  names(corr_table) <- c("featureCounts", "HTSeq", "derfinder", "srnadiff", "shortstack")
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
ss_DGE$is_de <- cut(ss_DGE$FDR,breaks = c(-Inf,0.05,Inf),labels =c("TRUE","FALSE"))

combi <- matrix(ncol = 7, nrow = 4173)
rownames(combi) <- Truth$transcriptid        
colnames(combi) <- c("ID","is.de","FC", "HT", "Der", "srna", "sstack")
combi[,1] <- Truth$transcriptid
combi[,2] <- Truth$DEstatus.2
combi[,3] <- FC_DE$is_de[match(rownames(combi), rownames(FC_DE))]
combi[,4] <- HT_DE$is_de[match(rownames(combi), rownames(HT_DE))]
combi[,5] <- Der_DE$is_de[match(rownames(combi), rownames(Der_DE))]
combi[,6] <- res_srna$is_de[match(rownames(combi), rownames(res_srna))]
combi[,7] <- ss_DGE$is_de[match(rownames(combi), rownames(ss_DGE))]



#combine all data into a single dataframe for plotting

All_data<-combi
replace<-ifelse(All_data[,2] == "TRUE","1","2")
colnames(All_data)<-c("id","truth","FeatCount","HTseq","DerFinder", "srnadiff", "shortstack")
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

sensitivity(table(All_data[,2],All_data[,7]))
specificity(table(All_data[,2],All_data[,7]))







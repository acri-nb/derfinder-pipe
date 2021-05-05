


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



######## For real data ############

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


setwd("D:/ALS_project/Analysis20200630")
FC_counts<- read.table("FC_Bowtie_gencode19_2020.txt", header = T, row.names = 1)
FC_counts<- FC_counts[,-c(1,2,3,4,5)]
FC_counts<- FC_counts[!(rownames(FC_counts) %in% supp_genes),]
FC_counts<- FC_counts[rowMeans(FC_counts) > 20,]
FC_DGE<- DGE_EdgeR(FC_counts)

HT_counts<- read.table("HT_gencode19_2020.txt", header = T, row.names = 1)
HT_counts<- HT_counts[!(rownames(HT_counts) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")),]
HT_counts<- HT_counts[!(rownames(HT_counts) %in% supp_genes),]
HT_counts<- HT_counts[rowMeans(HT_counts) > 20,]
HT_DGE<- DGE_EdgeR(HT_counts)

#Der_counts<- read.table("Der_counts_20200702.txt", header = T, sep= "\t")
Der_counts<- read.table("Derfinder_nocutoff_nonorm_20200707.txt", header = T, sep= "\t")
rownames(Der_counts)<- Der_counts$gene_id
Der_counts<- Der_counts[26:39]
Der_counts<- Der_counts[rowMeans(Der_counts) > 20,]
Der_DGE<- DGE_EdgeR(Der_counts)

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
srna_DGE<- DGE_EdgeR(srna_counts)

FC_DGE <- FC_DGE[order(rownames(FC_DGE)),]
HT_DGE <- HT_DGE[order(rownames(HT_DGE)),]
Der_DGE <- Der_DGE[order(rownames(Der_DGE)),]
srna_DGE <- srna_DGE[order(rownames(srna_DGE)),]

inter<- intersect(intersect(rownames(srna_DGE), rownames(Der_DGE)),rownames(FC_DGE))
grphic<- cbind(FC_DGE$logFC[rownames(FC_DGE) %in% inter], HT_DGE$logFC[rownames(HT_DGE) %in% inter], Der_DGE$logFC[rownames(Der_DGE) %in% inter], srna_DGE$logFC[rownames(srna_DGE) %in% inter])
grphic <- as.data.frame(grphic)
names(grphic) <- c("featureCounts", "HTSeq", "derfinder", "srnadiff")
ggpairs(grphic) #save 800x800

FC_counts <- FC_counts[order(rownames(FC_counts)) & rownames(FC_counts) %in% inter,]
HT_counts <- HT_counts[order(rownames(HT_counts)) & rownames(HT_counts) %in% inter,]
Der_counts <- Der_counts[order(rownames(Der_counts)) & rownames(Der_counts) %in% inter,]
srna_counts <- srna_counts[order(rownames(srna_counts)) & rownames(srna_counts) %in% inter,]




for (i in seq(1:14)){
  corr_table <- cbind(FC_counts[,i],HT_counts[,i], Der_counts[,i], srna_counts[,i])
  rownames(corr_table)<- rownames(srna_counts)
  corr_table<- log(corr_table)
  corr_table<- as.data.frame(corr_table)
  names(corr_table) <- c("featureCounts", "HTSeq", "derfinder", "srnadiff")
  print(cor(corr_table, method = "spearman"))
  
}

venn.diagram(
  x = list(rownames(FC_counts), rownames(HT_counts), rownames(Der_counts), rownames(srna_counts)),
  category.names = c("featureCounts" , "HTSeq" , "derfinder", "srnadiff"),
  filename = 'venn_diagramm_als.png',
  output=TRUE
)


#Plot RNA types abundances

gffdata$gene_type[gffdata$gene_type %in% c("Mt_tRNA", "Mt_rRNA")] <- "mitochondrial"
gffdata$gene_type[gffdata$gene_type %in% c("polymorphic_pseudogene")] <- "pseudogene"
gffdata$gene_type[gffdata$gene_type %in% c("sense_intronic", "sense_overlapping", "3prime_overlapping_ncrna")] <- "other"
gffdata$gene_type[gffdata$gene_type %in% c("IG_C_gene", "IG_C_pseudogene", "IG_D_gene", "IG_J_gene", "IG_J_pseudogene", "IG_V_gene", "IG_V_pseudogene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_J_pseudogene","TR_V_gene", "TR_V_pseudogene")] <- "other"

# Build a table for plotting containing gene-wise means and gene name
FC_means<- cbind(rownames(FC_counts), rowMeans(FC_counts))
FC_means<- FC_means[order(rownames(FC_means)),]
HT_means<- cbind(rownames(HT_counts), rowMeans(HT_counts))
HT_means<- HT_means[order(rownames(HT_means)),]
Der_means<- cbind(rownames(Der_counts), rowMeans(Der_counts))
Der_means<- Der_means[order(rownames(Der_means)),]
srna_means<- cbind(rownames(srna_counts), rowMeans(srna_counts))
srna_means<- srna_means[order(rownames(srna_means)),]


#Get the corresponding data for each gene in the gff annotation
gf_FC<- gffdata[gffdata$gene_id %in% FC_means,]
gf_FC<- gf_FC[order(gf_FC$gene_id),]
gf_HT<- gffdata[gffdata$gene_id %in% HT_means,]
gf_HT<- gf_HT[order(gf_HT$gene_id),]
gf_Der<- gffdata[gffdata$gene_id %in% Der_means,]
gf_Der<- gf_Der[order(gf_Der$gene_id),]
gf_srna<- gffdata[gffdata$gene_id %in% srna_means,]
gf_srna<- gf_srna[order(gf_srna$gene_id),]

# Add the gene_type data from the gff annotation to the means 
FC_means<- as.data.frame(cbind(FC_means, gf_FC$gene_type))
HT_means<- as.data.frame(cbind(HT_means, gf_HT$gene_type))
Der_means<- as.data.frame(cbind(Der_means, gf_Der$gene_type))
srna_means<- as.data.frame(cbind(srna_means, gf_srna$gene_type))

#Change means to numeric...

FC_means$V2 <- as.numeric(FC_means$V2)
HT_means$V2 <- as.numeric(HT_means$V2)
Der_means$V2 <- as.numeric(Der_means$V2)
srna_means$V2 <- as.numeric(srna_means$V2)

# aggregate counts by gene
FC_m<-aggregate(FC_means$V2, by=list(Players=FC_means$V3), FUN=sum)
HT_m<-aggregate(HT_means$V2, by=list(Players=HT_means$V3), FUN=sum)
Der_m<-aggregate(Der_means$V2, by=list(Players=Der_means$V3), FUN=sum)
srna_m<-aggregate(srna_means$V2, by=list(Players=srna_means$V3), FUN=sum)

#Construct final tables for plotting
final_type<- cbind(FC_m, HT_m$x, Der_m$x, srna_m$x)
names(final_type)<- c("RNA_Type", "featureCounts", "HTSeq", "derfinder", "srnadiff")







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
#normFactors(srnaExp)<-c(rep(1,14))



srnaExp    <- srnadiff(srnaExp,segMethod = "annotation")
#parameters(srnaExp) <- list("minSize"=0,"minGap"=0)
#srnaExp0 <- srnadiff(srnaExp,segMethod = "annotation")

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


# ####### For enrichment analysis ###############################

#Remove trailing version number from ensembl annotation
names<- colsplit(rownames(Der_DGE),"[.]",names = c("col1","col2"))
rownames(Der_DGE)<-names$col1
names<- colsplit(rownames(HT_DGE),"[.]",names = c("col1","col2"))
rownames(HT_DGE)<-names$col1
names<- colsplit(rownames(FC_DGE),"[.]",names = c("col1","col2"))
rownames(FC_DGE)<-names$col1
names<- colsplit(rownames(srna_DGE),"[.]",names = c("col1","col2"))
rownames(srna_DGE)<-names$col1


#For hypergeometric testing, we took only significant results from each tool, and used the ENTREZ ids as input

enrich_FC<-bitr(rownames(FC_DGE[FC_DGE$FDR<0.05,]), fromType = "ENSEMBL",toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop=TRUE)
enrich_HT<-bitr(rownames(HT_DGE[HT_DGE$FDR<0.05,]), fromType = "ENSEMBL",toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop=TRUE)
enrich_Der<-bitr(rownames(Der_DGE[Der_DGE$FDR<0.05,]), fromType = "ENSEMBL",toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop=TRUE)
enrich_srna<-bitr(rownames(srna_DGE[srna_DGE$FDR<0.05,]), fromType = "ENSEMBL",toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop=TRUE)

geneListfc = enrich_FC$ENTREZID
geneListfc = sort(geneListfc, decreasing = TRUE)
geneListht = enrich_HT$ENTREZID
geneListht = sort(geneListht, decreasing = TRUE)
geneListder = enrich_Der$ENTREZID
geneListder = sort(geneListder, decreasing = TRUE)
geneListsrna = enrich_srna$ENTREZID
geneListsrna = sort(geneListsrna, decreasing = TRUE)

pwefc <- enrichKEGG(gene=geneListfc, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5)
gofc <- enrichGO(gene=geneListfc, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5)
#head(summary(pwefc))
pweht <- enrichKEGG(gene=geneListht, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5)
goht <- enrichGO(gene=geneListht, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5)
#head(summary(pweht))
pweder <- enrichKEGG(gene=geneListder, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5)
goder <- enrichGO(gene=geneListder, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5)
#head(summary(pweder))
pwesrna <- enrichKEGG(gene=geneListsrna, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5)
gosrna <- enrichGO(gene=geneListsrna, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5)
#head(summary(pwesrna))

#merge dataframes for correlation analysis
m1<- merge(goder@result, gosrna@result, by.x="ID", by.y="ID", all.x = T, all.y = T)
m2<- merge(gofc@result, goht@result, by.x="ID", by.y="ID", all.x = T, all.y = T)
m3<- merge(m2, m1, by.x="ID", by.y="ID", all.x = T, all.y = T)
m4 <- m3[,c(1,6,15,24,33)]
names(m4) <- c("ID", "p_fc", "p_ht", "p_der", "p_srna")
cor(m4[,2:5], use = "pairwise.complete.obs", method = "spearman")

m1<- merge(pweder@result, pwesrna@result, by.x="ID", by.y="ID", all.x = T, all.y = T)
m2<- merge(pwefc@result, pweht@result, by.x="ID", by.y="ID", all.x = T, all.y = T)
m3<- merge(m2, m1, by.x="ID", by.y="ID", all.x = T, all.y = T)
m4 <- m3[,c(1,5,13,21,29)]
names(m4) <- c("ID", "p_fc", "p_ht", "p_der", "padj_srna")
cor(m4[,2:5], use = "pairwise.complete.obs", method = "spearman")




#For GSEA, each tool must have an ordered ranked list to submit to gseGO or gseKEGG
enrich_FCnames<-bitr(rownames(FC_DGE), fromType = "ENSEMBL",toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop=TRUE)
enrich_FCnames <- enrich_FCnames[!(duplicated(enrich_FCnames$ENSEMBL)),]
enrich_FC <- FC_DGE[rownames(FC_DGE) %in% enrich_FCnames$ENSEMBL,]
vect0 <- sign(enrich_FC$logFC)*(-log10(enrich_FC$FDR))
vect0<-setNames(vect0, enrich_FCnames$ENTREZID)
vect0 = sort(vect0, decreasing = TRUE)
gofc <- gseGO(geneList=vect0, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 1, minGSSize = 5, pAdjustMethod = "BH")
pwefc <- gseKEGG(geneList=vect0, organism = "hsa", pvalueCutoff = 1, minGSSize = 5, pAdjustMethod = "BH")

enrich_HTnames<-bitr(rownames(HT_DGE), fromType = "ENSEMBL",toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop=TRUE)
enrich_HTnames <- enrich_HTnames[!(duplicated(enrich_HTnames$ENSEMBL)),]
enrich_HT <- HT_DGE[rownames(HT_DGE) %in% enrich_HTnames$ENSEMBL,]
vect0 <- sign(enrich_HT$logFC)*(-log10(enrich_HT$FDR))
vect0<-setNames(vect0, enrich_HTnames$ENTREZID)
vect0 = sort(vect0, decreasing = TRUE)
goht <- gseGO(geneList=vect0, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 1, minGSSize = 5, pAdjustMethod = "BH")
pweht <- gseKEGG(geneList=vect0, organism = "hsa", pvalueCutoff = 1, minGSSize = 5, pAdjustMethod = "BH")

enrich_dernames<-bitr(rownames(Der_DGE), fromType = "ENSEMBL",toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop=TRUE)
enrich_dernames <- enrich_dernames[!(duplicated(enrich_dernames$ENSEMBL)),]
enrich_der <- Der_DGE[rownames(Der_DGE) %in% enrich_dernames$ENSEMBL,]
vect0 <- sign(enrich_der$logFC)*(-log10(enrich_der$FDR))
vect0<-setNames(vect0, enrich_dernames$ENTREZID)
vect0 = sort(vect0, decreasing = TRUE)
goder <- gseGO(geneList=vect0, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 1, minGSSize = 5, pAdjustMethod = "BH")
pweder <- gseKEGG(geneList=vect0, organism = "hsa", pvalueCutoff = 1, minGSSize = 5, pAdjustMethod = "BH")

enrich_srnanames<-bitr(rownames(srna_DGE), fromType = "ENSEMBL",toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop=TRUE)
enrich_srnanames <- enrich_srnanames[!(duplicated(enrich_srnanames$ENSEMBL)),]
enrich_srna <- srna_DGE[rownames(srna_DGE) %in% enrich_srnanames$ENSEMBL,]
vect0 <- sign(enrich_srna$logFC)*(-log10(enrich_srna$FDR))
vect0<-setNames(vect0, enrich_srnanames$ENTREZID)
vect0 = sort(vect0, decreasing = TRUE)
gosrna <- gseGO(geneList=vect0, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 1, minGSSize = 5, pAdjustMethod = "BH")
pwesrna <- gseKEGG(geneList=vect0, organism = "hsa", pvalueCutoff = 1, minGSSize = 5, pAdjustMethod = "BH")


#Merge dataframes for correlation analyses
m1<- merge(goder@result, gosrna@result, by.x="ID", by.y="ID", all.x = T, all.y = T)
m2<- merge(gofc@result, goht@result, by.x="ID", by.y="ID", all.x = T, all.y = T)
m3<- merge(m2, m1, by.x="ID", by.y="ID", all.x = T, all.y = T)
m4 <- m3[,c(1,7,18,29,40)]
names(m4) <- c("ID", "p_fc", "p_ht", "p_der", "p_srna")
cor(m4[,2:5], use = "pairwise.complete.obs", method = "spearman")

m1<- merge(pweder@result, pwesrna@result, by.x="ID", by.y="ID", all.x = T, all.y = T)
m2<- merge(pwefc@result, pweht@result, by.x="ID", by.y="ID", all.x = T, all.y = T)
m3<- merge(m2, m1, by.x="ID", by.y="ID", all.x = T, all.y = T)
m4 <- m3[,c(1,6,16,26,36)]
names(m4) <- c("ID", "p_fc", "p_ht", "p_der", "p_srna")
cor(m4[,2:5], use = "pairwise.complete.obs", method = "spearman")
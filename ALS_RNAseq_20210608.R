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

#Create count matrices and DGE results for each tool
##FeatureCounts
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

#for srnadiff
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
srna_counts <- as.data.frame(srna_counts)
srna_counts <- srna_counts[,names(FC_counts)]
CM <- CM[,names(FC_counts)]
srna_DGE<- DGE_EdgeR(srna_counts)


##shortstack##
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

#Get sum of counts per gene
summing<-function(x) sum(na.omit(sapply(strsplit(x,";"),as.integer))) 
for (i in names(ss)[3:17]){
  sums<-sapply(mt_dfr[[i]], summing)
  #replace with new values
  mt_dfr[[i]]<-sums
}
mt_dfr$contig<-as.character(mt_dfr$contig)
mt_dfr$start<-as.character(mt_dfr$start)
mt_dfr$end<-as.character(mt_dfr$end)

#DGE analysis
sstack<- mt_dfr
rownames(sstack) <- sstack$gene_id
sstack_counts <- sstack[,25:38]
sstack_counts <- as.matrix(sstack_counts)
sstack_counts <- sstack_counts[rowMeans(sstack_counts) > 20,]
sstack_counts <- as.data.frame(sstack_counts)
ss_DGE <- DGE_EdgeR(sstack_counts)



#order rownames

FC_DGE <- FC_DGE[order(rownames(FC_DGE)),]
HT_DGE <- HT_DGE[order(rownames(HT_DGE)),]
Der_DGE <- Der_DGE[order(rownames(Der_DGE)),]
srna_DGE <- srna_DGE[order(rownames(srna_DGE)),]
ss_DGE <- ss_DGE[order(rownames(ss_DGE)),]


#combine all DGE results into one DF
all_DE_genes<- unique(c(rownames(FC_DGE),rownames(HT_DGE),rownames(ss_DGE), rownames(srna_DGE), rownames(Der_DGE) ))
combi <- matrix(ncol = 5, nrow = 4680)
rownames(combi) <- all_DE_genes        
colnames(combi) <- c("FC", "HT", "Der", "srna", "sstack")
combi[,1] <- FC_DGE$logFC[match(rownames(combi), rownames(FC_DGE))]
combi[,2] <- HT_DGE$logFC[match(rownames(combi), rownames(HT_DGE))]
combi[,3] <- Der_DGE$logFC[match(rownames(combi), rownames(Der_DGE))]
combi[,4] <- srna_DGE$logFC[match(rownames(combi), rownames(srna_DGE))]
combi[,5] <- ss_DGE$logFC[match(rownames(combi), rownames(ss_DGE))]


# plot the effect sizes for combined dataframe

grphic <- as.data.frame(combi)
combi[is.na(combi)] <- 0
names(grphic) <- c("featureCounts", "HTSeq", "derfinder", "srnadiff", "shortstack")
library(GGally)
ggpairs(grphic) #save 800x800


FC_full <- matrix(ncol=14, nrow=length(all_DE_genes))
rownames(FC_full) <- all_DE_genes
HT_full <- matrix(ncol=14, nrow=length(all_DE_genes))
rownames(HT_full) <- all_DE_genes
Der_full <- matrix(ncol=14, nrow=length(all_DE_genes))
rownames(Der_full) <- all_DE_genes
srna_full <- matrix(ncol=14, nrow=length(all_DE_genes))
rownames(srna_full) <- all_DE_genes
ss_full <- matrix(ncol=14, nrow=length(all_DE_genes))
rownames(ss_full) <- all_DE_genes
for (i in seq(1:14)){
  FC_full[,i] <- FC_counts[,i][match(rownames(FC_full), rownames(FC_counts))]
  HT_full[,i] <- HT_counts[,i][match(rownames(HT_full), rownames(HT_counts))]
  Der_full[,i] <- Der_counts[,i][match(rownames(Der_full), rownames(Der_counts))]
  srna_full[,i] <- srna_counts[,i][match(rownames(srna_full), rownames(srna_counts))]
  ss_full[,i] <- sstack_counts[,i][match(rownames(ss_full), rownames(sstack_counts))]
}  

FC_full[is.na(FC_full)] <-0
HT_full[is.na(HT_full)] <-0
Der_full[is.na(Der_full)] <-0
srna_full[is.na(srna_full)] <-0
ss_full[is.na(ss_full)] <-0


## correlation analysis with spearman
for (i in seq(1:14)){
  corr_table <- cbind(FC_full[,i],HT_full[,i], Der_full[,i], srna_full[,i], ss_full[,i])
  rownames(corr_table)<- rownames(srna_full)
  #corr_table<- log(corr_table)
  corr_table<- as.data.frame(corr_table)
  names(corr_table) <- c("featureCounts", "HTSeq", "derfinder", "srnadiff", "shortstack")
  print(cor(corr_table,method = "spearman", use="complete.obs"))
  
}

##Venn diagrams
venn.diagram(
  x = list(rownames(FC_counts), rownames(HT_counts), rownames(Der_counts), rownames(srna_counts), rownames(sstack_counts)),
  category.names = c("featureCounts" , "HTSeq" , "derfinder", "srnadiff", "shortstack"),
  filename = 'venn_diagramm_als_ss2.png',
  output=TRUE
)


#Plot RNA types abundances

gffdata$gene_type[gffdata$gene_type %in% c("Mt_tRNA", "Mt_rRNA")] <- "mitochondrial"
gffdata$gene_type[gffdata$gene_type %in% c("polymorphic_pseudogene")] <- "pseudogene"
gffdata$gene_type[gffdata$gene_type %in% c("sense_intronic", "sense_overlapping", "3prime_overlapping_ncrna")] <- "other"
gffdata$gene_type[gffdata$gene_type %in% c("IG_C_gene", "IG_C_pseudogene", "IG_D_gene", "IG_J_gene", "IG_J_pseudogene", "IG_V_gene", "IG_V_pseudogene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_J_pseudogene","TR_V_gene", "TR_V_pseudogene")] <- "other"

# Build a table for plotting containing gene-wise means and gene name
FC_means<- cbind(rownames(FC_full), rowMeans(FC_full))
FC_means<- FC_means[order(rownames(FC_means)),]
HT_means<- cbind(rownames(HT_full), rowMeans(HT_full))
HT_means<- HT_means[order(rownames(HT_means)),]
Der_means<- cbind(rownames(Der_full), rowMeans(Der_full))
Der_means<- Der_means[order(rownames(Der_means)),]
srna_means<- cbind(rownames(srna_full), rowMeans(srna_full))
srna_means<- srna_means[order(rownames(srna_means)),]
ss_means <- cbind(rownames(data.frame(ss_full)), rowMeans(data.frame(ss_full)))
ss_means <- ss_means[order(rownames(ss_means)),]

#Get the corresponding data for each gene in the gff annotation
gf_FC<- gffdata[gffdata$gene_id %in% FC_means,]
gf_FC<- gf_FC[order(gf_FC$gene_id),]
gf_HT<- gffdata[gffdata$gene_id %in% HT_means,]
gf_HT<- gf_HT[order(gf_HT$gene_id),]
gf_Der<- gffdata[gffdata$gene_id %in% Der_means,]
gf_Der<- gf_Der[order(gf_Der$gene_id),]
gf_srna<- gffdata[gffdata$gene_id %in% srna_means,]
gf_srna<- gf_srna[order(gf_srna$gene_id),]
gf_ss <- gffdata[gffdata$gene_id %in% ss_means,]
gf_ss <- gf_ss[order(gf_ss$gene_id),]

# Add the gene_type data from the gff annotation to the means 
FC_means<- as.data.frame(cbind(FC_means, gf_FC$gene_type))
HT_means<- as.data.frame(cbind(HT_means, gf_HT$gene_type))
Der_means<- as.data.frame(cbind(Der_means, gf_Der$gene_type))
srna_means<- as.data.frame(cbind(srna_means, gf_srna$gene_type))
ss_means<- as.data.frame(cbind(ss_means, gf_ss$gene_type))

#Change means to numeric...

FC_means$V2 <- as.numeric(FC_means$V2)
HT_means$V2 <- as.numeric(HT_means$V2)
Der_means$V2 <- as.numeric(Der_means$V2)
srna_means$V2 <- as.numeric(srna_means$V2)
ss_means$V2 <- as.numeric(ss_means$V2)

# aggregate counts by gene
FC_m<-aggregate(FC_means$V2, by=list(Players=FC_means$V3), FUN=sum)
HT_m<-aggregate(HT_means$V2, by=list(Players=HT_means$V3), FUN=sum)
Der_m<-aggregate(Der_means$V2, by=list(Players=Der_means$V3), FUN=sum)
srna_m<-aggregate(srna_means$V2, by=list(Players=srna_means$V3), FUN=sum)
ss_m<-aggregate(ss_means$V2, by=list(Players=ss_means$V3), FUN=sum)

#Construct final tables for plotting
final_type<- cbind(FC_m, HT_m$x, Der_m$x, srna_m$x, ss_m$x)
names(final_type)<- c("RNA_Type", "featureCounts", "HTSeq", "derfinder", "srnadiff", "shortstack")
mlt<- melt(final_type)
names(mlt) <- c("RNA_Type", "Tools", "Counts")
mlt$RNA_Type<- factor(mlt$RNA_Type, levels = c("protein_coding", "miRNA", "misc_RNA", "lincRNA", "antisense", "pseudogene", "processed_transcript", "rRNA", "snoRNA", "snRNA","other"))


ggplot(mlt, aes(x=RNA_Type, y = Counts, fill = Tools)) + geom_bar(position="dodge", stat = "identity") + Eric_theme + theme(axis.text.x=element_text(angle = 45, hjust=1))






#Only significant results

combi <- matrix(ncol = 5, nrow = length(all_DE_genes))
rownames(combi) <- all_DE_genes        
colnames(combi) <- c("FC", "HT", "Der", "srna", "sstack")
combi[,1] <- FC_DGE$FDR[match(rownames(combi), rownames(FC_DGE))]
combi[,2] <- HT_DGE$FDR[match(rownames(combi), rownames(HT_DGE))]
combi[,3] <- Der_DGE$FDR[match(rownames(combi), rownames(Der_DGE))]
combi[,4] <- srna_DGE$FDR[match(rownames(combi), rownames(srna_DGE))]
combi[,5] <- ss_DGE$FDR[match(rownames(combi), rownames(ss_DGE))]


FC_sig <-combi[,1]
x<- FC_sig[FC_sig > 0.05]
FC_sig <- FC_sig[!(FC_sig %in% x)]

HT_sig <- combi[,2]
x<- HT_sig[HT_sig > 0.05]
HT_sig <- HT_sig[!(HT_sig %in% x)]

Der_sig <- combi[,3]
x<- Der_sig[Der_sig > 0.05]
Der_sig <- Der_sig[!(Der_sig %in% x)]

srna_sig <- combi[,4]
x<- srna_sig[srna_sig > 0.05]
srna_sig <- srna_sig[!(srna_sig %in% x)]

ss_sig <- combi[,5]
x<- ss_sig[ss_sig > 0.05]
ss_sig <- ss_sig[!(ss_sig %in% x)]


venn.diagram(
  x = list(names(FC_sig), names(HT_sig), names(Der_sig), names(srna_sig), names(ss_sig)),
  category.names = c("featureCounts" , "HTSeq" , "derfinder", "srnadiff", "shortstack"),
  filename = 'venn_diagramm_als_ss_sig.png',
  output=TRUE
)




#Unannotated analysis. 

der_ALS <- read.table("D:/ALS_project/Analysis20200630/regions_derfinder_20210531.txt", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
rownames(der_ALS)<- paste(der_ALS$contig, ":", der_ALS$start, "-", der_ALS$end, sep = "")
der_ALS<- der_ALS[4:17]
der_ALS <- der_ALS[rowMeans(der_ALS) > 20,]
der_ALS_DGE<- DGE_EdgeR(der_ALS)
names<- colsplit(rownames(der_ALS_DGE), ":", names=c("col1","col2"))
names1<- colsplit(names$col2, "-", names=c("col1","col2"))
der_ALS_DGE$chr<- names$col1
der_ALS_DGE$start<- names1$col1
der_ALS_DGE$end<- names1$col2
der_df<- makeGRangesFromDataFrame(der_ALS_DGE, keep.extra.columns = T)
der_df1<- der_df[der_df$FDR < 0.05,]


##shortstack##
ss<- read.table(file.choose(),header = T, sep = "\t", row.names = 1)

ssplit1<- colsplit(rownames(ss), ":", names=c("chr","start"))
ssplit2<- colsplit(ssplit1$start, "-", names=c("start","end"))
ss <- cbind(ssplit1$chr, ssplit2, ss[,3:16])
col.order <- c("ssplit1$chr", "start", "end", "ACRI40", "ACRI42", "ACRI62", "ALS16", "ALS1", "ALS2", "ALS3_cutadapt", "ALS4", "ALS5", "ALS6", "ALS7", "Control_orange", "Control_purple", "Control_red")
ss<- ss[,col.order]

#find overlapping regions with annotation
ss_gr<- makeGRangesFromDataFrame(ss, keep.extra.columns = T, seqnames.field = "ssplit1$chr")
mt_dfr<-cbind(as.character(seqnames(ss_gr)),start(ss_gr),end(ss_gr))
colnames(mt_dfr)<-c("contig","start","end")
mt_dfr<-cbind(mt_dfr,values(ss_gr))
mt_dfr$contig<-as.character(mt_dfr$contig)
mt_dfr$start<-as.character(mt_dfr$start)
mt_dfr$end<-as.character(mt_dfr$end)

#DGE analysis
sstack<- mt_dfr
rownames(sstack) <- paste(sstack$contig, ":", sstack$start,"-", sstack$end, sep = "")
sstack_counts <- sstack[,4:17]
sstack_counts <- as.matrix(sstack_counts)
sstack_counts <- sstack_counts[rowMeans(sstack_counts) > 20,]
sstack_counts <- as.data.frame(sstack_counts)
ss_DGE <- DGE_EdgeR(sstack_counts)
ssplit1<- colsplit(rownames(ss_DGE), ":", names=c("chr","start"))
ssplit2<- colsplit(ssplit1$start, "-", names=c("start","end"))
ss_DGE <- cbind(ssplit1$chr, ssplit2, ss_DGE[,1:5])
ss_regions<- makeGRangesFromDataFrame(ss_DGE, keep.extra.columns = T, seqnames.field = "ssplit1$chr")
ss_regions <- ss_regions[!(seqnames(ss_regions) %in% c("chrX", "chrY", "chrM")),]
ss_regions1 <- ss_regions[ss_regions$FDR < 0.05,]


####Making a new annotation for srnadiff#####
### Overlap of low-count-filtered derfinder and shortstack granges used for running srnadiff in annotation mode ###
new_annot<- subsetByOverlaps.keepAllMeta(der_df, ss_regions)
new_annot<- c(der_df, ss_regions)
new_annot <- unique(new_annot)
new_annot <- reduce(new_annot)
new_annot$type <- rep("gene", length(new_annot))
new_annot$gene_id <- paste(seqnames(new_annot),":", start(new_annot),"-", end(new_annot),sep = "")
setwd("D:/ALS_project")
export(new_annot, "derfinder_shortstack_regions.gtf")

SampleInfo<- matrix(nrow = 14, ncol = 3)
SampleInfo[,1]<-files
SampleInfo[,2]<-files
SampleInfo[,3]<-c(rep(0,3),rep(1,8), rep(0,3))
colnames(SampleInfo)<-c("FileName", "SampleName", "Condition")
SampleInfo <- as.data.frame(SampleInfo)
annotReg   <- readAnnotation("derfinder_shortstack_regions.gtf", feature="gene", source=NULL)

srnaExp    <- srnadiffExp(files, SampleInfo, annotReg)
parameters(srnaExp) <- list(minDepth=0, minSize=10, minGap = 0, minOverlap = 1)
#parameters(srnaExp) <- list(minDepth=5, minLogFC=0.1)
srnaExp    <- srnadiff(srnaExp,segMethod = "annotation")
#srnaExp    <- srnadiff(srnaExp,segMethod = "hmm") # of note, the choice of method will greatly affect the number of statistically significant res
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

regions <- regions(srnaExp, pvalue=1)

venn.diagram(
  x = list(names(der_df), names(ss_regions), rownames(srna_counts)),
  category.names = c( "derfinder", "shortstack", "srnadiff"),
  filename = 'venn_diagramm_noanno_sig.png',
  output=TRUE
)


#Get overlapping regions for downstream wet-lab targets#
#targets ID'ed by 2 tools were analyzed further. 
x<- subsetByOverlaps.keepAllMeta(der_df1, regions)
x2<- subsetByOverlaps.keepAllMeta(x, ss_regions1)
UA<-x2[!(x2 %over% gffdata),]


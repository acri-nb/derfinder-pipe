

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

#Set Working directory for Alignment files (.BAM)
# setwd("D:/ALS_project")
# files <- list.files(pattern = "\\.bam$") 
# 
# #Set Derfinder Cutoff for Quantifying expressed regions
# cutoff = 1 
# 
# #Set desired contigs for analysis
# chrs = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
# 
#chrs = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrM","chr1_KI270706v1_random","chr1_KI270707v1_random","chr1_KI270708v1_random","chr1_KI270709v1_random","chr1_KI270710v1_random","chr1_KI270711v1_random","chr1_KI270712v1_random","chr1_KI270713v1_random","chr1_KI270714v1_random","chr2_KI270715v1_random","chr2_KI270716v1_random","chr3_GL000221v1_random","chr4_GL000008v2_random","chr5_GL000208v1_random","chr9_KI270717v1_random","chr9_KI270718v1_random","chr9_KI270719v1_random","chr9_KI270720v1_random","chr11_KI270721v1_random","chr14_GL000009v2_random","chr14_GL000225v1_random","chr14_KI270722v1_random","chr14_GL000194v1_random","chr14_KI270723v1_random","chr14_KI270724v1_random","chr14_KI270725v1_random","chr14_KI270726v1_random","chr15_KI270727v1_random","chr16_KI270728v1_random","chr17_GL000205v2_random","chr17_KI270729v1_random","chr17_KI270730v1_random","chr22_KI270731v1_random","chr22_KI270732v1_random","chr22_KI270733v1_random","chr22_KI270734v1_random","chr22_KI270735v1_random","chr22_KI270736v1_random","chr22_KI270737v1_random","chr22_KI270738v1_random","chr22_KI270739v1_random","chrY_KI270740v1_random","chrUn_KI270302v1","chrUn_KI270304v1","chrUn_KI270303v1","chrUn_KI270305v1","chrUn_KI270322v1","chrUn_KI270320v1","chrUn_KI270310v1","chrUn_KI270316v1","chrUn_KI270315v1","chrUn_KI270312v1","chrUn_KI270311v1","chrUn_KI270317v1","chrUn_KI270412v1","chrUn_KI270411v1","chrUn_KI270414v1","chrUn_KI270419v1","chrUn_KI270418v1","chrUn_KI270420v1","chrUn_KI270424v1","chrUn_KI270417v1","chrUn_KI270422v1","chrUn_KI270423v1","chrUn_KI270425v1","chrUn_KI270429v1","chrUn_KI270442v1","chrUn_KI270466v1","chrUn_KI270465v1","chrUn_KI270467v1","chrUn_KI270435v1","chrUn_KI270438v1","chrUn_KI270468v1","chrUn_KI270510v1","chrUn_KI270509v1","chrUn_KI270518v1","chrUn_KI270508v1","chrUn_KI270516v1","chrUn_KI270512v1","chrUn_KI270519v1","chrUn_KI270522v1","chrUn_KI270511v1","chrUn_KI270515v1","chrUn_KI270507v1","chrUn_KI270517v1","chrUn_KI270529v1","chrUn_KI270528v1","chrUn_KI270530v1","chrUn_KI270539v1","chrUn_KI270538v1","chrUn_KI270544v1","chrUn_KI270548v1","chrUn_KI270583v1","chrUn_KI270587v1","chrUn_KI270580v1","chrUn_KI270581v1","chrUn_KI270579v1","chrUn_KI270589v1","chrUn_KI270590v1","chrUn_KI270584v1","chrUn_KI270582v1","chrUn_KI270588v1","chrUn_KI270593v1","chrUn_KI270591v1","chrUn_KI270330v1","chrUn_KI270329v1","chrUn_KI270334v1","chrUn_KI270333v1","chrUn_KI270335v1","chrUn_KI270338v1","chrUn_KI270340v1","chrUn_KI270336v1","chrUn_KI270337v1","chrUn_KI270363v1","chrUn_KI270364v1","chrUn_KI270362v1","chrUn_KI270366v1","chrUn_KI270378v1","chrUn_KI270379v1","chrUn_KI270389v1","chrUn_KI270390v1","chrUn_KI270387v1","chrUn_KI270395v1","chrUn_KI270396v1","chrUn_KI270388v1","chrUn_KI270394v1","chrUn_KI270386v1","chrUn_KI270391v1","chrUn_KI270383v1","chrUn_KI270393v1","chrUn_KI270384v1","chrUn_KI270392v1","chrUn_KI270381v1","chrUn_KI270385v1","chrUn_KI270382v1","chrUn_KI270376v1","chrUn_KI270374v1","chrUn_KI270372v1","chrUn_KI270373v1","chrUn_KI270375v1","chrUn_KI270371v1","chrUn_KI270448v1","chrUn_KI270521v1","chrUn_GL000195v1","chrUn_GL000219v1","chrUn_GL000220v1","chrUn_GL000224v1","chrUn_KI270741v1","chrUn_GL000226v1","chrUn_GL000213v1","chrUn_KI270743v1","chrUn_KI270744v1","chrUn_KI270745v1","chrUn_KI270746v1","chrUn_KI270747v1","chrUn_KI270748v1","chrUn_KI270749v1","chrUn_KI270750v1","chrUn_KI270751v1","chrUn_KI270752v1","chrUn_KI270753v1","chrUn_KI270754v1","chrUn_KI270755v1","chrUn_KI270756v1","chrUn_KI270757v1","chrUn_GL000214v1","chrUn_KI270742v1","chrUn_GL000216v2","chrUn_GL000218v1","chrEBV")
# #set Read Length
# Rlen = 22

#Function to get count matrices from derfinder output, Differentially expressed regions only with no annotations
Get_Count_Matrix<-function(chrs,cutoff, bams, Len){
  
  #Create fullCoverage object for derfinder
  fullCov <- fullCoverage(files = bams, chrs = chrs, verbose = F)
  filteredCov <- lapply(fullCov, filterData, cutoff = cutoff)
  
  
  #Get library sizes
  TM<-vector(mode="integer", length = 0)
  
  for (i in seq(1,length(bams))){
    Tmapped<-getTotalMapped(bams[[i]],chrs = chrs)
    TM<-c(TM,Tmapped[1])
  }
  
  #Get Expressed-region level counts for fullCoverage object
  regionMat <- regionMatrix(fullCov, cutoff = cutoff, L = Len, verbose = FALSE, totalMapped = TM, targetSize = mean(TM))
  
  #Extract data from RegionMatrix Object
  GRL<-GRangesList()
  
  #Create GrangesList object containing the information for each contig. 
  for (i in chrs) {
    values(regionMat[[i]]$regions)<-regionMat[[i]]$coverageMatrix
    GRL<-c(GRL,GRangesList(regionMat[[i]]$regions))
  }
  
  #Unlist this Grange list to get a Granges object. 
  UL<-unlist(GRL)
  
  #Create a dataframe object from the annotated subset
  mt_dfr<-cbind(as.character(seqnames(UL)),start(UL)-1,end(UL))
  colnames(mt_dfr)<-c("contig","start","end")
  mt_dfr<-cbind(mt_dfr,values(UL))
  
  #Split and sum the counts in the appropriate column. Only relevant for annotation-based analysis
  # summing<-function(x) sum(na.omit(sapply(strsplit(x,";"),as.integer)))
  # for (i in names(dfr)[4:length(bams)+3]){
  #   sums<-sapply(dfr[[i]], summing)
  #   #replace with new values
  #   dfr[[i]]<-sums
  # }
  mt_dfr$contig<-as.character(mt_dfr$contig)
  mt_dfr$start<-as.character(mt_dfr$start)
  mt_dfr$end<-as.character(mt_dfr$end)
  return(mt_dfr)
}


#example
#df_1<-Get_Count_Matrix(chrs,cutoff,files,Rlen)




################################## Generating Count Matrices with Annotation ###################################

#import annotation in GFF format. Analyses were done at the transcript or gene-level for comparing featurecounts, HTSeq and derfinder. 
#Headers must contain the following fields: gene_id; gene_type
#other header titles can be assigned to new variables

#Example:   gffdata$gene_id<- gffdata$ID
#           gffdata$gene_type<- gffdata$TYPE


# gffdata<-import.gff(file.choose())
# gffdata$score<-NULL
# seqlevelsStyle(gffdata) <- "UCSC"

# We recommend using a single feature-type annotation to simplify interpretation, in this case gene
# gffdata<- gffdata[gffdata$type == "gene",]


#Function to get count matrices from derfinder output with only annotated genes
Get_Classic_Annotated_Matrix<-function(chrs,cutoff, bams, Len, anno){
  
  #Create fullCoverage object for derfinder
  fullCov <- fullCoverage(files = bams, chrs = chrs, verbose = F)
  filteredCov <- lapply(fullCov, filterData, cutoff = cutoff)
  
  
  #Get library sizes
  TM<-vector(mode="integer", length = 0)

  for (i in seq(1,length(bams))){
    Tmapped<-getTotalMapped(bams[[i]],chrs = chrs)
    TM<-c(TM,Tmapped[1])
  }
  
  #Get Expressed-region level counts for fullCoverage object
  regionMat <- regionMatrix(fullCov, cutoff = cutoff, L = Len, verbose = FALSE, targetSize = mean(TM), totalMapped = TM)
  
  #Extract data from RegionMatrix Object
  GRL<-GRangesList()
  
  #Create GrangesList object containing the information for each contig. 
  for (i in chrs) {
    values(regionMat[[i]]$regions)<-regionMat[[i]]$coverageMatrix
    GRL<-c(GRL,GRangesList(regionMat[[i]]$regions))
  }
  
  #Unlist this Grange list to get an unlisted Granges object. 
  UL<-unlist(GRL)
  
  #Find the overlap between annotated regions and the results file. 
  OL<- subsetByOverlaps.keepAllMeta(anno, UL)
  
  #Find regions in the results corresponding to Unannotated regions of the genome, per our GTF file. 
  UA<-UL[!(UL %over% anno),]
  
  #Create a dataframe object from the annotated subset
  mt_dfr<-cbind(as.character(seqnames(OL)),start(OL)-1,end(OL))
  colnames(mt_dfr)<-c("contig","start","end")
  mt_dfr<-cbind(mt_dfr,values(OL))
  
  #Split and sum the counts in the appropriate column. Only relevant for annotation-based analysis
  summing<-function(x) sum(na.omit(sapply(strsplit(x,";"),as.integer))) #Value columns are not strings and need to be
  for (i in bams){
    sums<-sapply(mt_dfr[[i]], summing)
    #replace with new values
    mt_dfr[[i]]<-sums
  }
  
  mt_dfr$contig<-as.character(mt_dfr$contig)
  mt_dfr$start<-as.character(mt_dfr$start)
  mt_dfr$end<-as.character(mt_dfr$end)
  return(mt_dfr)
}

#Example of annotated Count Matrix with Gencode GFF
#df_1<-Get_Classic_Annotated_Matrix(chrs,cutoff,files,Rlen, gffdata)



################################### Generating stats and assigning to annotation #####################################

Get_Annotated_Matrix<-function(chrs,cutoff, bams, Len, anno){
  
  #Create fullCoverage object for derfinder
  print("Creating Full Coverage Object...")
  fullCov <- fullCoverage(files = bams, chrs = chrs, verbose = F)
  filteredCov <- lapply(fullCov, filterData, cutoff = cutoff)
  rm(fullCov)
  
  #Get library sizes
  TM<-vector(mode="integer", length = 0)
  
  for (i in seq(1,length(bams))){
    Tmapped<-getTotalMapped(bams[[i]],chrs = chrs)
    TM<-c(TM,Tmapped[1])
  }
  
  #Get Expressed-region level counts for fullCoverage object
  print("Extracting Count matrix from derfinder object...")
  regionMat <- list()
  regionMat <- regionMatrix(filteredCov, cutoff = cutoff, L = Len, verbose = FALSE, targetSize = mean(TM), totalMapped = TM)
  
  #Extract data from RegionMatrix Object
  GRL<-GRangesList()
  
  #Create GrangesList object containing the information for each contig. 
  for (i in chrs) {
    values(regionMat[[i]]$regions)<-regionMat[[i]]$coverageMatrix
    GRL<-c(GRL,GRangesList(regionMat[[i]]$regions))
  }
  
  #Unlist this Grange list to get an unlisted Granges object. 
  UL<-unlist(GRL)

  #find overlaps between the expressed regions and the target annotation
  hits1<-findOverlaps(UL,gffdata)
  
  #custom function to replace a vector of redundant strings with only unique values
  pasteu<-function(x){
    paste(unique(unlist(strsplit(x,";"))),collapse =';')
  }
  
  #convert to dataframe
  print("Finding Overlaps between annotation and Expressed Regions...")
  hits1<-data.frame(hits1)
  hits1<-as.data.table(hits1)
  
  #isolate gene IDs, RNA types and % coverage from the GTF annotation with associated indices
  hits1$id<-NA
  hits1$type<-NA
  hits1$perc_cov<-NA
  hits1$id<-sapply(gffdata$gene_id[(hits1$subjectHits)], `[[`, 1)
  hits1$type<-unlist(gffdata$gene_type[(hits1$subjectHits)])
  hits1$perc_cov<-round(width(UL)[hits1$queryHits]/width(gffdata)[(hits1$subjectHits)],2)
  
  #Aggregate findOverlaps hits by query indices while reducing the overlapping annotations to get only unique IDs
  hits2<-aggregate(hits1$id,list(hits1$queryHits),paste,collapse=';')
  reduced_ids<-lapply(hits2$x,pasteu)
  reduced_ids<-unlist(reduced_ids)
  hits2$x<-reduced_ids
  
  hits3<-aggregate(hits1$type,list(hits1$queryHits),paste,collapse=';')
  reduced_types<-lapply(hits3$x,pasteu)
  reduced_types<-unlist(reduced_types)
  hits3$x<-reduced_types
  
  hits4<-aggregate(hits1$perc_cov,list(hits1$queryHits),paste,collapse=';')
  reduced_cov<-lapply(hits4$x,pasteu)
  reduced_cov<-unlist(reduced_cov)
  hits4$x<-reduced_cov

  #Create a dataframe object from the annotated subset
  print("Building Final Tables")
  mt_dfr<-cbind(as.character(seqnames(UL)),start(UL)-1,end(UL))
  colnames(mt_dfr)<-c("contig","start","end")
  mt_dfr<-cbind(mt_dfr,values(UL))

  mt_dfr$contig<-as.character(mt_dfr$contig)
  mt_dfr$start<-as.character(mt_dfr$start)
  mt_dfr$end<-as.character(mt_dfr$end)

  #Append gene ids and RNA type info to expressed regions
  mt_dfr$ids<-NA
  mt_dfr$type<-NA
  mt_dfr$perc_cov<-NA
  mt_dfr$ids[hits2$Group.1]<-hits2$x
  mt_dfr$type[hits3$Group.1]<-hits3$x
  mt_dfr$perc_cov[hits4$Group.1]<-hits4$x
  
  #Return Complete Object
  return(mt_dfr)
}

#Example of annotated Count Matrix with Gencode GFF
# df_1<-Get_Annotated_Matrix(chrs,cutoff,files,Rlen, gffdata)
# write.table(df_1, file = "New_PCa_Agnost_20200225.txt", sep = "\t")










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

# DGE_EdgeR<-function(mat){
#   library("edgeR")
# 
#   design<- matrix(ncol=2,nrow=12)
#   rownames(design)<-names(Der_res)[25:36]
#   colnames(design)<-c("group","dummy")
#   design[,1]<-c(rep(0,6),rep(1,6))
# 
#   pheno <- new("AnnotatedDataFrame", data = as.data.frame(design))
#   factor1<-as.factor(pheno$group)
#   dge <- DGEList(mat)
#   dge <- calcNormFactors(dge, method = "TMM")
#   pseudo_TMM <- log2(scale(dge$counts,center=FALSE,scale=dge$samples$norm.factors)+1)
# 
# 
#   design_matrix <- model.matrix(~ factor1)
# 
# 
#   #Estimate dispertions, a necessary step to assess gene-level scatter in our data.
#   dge <- estimateGLMCommonDisp(dge, design_matrix)
#   dge <- estimateGLMTrendedDisp(dge, design_matrix)
#   dge <- estimateGLMTagwiseDisp(dge, design_matrix)
# 
#   # par(mfrow=c(1,1))
#   # plotBCV(dge, main = paste0("BCV plot"))
# 
#   fit<-glmFit(dge,design_matrix)
#   res <- glmLRT(fit, coef = 2)
# 
#   pvals<- data.frame(raw.pvalue = res$table$PValue,
#                    FDR=p.adjust(res$table$PValue, "BH"))
#   # hist(pvals$raw.pvalue, xlim=c(0,1))
# 
#   res<-topTags(res,nrow(dge$counts))
#   alpha=0.05
# 
#   #Quick check to see how many genes are flagges as DE and if there are any imbalances.
#   return(res$table)
# }
# 
# raw_counts <-as.matrix(Der_res[26:37])
# rownames(raw_counts)<-Der_res$transcript_id
# raw_counts_der <- raw_counts[rowMeans(raw_counts)>20,]
# Der_DE<-DGE_EdgeR(raw_counts_der)
# 
# FC_res<-read.table(file.choose(), header = T, stringsAsFactors = F)
# raw_counts <-as.matrix(FC_res[7:18])
# rownames(raw_counts)<-FC_res$Geneid
# raw_counts_FC <- raw_counts[rowMeans(raw_counts)>20,]
# raw_counts_FC<-raw_counts_FC[rownames(raw_counts_FC) %in% gffdata$transcript_id,]
# FC_DE<-DGE_EdgeR(raw_counts_FC)
# 
# HT_res<-read.table(file.choose(), header = F, stringsAsFactors = F)
# names(HT_res)<-c("Geneid",names(FC_res[7:18]))
# raw_counts <-as.matrix(HT_res[2:13])
# rownames(raw_counts)<-HT_res$Geneid
# raw_counts_HT <- raw_counts[rowMeans(raw_counts)>20,]
# raw_counts_HT<-raw_counts_HT[rownames(raw_counts_HT) %in% gffdata$transcript_id,]
# HT_DE<-DGE_EdgeR(raw_counts_HT)
# 
# ### ANALYSIS OF SIMULATED READS RESULTS ###
# Truth <-read.table(file.choose(), sep = "\t", stringsAsFactors = F, header = TRUE)
# FC_DE<-read.table(file.choose(), sep = "\t", stringsAsFactors = F, header = TRUE)
# HT_DE<-read.table(file.choose(), sep = "\t", stringsAsFactors = F, header = TRUE)
# Der_DE<-read.table(file.choose(), sep = "\t", stringsAsFactors = F, header = TRUE)
# 
# FC_DE$is_de<-cut(FC_DE$FDR,breaks = c(-Inf,0.05,Inf),labels =c("TRUE","FALSE"))
# HT_DE$is_de<-cut(HT_DE$FDR,breaks = c(-Inf,0.05,Inf),labels =c("TRUE","FALSE"))
# Der_DE$is_de<-cut(Der_DE$FDR,breaks = c(-Inf,0.05,Inf),labels =c("TRUE","FALSE"))
# 
# FC_DE<-FC_DE[rownames(FC_DE) %in% rownames(Der_DE),]
# HT_DE<-HT_DE[rownames(HT_DE) %in% rownames(Der_DE),]
# Der_DE<-Der_DE[rownames(Der_DE) %in% rownames(FC_DE),]
# 
# FC_DE<-FC_DE[order(rownames(FC_DE)),]
# HT_DE<-HT_DE[order(rownames(HT_DE)),]
# Der_DE<-Der_DE[order(rownames(Der_DE)),]
# 
# Truth<-Truth[Truth$transcriptid %in% rownames(Der_DE),]
# Truth<-Truth[order(Truth$transcriptid),]
# 
# All_data<-cbind(Truth$transcriptid,Truth$DEstatus.2,HT_DE$is_de,FC_DE$is_de, Der_DE$is_de)
# replace<-ifelse(All_data[,2] == "TRUE","1","2")
# colnames(All_data)<-c("id","truth","HTseq","FeatCount","DerFinder")
# All_data[,2]<-replace
# 
# 
# 
# sensitivity(table(All_data[,2],All_data[,3]))
# specificity(table(All_data[,2],All_data[,3]))
# 
# sensitivity(table(All_data[,2],All_data[,4]))
# specificity(table(All_data[,2],All_data[,4]))
# 
# sensitivity(table(All_data[,2],All_data[,5]))
# specificity(table(All_data[,2],All_data[,5]))
# 
# comp1<- cbind(Der_DE$logFC, FC_DE$logFC)
# colnames(comp1)<-c("derfinder", "featureCounts")
# comp1<- as.data.frame(comp1)
# comp2<- cbind(Der_DE$logFC, HT_DE$logFC)
# colnames(comp1)<-c("derfinder", "featureCounts")
# comp2<- as.data.frame(comp2)
# g<- ggplot(comp1, aes(x=derfinder, y=featureCounts))+ geom_point(alpha=0.1) +xlab("derfinder")+ ylab("featureCounts")+ Eric_theme
# g<- ggplot(comp2, aes(x=derfinder, y=featureCounts))+ geom_point(alpha=0.1) +xlab("derfinder")+ ylab("featureCounts")+ Eric_theme

# ####### For enrichment analysis ############################################################################################
# 
# #remove . from the end of ensembl IDs...
# 
# names<- colsplit(rownames(res_der),"[.]",names = c("col1","col2"))
# rownames(res_der)<-names$col1
# names<- colsplit(rownames(res_hts),"[.]",names = c("col1","col2"))
# rownames(res_hts)<-names$col1
# names<- colsplit(rownames(res_fc),"[.]",names = c("col1","col2"))
# rownames(res_fc)<-names$col1
# 
# Clustlist<-bitr(rownames(res_fc[res_fc$FDR<0.1,]), fromType = "ENSEMBL",toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop=TRUE)
# 
# geneList = Clustlist$ENTREZID
# geneList = sort(geneList, decreasing = TRUE)
# 
# pwe <- enrichKEGG(gene=geneList, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5)
# pwe1 <- enrichGO(gene=geneList, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5)
# head(summary(pwe))
# 
# write.table(pwe, file="KEGG_FC.txt", sep="\t")
# write.table(pwe1, file="GO_FC.txt", sep="\t")
# 
# 
# 
# ####### FOR COMPARISON WITH EPIGENETIC DATA
# 
# bed1="PMC_neuron.1.bed"
# bed2="PUT_neuron.1.bed"
# out.gr1 <- import(bed1,format = "bed")
# out.gr2 <- import(bed2,format= "bed")
# 
# subset_ATAC<- df_UA[df_UA$gene_names %in% rownames(res_der[res_der$FDR<0.05,]),]
# subset_ATAC.gr<- makeGRangesFromDataFrame(subset_H3K4, seqnames.field = "seqnames", start.field = "starts",end.field = "ends")
# dtn<- distanceToNearest(subset_ATAC.gr, out.gr1)
# dtn<- as.data.frame(dtn)
# hist(dtn$distance/1000, breaks = 1000, xlim = c(0,100))
# dtn[dtn$distance<1000,]
# subset_ATAC.gr[c(1,16,21,38),]
# 
# 
# 
# 
# 
# ############### Comparison with mirDeep2
# 
# mirdeep<- read.table(file.choose(), header = T, stringsAsFactors = F, sep = "\t")
# 
# names<- colsplit(AN_filtered$gene_id, "[.]", names=c("col1","col2"))
# AN_miR<- cbind(AN_filtered$ALS1_bowtie2.bam,AN_filtered$ALS1_bowtie2.bam)
# rownames(AN_miR)<- names$col1
# colnames(AN_miR)<-c("raw","total")
# AN_miR<- data.frame(AN_miR)
# 
# 
# #mirdeep1<-mirdeep[mirdeep$COUNT>20,]
# mirdeep1<- aggregate(mirdeep$COUNT~mirdeep$ENSG, mirdeep,FUN=max)
# mirdeep1<- mirdeep1[mirdeep1$`mirdeep$ENSG` %in% rownames(AN_miR),]
# mirdeep1<- mirdeep1[order(mirdeep1$`mirdeep$ENSG`),]
# 
# 
# subsetB<- AN_miR[rownames(AN_miR) %in% mirdeep1$`mirdeep$ENSG`,]
# subsetB<-subsetB[order(rownames(subsetB)),]
# 
# 
# plot(log(mirdeep1$`mirdeep$COUNT`),log(subsetB$raw))
# cor.test(log(mirdeep1$`mirdeep$COUNT`),log(subsetB$raw), method = "spearman")







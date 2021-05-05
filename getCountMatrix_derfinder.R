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


setwd("D:/ALS_project")
files <- list.files(pattern = "\\.bam$") 

#Set Derfinder Cutoff for Quantifying expressed regions
cutoff = 1

#Set desired contigs for analysis (autosomes in this case)
chrs = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")

#set median Read Length
Rlen = 22

################################## Generating Count Matrices with Annotation ###################################

#import annotation in GFF format. Analyses were done at the transcript or gene-level for comparing featurecounts, HTSeq and derfinder. 
#Headers must contain the following fields: gene_id; gene_type
#other header titles can be assigned to new variables

#Example:   gffdata$gene_id<- gffdata$ID
#           gffdata$gene_type<- gffdata$TYPE


gffdata<-import.gff("D:/ALS_project/Analysis20200109/gencode.v19.annotation.gtf")
gffdata$score<-NULL
seqlevelsStyle(gffdata) <- "UCSC"

# We used a single feature-type annotation to simplify benchmarking, in this case gene
gffdata<- gffdata[gffdata$type == "gene",]

#Keep a list of contigs and genes to exclude from downstream analyses (useful later, X, Y and M)
gff_sub <- gffdata[seqnames(gffdata) %in% c("chrX", "chrY", "chrM"),]
supp_genes <- unique(gff_sub$gene_id)

#Function to get count matrices from derfinder output with only annotated genes (gene-level counts)
Get_Classic_Gene_Exprs_Matrix<-function(chrs,cutoff, bams, Len, anno){
  
  #Create fullCoverage object for derfinder
  fullCov <- fullCoverage(files = bams, chrs = chrs, verbose = F)
  #filteredCov <- lapply(fullCov, filterData, cutoff = cutoff)
  
  
  #Get library sizes
  TM<-vector(mode="integer", length = 0)
  
  for (i in seq(1,length(bams))){
    Tmapped<-getTotalMapped(bams[[i]],chrs = chrs)
    TM<-c(TM,Tmapped[1])
  }
  
  #Get Expressed-region level counts for fullCoverage object
  regionMat <- list()
  regionMat <- regionMatrix(fullCov, cutoff = cutoff, L = Len, verbose = FALSE)
  
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
  #UA<-UL[!(UL %over% anno),]
  
  #Create a dataframe object from the annotated subset
  mt_dfr<-cbind(as.character(seqnames(OL)),start(OL)-1,end(OL))
  colnames(mt_dfr)<-c("contig","start","end")
  mt_dfr<-cbind(mt_dfr,values(OL))
  print(head(mt_dfr))
  #Split and sum the counts in the appropriate column. Only relevant for annotation-based analysis
  summing<-function(x) sum(na.omit(sapply(strsplit(x,";"),as.integer))) 
  for (i in bams){
    sums<-sapply(mt_dfr[[i]], summing)
    #replace with new values
    mt_dfr[[i]]<-sums
  }
  print(head(mt_dfr))
  mt_dfr$contig<-as.character(mt_dfr$contig)
  mt_dfr$start<-as.character(mt_dfr$start)
  mt_dfr$end<-as.character(mt_dfr$end)
  return(mt_dfr)
}


results1<- Get_Classic_Gene_Exprs_Matrix(chrs,cutoff,files,Rlen, gffdata)


################################## Generating Count Matrices with derfinder (region-level) ###################################


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
  
  #Get Expressed-region level counts for fullCoverage object, use total mapped and tagetsize only if you want to normalize, but we do TMM later, so raw is ok
  print("Extracting Count matrix from derfinder object...")
  regionMat <- list()
  regionMat <- regionMatrix(filteredCov, cutoff = cutoff, L = Len, verbose = FALSE)#, targetSize = mean(TM), totalMapped = TM) 
  
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
  
  #convert to dataframe (and datatable for speed)
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


results2<-Get_Annotated_Matrix(chrs,cutoff,files,Rlen, gffdata)






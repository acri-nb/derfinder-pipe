# derfinder-pipe
R functions for the annotation-agnostic approach.


### Two main functions are contained in the derfinder Rscript 
The Second - Get_Classic_Gene_Exprs_Matrix() is a function which uses derfinder to summarize counts from .bam files, however the expressed regions 
within a larger feature (ex. gene) are summed and counted as one feature, rather than mutliple smaller expressed regions. Expressed
regions outside known annotations are discarded, which returns a count matrix which ressembles classic count matrices produces by 
tools such as featureCounts. The output is in the following format. 

| Contig | Start | End | Sample1 | Sample2 | Sample3 | Sample4 | ID | Type |
| ------ | ----- | --- | ------- | ------- | ------- | ------- | -- | ---- |
| chr1 | 10 | 20 | 5 | 6 | 2 | 3 | ENSG1 | lncRNA |
| chr2 | 50 | 75 | 3 | 3 | 9 | 8 | ENSG2 | protein_coding |

Other tools usually directly output a similar count matrix or work with an input gff/gtf file for generating matrices.  

The second (and currently used in our lab) - Get_Annotated_Matrix() is a function which reports expressed regions rather than genes, much like Get_Count_Matrix(), however
Annotations are added post-hoc to each region, detailing all known annotations overlapping an expressed regions. 
This allows the user to see all annotations in a region of interest, while keeping all novel expressed RNAs. 
The format is as follows. 

| Contig | Start | End | Sample1 | Sample2 | Sample3 | Sample4 | ID | Type |
| ------ | ----- | --- | ------- | ------- | ------- | ------- | -- | ---- |
| chr4 | 11 | 20 | 5 | 6 | 2 | 3 | ENSG1;ENSG10 | lncRNA,protein_coding |
| chr5 | 12 | 79 | 1 | 1 | 9 | 9 | NA | Unknown |


Some supplementary code is commented out and is related to our initial comparisons with other software.

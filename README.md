# derfinder-pipe
R functions for the annotation-agnostic approach.


### Three main functions are contained in the Rscript 
The first - Get_count_matrix() is a 'quick and dirty' method which applies the derfinder 'regionMatrix' function to obtain 
Only expressed Regions from the input .bam files. This function simply takes bam files as input and outputs a dataframe in the
following format. 

| Contig | Start | End | Sample1 | Sample2 | Sample3 | Sample4 |
| ------ | ----- | --- | ------- | ------- | ------- | ------- |
| chr1 | 10 | 20 | 5 | 6 | 2 | 3 |
| chr3 | 50 | 75 | 3 | 3 | 9 | 8 |

The Second - Get_Classic_Annotated_Matrix() is a function which uses derfinder to summarize counts, however the expressed regions 
within a larger feature (ex. gene) are summed and counted as one feature, rather than mutliple smaller expressed regions. Expressed
regions outside known annotations are discarded, which returns a count matrix which ressembles classic count matrices produces by 
tools such as featureCounts. The output is in the following format. 

| Contig | Start | End | Sample1 | Sample2 | Sample3 | Sample4 | ID | Type |
| ------ | ----- | --- | ------- | ------- | ------- | ------- | -- | ---- |
| chr1 | 10 | 20 | 5 | 6 | 2 | 3 | ENSG1 | lncRNA |
| chr2 | 50 | 75 | 3 | 3 | 9 | 8 | ENSG2 | protein_coding |

The third - Get_Annotated_Matrix() is a function which reports expressed regions rather than genes, much like Get_Count_Matrix(), however
Annotations are added post-hoc to each region, detailing all known annotations overlapping an expressed regions. 
This allows the user to see all annotations in a region of interest, while keeping all novel expressed RNAs. 
The format is as follows. 

| Contig | Start | End | Sample1 | Sample2 | Sample3 | Sample4 | ID | Type |
| ------ | ----- | --- | ------- | ------- | ------- | ------- | -- | ---- |
| chr9 | 11 | 20 | 5 | 6 | 2 | 3 | ENSG1;ENSG10 | lncRNA,protein_coding |
| chr4 | 12 | 79 | 1 | 1 | 9 | 9 | NA | Unknown |


Some supplementary code is commented out and is related to our initial comparisons with other software.

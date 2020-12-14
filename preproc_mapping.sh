##pipeline developed for workstation 10.253.0.20, but it can be addapted in other machines
## recomended usage of fastqc to check qualities of fastqs before and after mapping
##the end of this pipeline you shoud have a bam file sorted and a index of the bam file to use as input of DERFINDER pipeline

### Trimming
#trilink ion proton recommended
~/bin/trim_galore/Trim_galoreV0.6.5/trim_galore -a TCACCGACTGCCATAGAG file.fastq.gz
#nextflex triming recommended 
~/bin/trim_galore/Trim_galoreV0.6.5/trim_galore -a TGGAATTCTCGGGTGCCAAGG --clip_R1 4 --three_prime_clip_R1 4 --length 23 file.fastq.gz
#comboseq triming recommended 
~/bin/trim_galore/Trim_galoreV0.6.5/trim_galore -a AAAAAAAA --length 15 file.fastq.gz

###STAR Mapping

##reference genome in: acri@192.168.0.200:/mnt/user/LTS/backup2/Reference_data/References_genome/
##before this step recommended the use of :
##sshfs acri@192.168.0.200:/mnt/user/LTS/backup2/Reference_data/References_genome/ /References_data/References_genome/acri/
#star for trilink and netflex mapping 
/home/iarc/bin/STAR/STAR-2.7.0f/bin/Linux_x86_64/STAR   --runThreadN 20   --genomeDir /References_data/References_genome/acri/Homo_sapiens/hg38/STAR/   --readFilesIn Neo-02-3K_trimmed.fq      --outFileNamePrefix Neo-02-3K   --outFilterScoreMinOverLread 0   --outFilterMatchNmin 16   --outFilterMatchNminOverLread 0   --outFilterMismatchNoverLmax 0.05   --alignIntronMax 1

#star for comboflex
/home/iarc/bin/STAR/STAR-2.7.0f/bin/Linux_x86_64/STAR   --runThreadN 20   --genomeDir /References_data/References_genome/acri/Homo_sapiens/hg38/STAR/   --readFilesIn Neo-02-3K_trimmed.fq      --outFileNamePrefix Neo-02-3K   --outFilterScoreMinOverLread 0   --outFilterMatchNmin 15   --outFilterMatchNminOverLread 0.9  --outFilterMismatchNmax 1 --outFilterMismatchNoverLmax 0.3   --alignIntronMax 1


###Samtools sorting of Star alignment  and indexing
samtools sort -o output.sorted.bam Aligned_Star.sam; samtools index output.sorted.bam



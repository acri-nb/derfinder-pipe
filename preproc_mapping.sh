

###STAR Mapping
#reference genome in: acri@192.168.0.200:/mnt/user/LTS/backup2/Reference_data/References_genome/
#before this step recommended the use of :
#sshfs acri@192.168.0.200:/mnt/user/LTS/backup2/Reference_data/References_genome/ /References_data/References_genome/acri/
/home/iarc/bin/STAR/STAR-2.7.0f/bin/Linux_x86_64/STAR   --runThreadN 20   --genomeDir /References_data/References_genome/acri/Homo_sapiens/hg38/STAR/   --readFilesIn Neo-02-3K_trimmed.fq      --outFileNamePrefix Neo-02-3K   --outFilterScoreMinOverLread 0   --outFilterMatchNmin 16   --outFilterMatchNminOverLread 0   --outFilterMismatchNoverLmax 0.05   --alignIntronMax 1

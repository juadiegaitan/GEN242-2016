###################################################
## Exercises for Common Bioinformatics Use Cases ##
###################################################

## Download Halobacterium proteome and inspect it
module load ncbi-blast/2.2.26
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/Halobacterium_salinarum/representative/GCA_000069025.1_ASM6902v1/GCA_000069025.1_ASM6902v1_protein.faa.gz
gunzip GCA_000069025.1_ASM6902v1_protein.faa.gz
mv GCA_000069025.1_ASM6902v1_protein.faa halobacterium.faa

less halobacterium.faa # press q to quit

## How many protein sequences are stored in the downloaded file?
grep '>' halobacterium.faa | wc
grep '^>' halobacterium.faa --count

## How many proteins contain the pattern "WxHxxH" or "WxHxxHH"?
egrep 'W.H..H{1,2}' halobacterium.faa --count

## Use less to find IDs for pattern matches or use awk
awk --posix -v RS='>' '/W.H..(H){1,2}/ { print ">" $0;}' halobacterium.faa | less
awk --posix -v RS='>' '/W.H..(H){1,2}/ { print ">" $0;}' halobacterium.faa | grep '^>' | cut -c 2- | cut -f 1 -d\ > myIDs

## Create a BLASTable database with formatdb
formatdb -i halobacterium.faa -p T -o

## Query BLASTable database by IDs stored in a file (e.g. myIDs)
fastacmd -d halobacterium.faa -i myIDs > myseq.fasta

## Run BLAST search for sequences stored in myseq.fasta
blastall -p blastp -i myseq.fasta -d halobacterium.faa -o blastp.out -e 1e-6 -v 10 -b 10
blastall -p blastp -i myseq.fasta -d halobacterium.faa -m 8 -e 1e-6 > blastp.tab

## More exercises in Linux Manual 
##  http://manuals.bioinformatics.ucr.edu/home/linux-basics#TOC-Exercises

## Submit job to queuing system of cluster

## (i) Add these two lines to beginning of your script
    ## #!/bin/bash
    ## cd $PBS_O_WORKDIR
    
## (ii) Submit script to cluster as follows
    ## qsub -l nodes=1:ppn=1,mem=1gb,walltime=4:00:00 Linux.sh

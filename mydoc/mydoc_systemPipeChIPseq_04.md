---
title: Alignments
keywords: 
last_updated: Tue May  3 13:38:05 2016
---

## Read mapping with `Bowtie2` 

The NGS reads of this project will be aligned with `Bowtie2` against the
reference genome sequence (Langmead et al., 2012). The parameter settings of the
aligner are defined in the `bowtieSE.param` file. In ChIP-Seq experiments it is
usually more appropriate to eliminate reads mapping to multiple locations. To
achieve this, users want to remove the argument setting `-k 50 –non-deterministic` 
in the `bowtieSE.param` file.


{% highlight r %}
args <- systemArgs(sysma="param/bowtieSE.param", mytargets="targets_chip_trim.txt")
sysargs(args)[1] # Command-line parameters for first FASTQ file
moduleload(modules(args)) # Skip if a module system is not used
system("bowtie2-build ./data/tair10.fasta ./data/tair10.fasta") # Indexes reference genome
runCommandline(args)
writeTargetsout(x=args, file="targets_bam.txt", overwrite=TRUE)
{% endhighlight %}

Check whether all BAM files have been created

{% highlight r %}
file.exists(outpaths(args))
{% endhighlight %}

## Read and alignment stats

The following provides an overview of the number of reads in each sample
and how many of them aligned to the reference.


{% highlight r %}
read_statsDF <- alignStats(args=args) 
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
read.delim("results/alignStats.xls")
{% endhighlight %}

## Create symbolic links for viewing BAM files in IGV

The `symLink2bam` function creates symbolic links to view the BAM alignment files in a
genome browser such as IGV without moving these large files to a local
system. The corresponding URLs are written to a file with a path
specified under `urlfile`, here `IGVurl.txt`.


{% highlight r %}
symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"), 
            urlbase="http://biocluster.ucr.edu/~tgirke/", 
            urlfile="./results/IGVurl.txt")
{% endhighlight %}


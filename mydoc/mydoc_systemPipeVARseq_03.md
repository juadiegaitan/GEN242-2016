---
title: Read preprocessing
keywords: 
last_updated: Thu May  5 11:37:34 2016
---

## Read quality filtering and trimming

The following removes reads with low quality base calls (here Phred
scores below 20) from all FASTQ files.


{% highlight r %}
args <- systemArgs(sysma="param/trimPE.param", mytargets="targetsPE.txt")[1:4] # Note: subsetting!
filterFct <- function(fq, cutoff=20, Nexceptions=0) {
    qcount <- rowSums(as(quality(fq), "matrix") <= cutoff)
    fq[qcount <= Nexceptions] # Retains reads where Phred scores are >= cutoff with N exceptions
}
preprocessReads(args=args, Fct="filterFct(fq, cutoff=20, Nexceptions=0)", batchsize=100000)
writeTargetsout(x=args, file="targets_PEtrim.txt", overwrite=TRUE)
{% endhighlight %}


## FASTQ quality report

The following `seeFastq` and `seeFastqPlot` functions generate and plot a series of 
useful quality statistics for a set of FASTQ files including per cycle quality box
plots, base proportions, base-level quality trends, relative k-mer
diversity, length and occurrence distribution of reads, number of reads
above quality cutoffs and mean quality distribution. The results are
written to a PDF file named `fastqReport.pdf`.


{% highlight r %}
args <- systemArgs(sysma="param/tophat.param", mytargets="targets.txt")
fqlist <- seeFastq(fastq=infile1(args), batchsize=100000, klength=8)
pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
seeFastqPlot(fqlist)
dev.off()
{% endhighlight %}

![](systemPipeVARseq_files/fastqReport.png)
<div align="center">Figure 1: FASTQ quality report for 18 samples</div>


---
title: Alignments
keywords: 
last_updated: Sun May  1 20:44:33 2016
---

## Read mapping with `BWA-MEM` 

The NGS reads of this project are aligned against the reference genome
sequence using the highly variant tolerant short read aligner `BWA-MEM`
(Li , 2013; Li et al., 2009). The parameter settings of the aligner are
defined in the `bwa.param` file.


{% highlight r %}
args <- systemArgs(sysma="bwa.param", mytargets="targets.txt")
sysargs(args)[1] # Command-line parameters for first FASTQ file
{% endhighlight %}


Runs the alignments sequentially (_e.g._ on a single machine)


{% highlight r %}
bampaths <- runCommandline(args=args)
{% endhighlight %}

Alternatively, the alignment jobs can be submitted to a compute cluster,
here using 72 CPU cores (18 `qsub` processes each with 4 CPU cores).


{% highlight r %}
moduleload(modules(args))
system("bwa index -a bwtsw ./data/tair10.fasta")
resources <- list(walltime="20:00:00", nodes=paste0("1:ppn=", cores(args)), memory="10gb")
reg <- clusterRun(args, conffile=".BatchJobs.R", template="torque.tmpl", Njobs=18, runid="01", 
                  resourceList=resources)
waitForJobs(reg)
{% endhighlight %}

Check whether all BAM files have been created


{% highlight r %}
file.exists(outpaths(args))
{% endhighlight %}

## Read mapping with `gsnap` 

An alternative variant tolerant aligner is `gsnap` from the `gmapR` package
(Wu et al., 2010). The following code shows how to run this aligner on
multiple nodes of a computer cluster that uses Torque as scheduler.


{% highlight r %}
library(gmapR); library(BiocParallel); library(BatchJobs)
gmapGenome <- GmapGenome(reference(args), directory="data", name="gmap_tair10chr", create=TRUE)
args <- systemArgs(sysma="gsnap.param", mytargets="targetsPE.txt")
f <- function(x) {
    library(gmapR); library(systemPipeR)
    args <- systemArgs(sysma="gsnap.param", mytargets="targetsPE.txt")
    gmapGenome <- GmapGenome(reference(args), directory="data", name="gmap_tair10chr", create=FALSE)
    p <- GsnapParam(genome=gmapGenome, unique_only=TRUE, molecule="DNA", max_mismatches=3)
    o <- gsnap(input_a=infile1(args)[x], input_b=infile2(args)[x], params=p, output=outfile1(args)[x])
}
funs <- makeClusterFunctionsTorque("torque.tmpl")
param <- BatchJobsParam(length(args), resources=list(walltime="20:00:00", nodes="1:ppn=1", memory="6gb"), cluster.functions=funs)
register(param)
d <- bplapply(seq(along=args), f)
writeTargetsout(x=args, file="targets_gsnap_bam.txt")
{% endhighlight %}


## Read and alignment stats

The following generates a summary table of the number of reads in each
sample and how many of them aligned to the reference.


{% highlight r %}
read_statsDF <- alignStats(args=args) 
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
{% endhighlight %}


## Create symbolic links for viewing BAM files in IGV

The `symLink2bam` function creates symbolic links to view the BAM alignment files in a
genome browser such as IGV. The corresponding URLs are written to a file
with a path specified under `urlfile`, here `IGVurl.txt`.


{% highlight r %}
symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"), 
            urlbase="http://biocluster.ucr.edu/~tgirke/", 
            urlfile="./results/IGVurl.txt")
{% endhighlight %}



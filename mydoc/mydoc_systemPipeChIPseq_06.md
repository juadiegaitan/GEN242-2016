---
title: Annotate peaks with genomic context
keywords: 
last_updated: Sun May  1 17:34:29 2016
---

## Annotation with `ChIPpeakAnno` package

The following annotates the identified peaks with genomic context
information using the `ChIPpeakAnno` and `ChIPseeker` packages, respectively
(Zhu et al., 2010; Yu et al., 2015).


{% highlight r %}
library(ChIPpeakAnno); library(GenomicFeatures)
args <- systemArgs(sysma="param/annotate_peaks.param", mytargets="targets_macs.txt")
txdb <- loadDb("./data/tair10.sqlite")
ge <- genes(txdb, columns=c("tx_name", "gene_id", "tx_type")) 
for(i in seq(along=args)) {
    peaksGR <- as(read.delim(infile1(args)[i], comment="#"), "GRanges")
    annotatedPeak <- annotatePeakInBatch(peaksGR, AnnotationData=genes(txdb))
    df <- data.frame(as.data.frame(annotatedPeak), as.data.frame(values(ge[values(annotatedPeak)$feature,])))
    write.table(df, outpaths(args[i]), quote=FALSE, row.names=FALSE, sep="\t")
}
writeTargetsout(x=args, file="targets_peakanno.txt", overwrite=TRUE)
{% endhighlight %}




The peak annotation results are written for each peak set to separate
files in the `results` directory. They are named after the corresponding peak
files with extensions specified in the `annotate_peaks.param` file, 
here `*.peaks.annotated.xls`.

## Annotation with `ChIPseeker` package

Same as in previous step but using the `ChIPseeker` package for annotating the peaks.


{% highlight r %}
library(ChIPseeker)
txdb <- loadDb("./data/tair10.sqlite")
for(i in seq(along=args)) {
    peakAnno <- annotatePeak(infile1(args)[i], TxDb=txdb, verbose=FALSE)
    df <- as.data.frame(peakAnno)
    write.table(df, outpaths(args[i]), quote=FALSE, row.names=FALSE, sep="\t")
}
writeTargetsout(x=args, file="targets_peakanno.txt", overwrite=TRUE)
{% endhighlight %}

Summary plots provided by the `ChIPseeker` package. Here applied only to one sample
for demonstration purposes.


{% highlight r %}
peak <- readPeakFile(infile1(args)[1])
covplot(peak, weightCol="X.log10.pvalue.")
peakHeatmap(outpaths(args)[1], TxDb=txdb, upstream=1000, downstream=1000, color="red")
plotAvgProf2(outpaths(args)[1], TxDb=txdb, upstream=1000, downstream=1000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
{% endhighlight %}


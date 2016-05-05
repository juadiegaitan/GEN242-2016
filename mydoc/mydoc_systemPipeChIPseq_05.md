---
title: Utilities for coverage data
keywords: 
last_updated: Wed May  4 17:24:19 2016
---

The following introduces several utilities useful for ChIP-Seq data. They are not part of the actual
workflow.

## Rle object stores coverage information

{% highlight r %}
library(rtracklayer); library(GenomicRanges); library(Rsamtools); library(GenomicAlignments)
aligns <- readGAlignments(outpaths(args)[1])
cov <- coverage(aligns)
cov
{% endhighlight %}

## Resizing aligned reads

{% highlight r %}
trim(resize(as(aligns, "GRanges"), width = 200))
{% endhighlight %}

## Naive peak calling

{% highlight r %}
islands <- slice(cov, lower = 15)
islands[[1]]
{% endhighlight %}

## Plot coverage for defined region

{% highlight r %}
library(ggbio)
myloc <- c("Chr1", 1, 100000)
ga <- readGAlignments(outpaths(args)[1], use.names=TRUE, param=ScanBamParam(which=GRanges(myloc[1], IRanges(as.numeric(myloc[2]), as.numeric(myloc[3])))))
autoplot(ga, aes(color = strand, fill = strand), facets = strand ~ seqnames, stat = "coverage")
{% endhighlight %}


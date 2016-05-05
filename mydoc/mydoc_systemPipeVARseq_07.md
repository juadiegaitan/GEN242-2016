---
title: Annotate filtered variants
keywords: 
last_updated: Thu May  5 11:34:28 2016
---

The function `variantReport` generates a variant report using
utilities provided by the `VariantAnnotation` package. The report for
each sample is written to a tabular file containing genomic context annotations
(_e.g._ coding or non-coding SNPs, amino acid changes, IDs of affected
genes, etc.) along with confidence statistics for each variant. The parameter
file `annotate_vars.param` defines the paths to the input and output
files which are stored in a new `SYSargs` instance. 

## Basics of annotating variants

Variants overlapping with common annotation features can be identified with `locateVariants`.

{% highlight r %}
library("GenomicFeatures")
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_gatk_filtered.txt")
txdb <- loadDb("./data/tair10.sqlite")
vcf <- readVcf(infile1(args)[1], "A. thaliana")
locateVariants(vcf, txdb, CodingVariants())
{% endhighlight %}
Synonymous/non-synonymous variants of coding sequences are computed by the predictCoding function for variants overlapping with coding regions.


{% highlight r %}
fa <- FaFile(systemPipeR::reference(args))
predictCoding(vcf, txdb, seqSource=fa)
{% endhighlight %}

## Annotate filtered variants called by `GATK`


{% highlight r %}
library("GenomicFeatures")
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_gatk_filtered.txt")
txdb <- loadDb("./data/tair10.sqlite")
fa <- FaFile(systemPipeR::reference(args))
suppressAll(variantReport(args=args, txdb=txdb, fa=fa, organism="A. thaliana"))
{% endhighlight %}

## Annotate filtered variants called by `BCFtools`


{% highlight r %}
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_sambcf_filtered.txt")
txdb <- loadDb("./data/tair10.sqlite")
fa <- FaFile(systemPipeR::reference(args))
suppressAll(variantReport(args=args, txdb=txdb, fa=fa, organism="A. thaliana"))
{% endhighlight %}

## Annotate filtered variants called by `VariantTools`


{% highlight r %}
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_vartools_filtered.txt")
txdb <- loadDb("./data/tair10.sqlite")
fa <- FaFile(systemPipeR::reference(args))
suppressAll(variantReport(args=args, txdb=txdb, fa=fa, organism="A. thaliana"))
{% endhighlight %}

View annotation result for single sample

{% highlight r %}
read.delim(outpaths(args)[1])[38:40,]
{% endhighlight %}


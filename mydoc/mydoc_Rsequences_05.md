---
title: Range Operations  
keywords: 
last_updated: Wed Apr 20 19:57:09 2016
---

## Important Data Objects for Range Operations

* `IRanges`: stores range data only (IRanges library)
* `GRanges`: stores ranges and annotations (GenomicRanges library)
* `GRangesList`: list version of GRanges container (GenomicRanges library)

## Range Data Are Stored in `IRanges` and `GRanges` Containers

### Construct `GRanges` Object


{% highlight r %}
library(GenomicRanges); library(rtracklayer)
gr <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)), ranges = IRanges(1:10, end = 7:16, names = head(letters, 10)), strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)), score = 1:10, GC = seq(1, 0, length = 10)) # Example of creating a GRanges object with its constructor function.
{% endhighlight %}

### Import GFF into `GRanges` Object

{% highlight r %}
gff <- import.gff("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/gff3.gff") # Imports a simplified GFF3 genome annotation file.
seqlengths(gff) <- end(ranges(gff[which(values(gff)[,"type"]=="chromosome"),])) 
names(gff) <- 1:length(gff) # Assigns names to corresponding slot.
gff[1:4,]
{% endhighlight %}

{% highlight txt %}
## GRanges object with 4 ranges and 10 metadata columns:
##     seqnames           ranges strand |   source       type     score     phase                  ID
##        <Rle>        <IRanges>  <Rle> | <factor>   <factor> <numeric> <integer>         <character>
##   1     Chr1 [   1, 30427671]      + |   TAIR10 chromosome      <NA>      <NA>                Chr1
##   2     Chr1 [3631,     5899]      + |   TAIR10       gene      <NA>      <NA>           AT1G01010
##   3     Chr1 [3631,     5899]      + |   TAIR10       mRNA      <NA>      <NA>         AT1G01010.1
##   4     Chr1 [3760,     5630]      + |   TAIR10    protein      <NA>      <NA> AT1G01010.1-Protein
##            Name                Note          Parent       Index Derives_from
##     <character>     <CharacterList> <CharacterList> <character>  <character>
##   1        Chr1                                            <NA>         <NA>
##   2   AT1G01010 protein_coding_gene                        <NA>         <NA>
##   3 AT1G01010.1                           AT1G01010           1         <NA>
##   4 AT1G01010.1                                            <NA>  AT1G01010.1
##   -------
##   seqinfo: 7 sequences from an unspecified genome
{% endhighlight %}


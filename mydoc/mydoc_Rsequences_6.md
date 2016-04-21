---
title: Transcript Ranges
keywords: 
last_updated: Wed Apr 20 20:28:41 2016
---

Storing annotation ranges in `TranscriptDb` databases makes many operations more robust and convenient.

{% highlight r %}
library(GenomicFeatures)
download.file("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/gff3.gff", "data/gff3.gff")
txdb <- makeTxDbFromGFF(file="data/gff3.gff", format="gff", dataSource="TAIR", organism="Arabidopsis thaliana")
{% endhighlight %}

{% highlight txt %}
## Warning in .extract_exons_from_GRanges(cds_IDX, gr, ID, Name, Parent, feature = "cds", : The following orphan CDS were dropped (showing only the 6 first):
##   seqid start  end strand   ID              Parent Name
## 1  Chr1  3760 3913      + <NA> AT1G01010.1-Protein <NA>
## 2  Chr1  3996 4276      + <NA> AT1G01010.1-Protein <NA>
## 3  Chr1  4486 4605      + <NA> AT1G01010.1-Protein <NA>
## 4  Chr1  4706 5095      + <NA> AT1G01010.1-Protein <NA>
## 5  Chr1  5174 5326      + <NA> AT1G01010.1-Protein <NA>
## 6  Chr1  5439 5630      + <NA> AT1G01010.1-Protein <NA>
{% endhighlight %}

{% highlight txt %}
## Warning: closing unused connection 6 (data/SRR038845.fastq)
{% endhighlight %}

{% highlight r %}
saveDb(txdb, file="./data/TAIR10.sqlite")
{% endhighlight %}

{% highlight txt %}
## TxDb object:
## # Db type: TxDb
## # Supporting package: GenomicFeatures
## # Data source: TAIR
## # Organism: Arabidopsis thaliana
## # Taxonomy ID: 3702
## # miRBase build ID: NA
## # Genome: NA
## # transcript_nrow: 28
## # exon_nrow: 113
## # cds_nrow: 99
## # Db created by: GenomicFeatures package from Bioconductor
## # Creation time: 2016-04-20 20:19:14 -0700 (Wed, 20 Apr 2016)
## # GenomicFeatures version at creation time: 1.22.6
## # RSQLite version at creation time: 1.0.0
## # DBSCHEMAVERSION: 1.1
{% endhighlight %}

{% highlight r %}
txdb <- loadDb("./data/TAIR10.sqlite")
tr <- transcripts(txdb)
GRList <- transcriptsBy(txdb, by = "gene")
{% endhighlight %}

## `txdb` from BioMart

Alternative sources for creating `txdb` databases are BioMart, Bioc annotation packages, UCSC, etc. The following shows how to create a `txdb` from BioMart.

{% highlight r %}
library(GenomicFeatures); library("biomaRt")
txdb <- makeTranscriptDbFromBiomart(biomart = "plants_mart_25", dataset = "athaliana_eg_gene")
{% endhighlight %}

The following steps are useful to find out what is availble in BioMart. 

{% highlight r %}
listMarts() # Lists BioMart databases
mymart <- useMart("plants_mart_25") # Select one, here plants_mart_25
listDatasets(mymart) # List datasets available in the selected BioMart database
mymart <- useMart("plants_mart_25", dataset="athaliana_eg_gene")
listAttributes(mymart) # List available features 
getBM(attributes=c("ensembl_gene_id", "description"), mart=mymart)[1:4,]
{% endhighlight %}

## Efficient Sequence Parsing with `getSeq`

The following parses all annotation ranges provided by `GRanges` object (e.g. `gff`) from a genome sequence stored in a local file.

{% highlight r %}
gff <- gff[values(gff)$type != "chromosome"] # Remove chromosome ranges
rand <- DNAStringSet(sapply(unique(as.character(seqnames(gff))), function(x) paste(sample(c("A","T","G","C"), 200000, replace=T), collapse="")))
writeXStringSet(DNAStringSet(rand), "./data/test")
getSeq(FaFile("./data/test"), gff)
{% endhighlight %}

{% highlight txt %}
##   A DNAStringSet instance of length 442
##       width seq                                                                 names               
##   [1]  2269 TAATCCCGCCAGGGGTACCTTACAATGAGACG...CGTCGCAGCTAAACCTAGTATGTTCTATATAC Chr1
##   [2]  2269 TAATCCCGCCAGGGGTACCTTACAATGAGACG...CGTCGCAGCTAAACCTAGTATGTTCTATATAC Chr1
##   [3]  1871 TTCAACGTGAGTCCCCAAACTCCGTCGTCAAG...TCTACGAAATGTGTGCGGAGGACCATATCCGT Chr1
##   [4]   283 TAATCCCGCCAGGGGTACCTTACAATGAGACG...ATGCAGTAATTTCCGCTGTCTCTAGCCTAAAT Chr1
##   [5]   129 TAATCCCGCCAGGGGTACCTTACAATGAGACG...ACTGCTCATACCCGGCACTTTATCCCCCGCAG Chr1
##   ...   ... ...
## [438]   324 TCCCCTCCCAGCCTCCTAGAGGTGACGAGGAT...GACGGGGGGAAATGCTGTAGTCTCAAATGACC ChrM
## [439]   324 TCCCCTCCCAGCCTCCTAGAGGTGACGAGGAT...GACGGGGGGAAATGCTGTAGTCTCAAATGACC ChrM
## [440]   324 TCCCCTCCCAGCCTCCTAGAGGTGACGAGGAT...GACGGGGGGAAATGCTGTAGTCTCAAATGACC ChrM
## [441]   324 TCCCCTCCCAGCCTCCTAGAGGTGACGAGGAT...GACGGGGGGAAATGCTGTAGTCTCAAATGACC ChrM
## [442]   324 TCCCCTCCCAGCCTCCTAGAGGTGACGAGGAT...GACGGGGGGAAATGCTGTAGTCTCAAATGACC ChrM
{% endhighlight %}


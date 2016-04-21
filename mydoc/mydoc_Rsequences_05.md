---
title: NGS Sequences
keywords: 
last_updated: Thu Apr 21 12:28:13 2016
---

## Sequence and Quality Data: FASTQ Format

Four lines per sequence:

1. ID
2. Sequence
3. ID
4. Base call qualities (Phred scores) as ASCII characters

The following gives an example of 3 Illumina reads in a FASTQ file. The numbers
at the beginning of each line are not part of the FASTQ format. They have been added 
solely for illustration purposes.

{% highlight txt %}
1. @SRR038845.3 HWI-EAS038:6:1:0:1938 length=36
2. CAACGAGTTCACACCTTGGCCGACAGGCCCGGGTAA
3. +SRR038845.3 HWI-EAS038:6:1:0:1938 length=36
4. BA@7>B=>:>>7@7@>>9=BAA?;>52;>:9=8.=A
1. @SRR038845.41 HWI-EAS038:6:1:0:1474 length=36
2. CCAATGATTTTTTTCCGTGTTTCAGAATACGGTTAA
3. +SRR038845.41 HWI-EAS038:6:1:0:1474 length=36
4. BCCBA@BB@BBBBAB@B9B@=BABA@A:@693:@B=
1. @SRR038845.53 HWI-EAS038:6:1:1:360 length=36
2. GTTCAAAAAGAACTAAATTGTGTCAATAGAAAACTC
3. +SRR038845.53 HWI-EAS038:6:1:1:360 length=36
4. BBCBBBBBB@@BAB?BBBBCBC>BBBAA8>BBBAA@
{% endhighlight %}

## Sequence and Quality Data: `QualityScaleXStringSet`

Phred quality scores are integers from 0-50 that are
stored as ASCII characters after adding 33. The basic R functions `rawToChar` and
`charToRaw` can be used to interconvert among their representations.

Phred score interconversion

{% highlight r %}
phred <- 1:9
phreda <- paste(sapply(as.raw((phred)+33), rawToChar), collapse="")
phreda
{% endhighlight %}

{% highlight txt %}
## [1] "\"#$%&'()*"
{% endhighlight %}

{% highlight r %}
as.integer(charToRaw(phreda))-33 
{% endhighlight %}

{% highlight txt %}
## [1] 1 2 3 4 5 6 7 8 9
{% endhighlight %}

Construct `QualityScaledDNAStringSet` from scratch

{% highlight r %}
dset <- DNAStringSet(sapply(1:100, function(x) paste(sample(c("A","T","G","C"), 20, replace=T), collapse=""))) # Creates random sample sequence.
myqlist <- lapply(1:100, function(x) sample(1:40, 20, replace=T)) # Creates random Phred score list.
myqual <- sapply(myqlist, function(x) toString(PhredQuality(x))) # Converts integer scores into ASCII characters.
myqual <- PhredQuality(myqual) # Converts to a PhredQuality object.
dsetq1 <- QualityScaledDNAStringSet(dset, myqual) # Combines DNAStringSet and quality data in QualityScaledDNAStringSet object.
dsetq1[1:2]
{% endhighlight %}

{% highlight txt %}
##   A QualityScaledDNAStringSet instance containing:
## 
##   A DNAStringSet instance of length 2
##     width seq
## [1]    20 TGAGGTCTGGCCCATTGAAT
## [2]    20 ACGGCGGTAAGTCCTGCAAA
## 
##   A PhredQuality instance of length 2
##     width seq
## [1]    20 ;;A1HH/"@2E/2D"9/(DD
## [2]    20 A2+7.IG:->0*3I7"5H-%
{% endhighlight %}

## Processing FASTQ Files with ShortRead

The following expains the basic usage of `ShortReadQ` objects. To make the sample code work, 
download and unzip this [file](http://faculty.ucr.edu/~tgirke/HTML_Presentations/Manuals/Workshop_Dec_6_10_2012/Rsequences/data.zip) to your current working directory.
The following code performs the download for you.


{% highlight r %}
library(ShortRead)
download.file("http://faculty.ucr.edu/~tgirke/HTML_Presentations/Manuals/Workshop_Dec_6_10_2012/Rsequences/data.zip", "data.zip")
unzip("data.zip")
{% endhighlight %}

Important utilities for accessing FASTQ files

{% highlight r %}
fastq <- list.files("data", "*.fastq$"); fastq <- paste("data/", fastq, sep="")
names(fastq) <- paste("flowcell6_lane", 1:length(fastq), sep="_") 
(fq <- readFastq(fastq[1])) # Imports first FASTQ file
{% endhighlight %}

{% highlight txt %}
## class: ShortReadQ
## length: 1000 reads; width: 36 cycles
{% endhighlight %}

{% highlight r %}
countLines(dirPath="./data", pattern=".fastq$")/4 # Counts numbers of reads in FASTQ files
{% endhighlight %}

{% highlight txt %}
## SRR038845.fastq SRR038846.fastq SRR038848.fastq SRR038850.fastq 
##            1000            1000            1000            1000
{% endhighlight %}

{% highlight r %}
id(fq)[1] # Returns ID field
{% endhighlight %}

{% highlight txt %}
##   A BStringSet instance of length 1
##     width seq
## [1]    43 SRR038845.3 HWI-EAS038:6:1:0:1938 length=36
{% endhighlight %}

{% highlight r %}
sread(fq)[1] # Returns sequence
{% endhighlight %}

{% highlight txt %}
##   A DNAStringSet instance of length 1
##     width seq
## [1]    36 CAACGAGTTCACACCTTGGCCGACAGGCCCGGGTAA
{% endhighlight %}

{% highlight r %}
quality(fq)[1] # Returns Phred scores 
{% endhighlight %}

{% highlight txt %}
## class: FastqQuality
## quality:
##   A BStringSet instance of length 1
##     width seq
## [1]    36 BA@7>B=>:>>7@7@>>9=BAA?;>52;>:9=8.=A
{% endhighlight %}

{% highlight r %}
as(quality(fq), "matrix")[1:4,1:12] # Coerces Phred scores to numeric matrix
{% endhighlight %}

{% highlight txt %}
##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
## [1,]   33   32   31   22   29   33   28   29   25    29    29    22
## [2,]   33   34   34   33   32   31   33   33   31    33    33    33
## [3,]   33   33   34   33   33   33   33   33   33    31    31    33
## [4,]   33   33   33   33   31   33   28   31   28    32    33    33
{% endhighlight %}

{% highlight r %}
ShortReadQ(sread=sread(fq), quality=quality(fq), id=id(fq)) # Constructs a ShortReadQ from components
{% endhighlight %}

{% highlight txt %}
## class: ShortReadQ
## length: 1000 reads; width: 36 cycles
{% endhighlight %}

## FASTQ Quality Reports 

### Using `systemPipeR`

The following `seeFastq` and `seeFastqPlot` functions generate and plot a series of useful quality statistics for a set of FASTQ files.


{% highlight r %}
library(systemPipeR)
fqlist <- seeFastq(fastq=fastq, batchsize=800, klength=8) # For real data set batchsize to at least 10^5 
seeFastqPlot(fqlist)
{% endhighlight %}

![](Rsequences_files/see_fastq-1.png)

Handles many samples in one PDF file. For more details see [here](http://bioconductor.org/packages/devel/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html) 

### Using `ShortRead`

The `ShortRead` package contains several FASTQ quality reporting functions.

{% highlight r %}
sp <- SolexaPath(system.file('extdata', package='ShortRead'))
fl <- file.path(analysisPath(sp), "s_1_sequence.txt") 
fls <- c(fl, fl) 
coll <- QACollate(QAFastqSource(fls), QAReadQuality(), QAAdapterContamination(), 
	    QANucleotideUse(), QAQualityUse(), QASequenceUse(), QAFrequentSequence(n=10), 
		QANucleotideByCycle(), QAQualityByCycle())
x <- qa2(coll, verbose=TRUE)
res <- report(x)
if(interactive())
browseURL(res) 
{% endhighlight %}

## Filtering and Trimming FASTQ Files with ShortRead 

### Adaptor trimming


{% highlight r %}
fqtrim <- trimLRPatterns(Rpattern="GCCCGGGTAA", subject=fq)
sread(fq)[1:2] # Before trimming
{% endhighlight %}

{% highlight txt %}
##   A DNAStringSet instance of length 2
##     width seq
## [1]    36 CAACGAGTTCACACCTTGGCCGACAGGCCCGGGTAA
## [2]    36 CCAATGATTTTTTTCCGTGTTTCAGAATACGGTTAA
{% endhighlight %}

{% highlight r %}
sread(fqtrim)[1:2] # After trimming
{% endhighlight %}

{% highlight txt %}
##   A DNAStringSet instance of length 2
##     width seq
## [1]    26 CAACGAGTTCACACCTTGGCCGACAG
## [2]    36 CCAATGATTTTTTTCCGTGTTTCAGAATACGGTTAA
{% endhighlight %}
### Read counting and duplicate removal


{% highlight r %}
tables(fq)$distribution # Counts read occurences
{% endhighlight %}

{% highlight txt %}
##   nOccurrences nReads
## 1            1    948
## 2            2     26
{% endhighlight %}

{% highlight r %}
sum(srduplicated(fq)) # Identifies duplicated reads
{% endhighlight %}

{% highlight txt %}
## [1] 26
{% endhighlight %}

{% highlight r %}
fq[!srduplicated(fq)]
{% endhighlight %}

{% highlight txt %}
## class: ShortReadQ
## length: 974 reads; width: 36 cycles
{% endhighlight %}

### Trimming low quality tails


{% highlight r %}
cutoff <- 30
cutoff <- rawToChar(as.raw(cutoff+33))
sread(trimTails(fq, k=2, a=cutoff, successive=FALSE))[1:2]
{% endhighlight %}

{% highlight txt %}
##   A DNAStringSet instance of length 2
##     width seq
## [1]     4 CAAC
## [2]    20 CCAATGATTTTTTTCCGTGT
{% endhighlight %}

### Removal of reads with Phred scores below a threshold value


{% highlight r %}
cutoff <- 30
qcount <- rowSums(as(quality(fq), "matrix") <= 20) 
fq[qcount == 0] # Number of reads where all Phred scores >= 20
{% endhighlight %}

{% highlight txt %}
## class: ShortReadQ
## length: 349 reads; width: 36 cycles
{% endhighlight %}

### Removal of reads with x Ns and/or low complexity segments


{% highlight r %}
filter1 <- nFilter(threshold=1) # Keeps only reads without Ns
filter2 <- polynFilter(threshold=20, nuc=c("A","T","G","C")) # Removes reads with >=20 of one nucleotide
filter <- compose(filter1, filter2)
fq[filter(fq)]
{% endhighlight %}

{% highlight txt %}
## class: ShortReadQ
## length: 989 reads; width: 36 cycles
{% endhighlight %}

## Memory Efficient FASTQ Processing

Streaming through FASTQ files with `FastqStreamer` and random sampling reads with `FastqSampler`


{% highlight r %}
fq <- yield(FastqStreamer(fastq[1], 50)) # Imports first 50 reads 
fq <- yield(FastqSampler(fastq[1], 50)) # Random samples 50 reads 
{% endhighlight %}

{% highlight txt %}
## Warning: closing unused connection 5 (data/SRR038845.fastq)
{% endhighlight %}

Streaming through a FASTQ file while applying filtering/trimming functions and writing the results to a new file
 here `SRR038845.fastq_sub` in `data` directory.


{% highlight r %}
f <- FastqStreamer(fastq[1], 50) 
while(length(fq <- yield(f))) {
	fqsub <- fq[grepl("^TT", sread(fq))] 
	writeFastq(fqsub, paste(fastq[1], "sub", sep="_"), mode="a", compress=FALSE)
}
close(f)
{% endhighlight %}


---
title: Sequences in Bioconductor
keywords: 
last_updated: Wed Apr 20 19:57:09 2016
---

## Important Data Objects in Biostrings

### `XString` for single sequence

* `DNAString`: for DNA
* `RNAString`: for RNA
* `AAString`: for amino acid 
* `BString`: for any string

### `XStringSet` for many sequences
        
* `DNAStringSet``: for DNA
* `RNAStringSet`: for RNA
* `AAStringSet`: for amino acid 
* `BStringSet`: for any string

### `QualityScaleXStringSet` for sequences with quality data

* `QualityScaledDNAStringSet`: for DNA
* `QualityScaledRNAStringSet`: for RNA
* `QualityScaledAAStringSet`: for amino acid 
* `QualityScaledBStringSet`: for any string

## Sequence Import and Export}

Download the following sequences to your current working directory and then import them into R: 
[ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.ffn](ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.ffn)


{% highlight r %}
dir.create("data")
{% endhighlight %}

{% highlight txt %}
## Warning in dir.create("data"): 'data' already exists
{% endhighlight %}

{% highlight r %}
# system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.ffn")
download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.ffn", "data/AE004437.ffn")
myseq <- readDNAStringSet("data/AE004437.ffn")
myseq[1:3]
{% endhighlight %}

{% highlight txt %}
##   A DNAStringSet instance of length 3
##     width seq                                                                   names               
## [1]  1206 ATGACTCGGCGGTCTCGTGTCGGTGCCGGCCTC...GTCGTCGTTGTTCGACGCTGGCGGAACCCATGA gi|12057215|gb|AE...
## [2]   666 ATGAGCATCATCGAACTCGAAGGCGTGGTCAAA...GTCAACCTCGTCGATGGGGTGTTACACACGTGA gi|12057215|gb|AE...
## [3]  1110 ATGGCGTGGCGGAACCTCGGGCGGAACCGCGTG...AACGATCCGCCCGTCGAGGCGCTCGGCGAATGA gi|12057215|gb|AE...
{% endhighlight %}

{% highlight r %}
sub <- myseq[grep("99.*", names(myseq))]
length(sub)
{% endhighlight %}

{% highlight txt %}
## [1] 170
{% endhighlight %}

{% highlight r %}
writeXStringSet(sub, file="./data/AE004437sub.ffn", width=80)
{% endhighlight %}

Now inspect exported sequence file `AE004437sub.ffn` in a text editor
	
    
## Working with `XString` Containers

The `XString` stores the different types of biosequences in dedicated containers

{% highlight r %}
library(Biostrings)
d <- DNAString("GCATAT-TAC")
d
{% endhighlight %}

{% highlight txt %}
##   10-letter "DNAString" instance
## seq: GCATAT-TAC
{% endhighlight %}

{% highlight r %}
d[1:4]
{% endhighlight %}

{% highlight txt %}
##   4-letter "DNAString" instance
## seq: GCAT
{% endhighlight %}

{% highlight r %}
r <- RNAString("GCAUAU-UAC") 
r <- RNAString(d) # Converts d to RNAString object
p <- AAString("HCWYHH")
b <- BString("I store any set of characters. Other XString objects store only the IUPAC characters.")
{% endhighlight %}

## Working with `XStringSet` Containers

`XStringSet` containers allow to store many biosequences in one object:

{% highlight r %}
dset <- DNAStringSet(c("GCATATTAC", "AATCGATCC", "GCATATTAC")) 
names(dset) <- c("seq1", "seq2", "seq3") # Assigns names
dset[1:2]
{% endhighlight %}

{% highlight txt %}
##   A DNAStringSet instance of length 2
##     width seq                                                                   names               
## [1]     9 GCATATTAC                                                             seq1
## [2]     9 AATCGATCC                                                             seq2
{% endhighlight %}

{% highlight r %}
width(dset) # Returns the length of each sequences
{% endhighlight %}

{% highlight txt %}
## [1] 9 9 9
{% endhighlight %}

{% highlight r %}
d <- dset[[1]] # The [[ subsetting operator returns a single entry as XString object
dset2 <- c(dset, dset) # Appends/concatenates two XStringSet objects
dsetchar <- as.character(dset) # Converts XStringSet to named vector 
dsetone <- unlist(dset) # Collapses many sequences to a single one stored in a DNAString container
{% endhighlight %}
Sequence subsetting by positions:

{% highlight r %}
DNAStringSet(dset, start=c(1,2,3), end=c(4,8,5)) 
{% endhighlight %}

{% highlight txt %}
##   A DNAStringSet instance of length 3
##     width seq                                                                   names               
## [1]     4 GCAT                                                                  seq1
## [2]     7 ATCGATC                                                               seq2
## [3]     3 ATA                                                                   seq3
{% endhighlight %}

## Multiple Alignment Class

The `XMultipleAlignment` class stores the different types of multiple sequence alignments:


{% highlight r %}
origMAlign <- readDNAMultipleAlignment(filepath = system.file("extdata",
              "msx2_mRNA.aln", package = "Biostrings"), format = "clustal")
origMAlign
{% endhighlight %}

{% highlight txt %}
## DNAMultipleAlignment with 8 rows and 2343 columns
##      aln                                                                        names               
## [1] -----TCCCGTCTCCGCAGCAAAAAAGTTTGAGTCG...TTGTCCAAACTCACAATTAAAAAAAAAAAAAAAAA gi|84452153|ref|N...
## [2] ------------------------------------...----------------------------------- gi|208431713|ref|...
## [3] ------------------------------------...----------------------------------- gi|118601823|ref|...
## [4] ----------------------AAAAGTTGGAGTCT...----------------------------------- gi|114326503|ref|...
## [5] ------------------------------------...----------------------------------- gi|119220589|ref|...
## [6] ------------------------------------...----------------------------------- gi|148540149|ref|...
## [7] --------------CGGCTCCGCAGCGCCTCACTCG...----------------------------------- gi|45383056|ref|N...
## [8] GGGGGAGACTTCAGAAGTTGTTGTCCTCTCCGCTGA...----------------------------------- gi|213515133|ref|...
{% endhighlight %}

## Basic Sequence Manipulations

### Reversed & Complement


{% highlight r %}
randset <- DNAStringSet(rand)
complement(randset[1:2])
{% endhighlight %}

{% highlight txt %}
##   A DNAStringSet instance of length 2
##     width seq
## [1]    10 TCTAGATCGC
## [2]    16 ATTGTAGCCTAAGATT
{% endhighlight %}

{% highlight r %}
reverse(randset[1:2])
{% endhighlight %}

{% highlight txt %}
##   A DNAStringSet instance of length 2
##     width seq
## [1]    10 GCGATCTAGA
## [2]    16 AATCTTAGGCTACAAT
{% endhighlight %}

{% highlight r %}
reverseComplement(randset[1:2])
{% endhighlight %}

{% highlight txt %}
##   A DNAStringSet instance of length 2
##     width seq
## [1]    10 CGCTAGATCT
## [2]    16 TTAGAATCCGATGTTA
{% endhighlight %}

## Translate DNA into Protein

{% highlight r %}
translate(randset[1:2])
{% endhighlight %}

{% highlight txt %}
## Warning in .Call2("DNAStringSet_translate", x, skip_code, dna_codes[codon_alphabet], : in 'x[[1]]':
## last base was ignored
{% endhighlight %}

{% highlight txt %}
## Warning in .Call2("DNAStringSet_translate", x, skip_code, dna_codes[codon_alphabet], : in 'x[[2]]':
## last base was ignored
{% endhighlight %}

{% highlight txt %}
##   A AAStringSet instance of length 2
##     width seq
## [1]     3 RSS
## [2]     5 *HRIL
{% endhighlight %}

## Pattern Matching

### Pattern matching with mismatches

Find pattern matches in reference 

{% highlight r %}
myseq1 <- readDNAStringSet("./data/AE004437.ffn") 
mypos <- matchPattern("ATGGTG", myseq1[[1]], max.mismatch=1) 
{% endhighlight %}

Count only the corresponding matches

{% highlight r %}
countPattern("ATGGCT", myseq1[[1]], max.mismatch=1) 
{% endhighlight %}

{% highlight txt %}
## [1] 3
{% endhighlight %}

Count only the matches in many sequences

{% highlight r %}
vcountPattern("ATGGCT", myseq1, max.mismatch=1)[1:20]
{% endhighlight %}

{% highlight txt %}
##  [1] 3 0 5 4 1 2 2 1 4 3 0 0 1 2 0 1 4 0 0 1
{% endhighlight %}

Results shown in DNAStringSet object

{% highlight r %}
tmp <- c(DNAStringSet("ATGGTG"), DNAStringSet(mypos)) 
{% endhighlight %}

Return a consensus  matrix for query and hits.

{% highlight r %}
consensusMatrix(tmp)[1:4,] 
{% endhighlight %}

{% highlight txt %}
##   [,1] [,2] [,3] [,4] [,5] [,6]
## A    3    0    0    0    0    0
## C    1    1    0    0    0    0
## G    0    0    4    4    1    4
## T    0    3    0    0    3    0
{% endhighlight %}

Find all pattern matches in reference

{% highlight r %}
myvpos <- vmatchPattern("ATGGCT", myseq1, max.mismatch=1) 
myvpos # The results are stored as MIndex object.
{% endhighlight %}

{% highlight txt %}
## MIndex object of length 2058
## $`gi|12057215|gb|AE004437.1|:248-1453 Halobacterium sp. NRC-1, complete genome`
## IRanges of length 3
##     start end width
## [1]     1   6     6
## [2]   383 388     6
## [3]   928 933     6
## 
## $`gi|12057215|gb|AE004437.1|:1450-2115 Halobacterium sp. NRC-1, complete genome`
## IRanges of length 0
## 
## $`gi|12057215|gb|AE004437.1|:2145-3254 Halobacterium sp. NRC-1, complete genome`
## IRanges of length 5
##     start end width
## [1]     1   6     6
## [2]    94  99     6
## [3]   221 226     6
## [4]   535 540     6
## [5]   601 606     6
## 
## ...
## <2055 more elements>
{% endhighlight %}

{% highlight r %}
Views(myseq1[[1]], start(myvpos[[1]]), end(myvpos[[1]])) # Retrieves the result for single entry
{% endhighlight %}

{% highlight txt %}
##   Views on a 1206-letter DNAString subject
## subject: ATGACTCGGCGGTCTCGTGTCGGTGCCGGCCTCGCAGCCATTGT...TTGCGATCGTCGTCGTCGTTGTTCGACGCTGGCGGAACCCATGA
## views:
##     start end width
## [1]     1   6     6 [ATGACT]
## [2]   383 388     6 [ATGGCA]
## [3]   928 933     6 [ATGACT]
{% endhighlight %}

Return all matches

{% highlight r %}
sapply(seq(along=myseq1), function(x) 
       as.character(Views(myseq1[[x]], start(myvpos[[x]]), end(myvpos[[x]])))) 
{% endhighlight %}

### Pattern matching with regular expression support


{% highlight r %}
myseq <- DNAStringSet(c("ATGCAGACATAGTG", "ATGAACATAGATCC", "GTACAGATCAC"))
myseq[grep("^ATG", myseq, perl=TRUE)] # String searching with regular expression support
{% endhighlight %}

{% highlight txt %}
##   A DNAStringSet instance of length 2
##     width seq
## [1]    14 ATGCAGACATAGTG
## [2]    14 ATGAACATAGATCC
{% endhighlight %}

{% highlight r %}
pos1 <- regexpr("AT", myseq) # Searches 'myseq' for first match of pattern "AT"
as.numeric(pos1); attributes(pos1)$match.length # Returns position information of matches
{% endhighlight %}

{% highlight txt %}
## [1] 1 1 7
{% endhighlight %}

{% highlight txt %}
## [1] 2 2 2
{% endhighlight %}

{% highlight r %}
pos2 <- gregexpr("AT", myseq) # Searches 'myseq' for all matches of pattern "AT"
as.numeric(pos2[[1]]); attributes(pos2[[1]])$match.length # Match positions in first sequence
{% endhighlight %}

{% highlight txt %}
## [1] 1 9
{% endhighlight %}

{% highlight txt %}
## [1] 2 2
{% endhighlight %}

{% highlight r %}
DNAStringSet(gsub("^ATG", "NNN", myseq)) # String substitution with regular expression support
{% endhighlight %}

{% highlight txt %}
##   A DNAStringSet instance of length 3
##     width seq
## [1]    14 NNNCAGACATAGTG
## [2]    14 NNNAACATAGATCC
## [3]    11 GTACAGATCAC
{% endhighlight %}

## PWM Viewing and Searching


{% highlight r %}
pwm <- PWM(DNAStringSet(c("GCT", "GGT", "GCA"))) 
library(seqLogo); seqLogo(t(t(pwm) * 1/colSums(pwm)))
{% endhighlight %}

{% highlight txt %}
## Loading required package: grid
{% endhighlight %}

![](Rsequences_files/pwm_logo-1.png)

{% highlight r %}
chr <- DNAString("AAAGCTAAAGGTAAAGCAAAA") 
matchPWM(pwm, chr, min.score=0.9) # Searches sequence for PWM matches with score better than min.score.
{% endhighlight %}

{% highlight txt %}
##   Views on a 21-letter DNAString subject
## subject: AAAGCTAAAGGTAAAGCAAAA
## views:
##     start end width
## [1]     4   6     3 [GCT]
## [2]    10  12     3 [GGT]
## [3]    16  18     3 [GCA]
{% endhighlight %}


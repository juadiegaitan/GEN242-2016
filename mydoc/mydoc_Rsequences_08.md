---
title: Homework 6
keywords: 
last_updated: Wed Apr 20 19:57:09 2016
---

## HW 6a - Demultiplexing
	
Write a demultiplexing function that accepts any number of
barcodes and splits a FASTQ file into as many subfiles as there are barcodes.
At the same time the function should remove low quality tails from the reads.
The following function accomplishes the first step. Expand this function so
that it performs the second step as well.  \end{itemize} 


{% highlight r %}
demultiplex <- function(x, barcode, nreads) {
	f <- FastqStreamer(x, nreads) 
	while(length(fq <- yield(f))) {
		for(i in barcode) {
			pattern <- paste("^", i, sep="")
			fqsub <- fq[grepl(pattern, sread(fq))] 
			if(length(fqsub) > 0) {
				writeFastq(fqsub, paste(x, i, sep="_"), mode="a", compress=FALSE)
			}
		}
	}
	close(f)
}
demultiplex(x=fastq[1], barcode=c("TT", "AA", "GG"), nreads=50)
{% endhighlight %}

## HW6b - Sequence Parsing 

* Download `GFF` from _Halobacterium sp_  [here](ftp://ftp.ncbi.nih.gov/genbank/genomes/Bacteria/Halobacterium\_sp\_uid217/AE004437.gff)
* Download genome sequence from -Halobacterium sp_ [here](ftp://ftp.ncbi.nih.gov/genbank/genomes/Bacteria/Halobacterium\_sp\_uid217/AE004437.fna)
* __Task 1__ Extract gene ranges, parse their sequences from genome and translate them into proteins
* __Task 2__ Reduce overlapping genes and parse their sequences from genome
* __Task 3__ Generate intergenic ranges and parse their sequences from genome

__Useful commands__

{% highlight r %}
download.file("ftp://ftp.ncbi.nih.gov/genbank/genomes/Bacteria/Halobacterium_sp_uid217/AE004437.gff", "data/AE004437.gff")
download.file("ftp://ftp.ncbi.nih.gov/genbank/genomes/Bacteria/Halobacterium_sp_uid217/AE004437.fna", "data/AE004437.fna")
chr <- readDNAStringSet("data/AE004437.fna")
gff <- import("data/AE004437.gff", asRangedData=FALSE)
gffgene <- gff[values(gff)[,"type"]=="gene"]
gene <- DNAStringSet(Views(chr[[1]], IRanges(start(gffgene), end(gffgene))))
names(gene) <- values(gffgene)[,"locus_tag"]
pos <- values(gffgene[strand(gffgene) == "+"])[,"locus_tag"]
p1 <- translate(gene[names(gene) %in% pos])
names(p1) <- names(gene[names(gene) %in% pos])
neg <- values(gffgene[strand(gffgene) == "-"])[,"locus_tag"]
p2 <- translate(reverseComplement(gene[names(gene) %in% neg]))
names(p2) <- names(gene[names(gene) %in% neg])
writeXStringSet(c(p1, p2), "mypep.fasta")
{% endhighlight %}

\end{changemargin}
\end{frame}


## Homework submission
Submit the homework results in one well structured and annotated R script to the instructor. The script should include instructions on how to use the functions.

## Due date

This homework is due on Thu, April 28th at 6:00 PM.



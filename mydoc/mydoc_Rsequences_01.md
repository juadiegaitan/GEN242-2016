---
title: Sequence Analysis 
keywords: 
last_updated: Wed Apr 20 19:57:09 2016
---
Author: Thomas Girke

Last update: 20 April, 2016 

Alternative formats of this vignette:
[`Single-page .Rmd HTML`](https://htmlpreview.github.io/?https://github.com/tgirke/GEN242/blob/master/vignettes/09_Programming_in_R/Programming_in_R.html),
[`.Rmd`](https://raw.githubusercontent.com/tgirke/GEN242/master/vignettes/09_Programming_in_R/Programming_in_R.Rmd),
[`.R`](https://raw.githubusercontent.com/tgirke/GEN242/master/vignettes/09_Programming_in_R/Programming_in_R.R)
[Old Slides `.pdf`](https://drive.google.com/file/d/0B-lLYVUOliJFSnplWlZWcllYY0E/view?usp=sharing)

## Overview

## Sequence Analysis in R and Bioconductor

__R Base__

* Some basic string handling utilities. Wide spectrum of numeric data analysis tools.

__Bioconductor__

Bioconductor packages provide much more sophisticated string handling utilities for sequence analysis.

* [Biostrings](http://bioconductor.org/packages/release/bioc/html/Biostrings.html): general sequence analysis environment
* [ShortRead](http://bioconductor.org/packages/release/bioc/html/ShortRead.html): pipeline for short read data
* [IRanges](http://bioconductor.org/packages/release/bioc/html/IRanges.html): low-level infrastructure for range data
* [GenomicRanges](http://bioconductor.org/packages/release/bioc/html/GenomicRanges.html): high-level infrastructure for range data
* [GenomicFeatures](http://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html): managing transcript centric annotations
* [GenomicAlignments](http://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html): handling short genomic alignments
* [Rsamtools](http://bioconductor.org/packages/release/bioc/html/Rsamtools.html): interface to  `samtools`, `bcftools` and `tabix` 
* [BSgenome](http://bioconductor.org/packages/release/bioc/html/BSgenome.html): genome annotation data
* [biomaRt](http://bioconductor.org/packages/release/bioc/html/biomaRt.html): interface to BioMart annotations
* [rtracklayer](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html): Annotation imports, interface to online genome browsers


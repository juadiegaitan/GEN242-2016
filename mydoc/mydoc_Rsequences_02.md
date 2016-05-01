---
title: Package Requirements
keywords: 
last_updated: Sat Apr 30 19:54:48 2016
---

Several Bioconductor packages are required for this tutorial. To install them, execute
the following lines in the R console. Please also make sure that you have a recent R version
installed on your system. R versions `3.2.x` or higher are recommended.


{% highlight r %}
source("https://bioconductor.org/biocLite.R")
biocLite(c("Biostrings", "GenomicRanges", "GenomicRanges", "rtracklayer", "systemPipeR", "seqLogo", "ShortRead"))
{% endhighlight %}


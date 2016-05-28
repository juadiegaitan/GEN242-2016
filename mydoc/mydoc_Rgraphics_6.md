---
title: Specialty Graphics
keywords: 
last_updated: Fri May 27 17:09:46 2016
---

## Venn Diagrams 


{% highlight r %}
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R")
setlist5 <- list(A=sample(letters, 18), B=sample(letters, 16), C=sample(letters, 20), D=sample(letters, 22), E=sample(letters, 18))
OLlist5 <- overLapper(setlist=setlist5, sep="_", type="vennsets")
counts <- sapply(OLlist5$Venn_List, length)
vennPlot(counts=counts, ccol=c(rep(1,30),2), lcex=1.5, ccex=c(rep(1.5,5), rep(0.6,25),1.5))
{% endhighlight %}

![](Rgraphics_files/specgraph_venn-1.png)

## Compound Structures 

Plots depictions of small molecules with `ChemmineR` package


{% highlight r %}
library(ChemmineR)
{% endhighlight %}

{% highlight txt %}
## 
## Attaching package: 'ChemmineR'
{% endhighlight %}

{% highlight txt %}
## The following object is masked from 'package:ShortRead':
## 
##     view
{% endhighlight %}

{% highlight txt %}
## The following object is masked from 'package:S4Vectors':
## 
##     fold
{% endhighlight %}

{% highlight r %}
data(sdfsample)
plot(sdfsample[1], print=FALSE)
{% endhighlight %}

![](Rgraphics_files/specgraph_structure-1.png)

## ROC Plots

A variety of libraries are available for plotting receiver operating characteristic (ROC) curves in R:

+ [ROCR](http://rocr.bioinf.mpi-sb.mpg.de/)
+ [ROC](http://bioconductor.org/packages/release/bioc/html/ROC.html)
+ [pROC](http://web.expasy.org/pROC/)
+ [ggplot2](http://largedata.blogspot.com/2011/07/plotting-roc-curves-in-ggplot2.html)

## Trees 

The [`ape`](http://ape-package.ird.fr/ape_screenshots.html) package provides many useful utilities for phylogenetic analysis and tree plotting. Another useful package for 
plotting trees is [`ggtree`](http://bioconductor.org/packages/release/bioc/html/ggtree.html). The following example plots two trees face to face with links to identical
leaf labels.


{% highlight r %}
library(ape)
tree1 <- rtree(40)
tree2 <- rtree(20)
association <- cbind(tree2$tip.label, tree2$tip.label)
cophyloplot(tree1, tree2, assoc = association,
            length.line = 4, space = 28, gap = 3)
{% endhighlight %}

![](Rgraphics_files/trees_ape1-1.png)




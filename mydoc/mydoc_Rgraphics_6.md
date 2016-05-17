---
title: Specialty Graphics
keywords: 
last_updated: Tue May 17 09:39:48 2016
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


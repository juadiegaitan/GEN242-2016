---
title: Reading and Writing External Data
keywords: 
last_updated: Thu Apr  7 10:17:06 2016
---
## Import data from tabular files into R


{% highlight r %}
myDF <- read.delim("myData.xls", sep="\t")
{% endhighlight %}

## Export data from R to tabular files

{% highlight r %}
write.table(myDF, file="myfile.xls", sep="\t", quote=FALSE, col.names=NA)
{% endhighlight %}

## Copy and paste into R

On Windows/Linux systems

{% highlight r %}
read.delim("clipboard") 
{% endhighlight %}
On Mac OS X systems

{% highlight r %}
read.delim(pipe("pbpaste")) 
{% endhighlight %}

## Copy and paste from R 

On Windows/Linux systems

{% highlight r %}
write.table(iris, "clipboard", sep="\t", col.names=NA, quote=F) 
{% endhighlight %}

On Mac OS X systems

{% highlight r %}
zz <- pipe('pbcopy', 'w')
write.table(iris, zz, sep="\t", col.names=NA, quote=F)
close(zz) 
{% endhighlight %}

## Homework 3A 

Homework 3A: [Object Subsetting Routines and Import/Export](http://girke.bioinformatics.ucr.edu/GEN242/mydoc/mydoc_homework_03.html)



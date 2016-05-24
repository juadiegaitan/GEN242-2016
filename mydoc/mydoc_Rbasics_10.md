---
title: Reading and Writing External Data
keywords: 
last_updated: Mon May 23 21:13:38 2016
---
## Import of tabular data

Import of a tab-delimited tabular file

{% highlight r %}
myDF <- read.delim("myData.xls", sep="\t")
{% endhighlight %}

Import of Excel file. Note: working with tab- or comma-delimited files is more flexible and preferred.

{% highlight r %}
library(gdata)
myDF <- read.xls"myData.xls")
{% endhighlight %}

## Export of tabular data

{% highlight r %}
write.table(myDF, file="myfile.xls", sep="\t", quote=FALSE, col.names=NA)
{% endhighlight %}

## Line-wise import

{% highlight r %}
myDF <- readLines("myData.txt")
{% endhighlight %}

## Line-wise export

{% highlight r %}
writeLines(month.name, "myData.txt")
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


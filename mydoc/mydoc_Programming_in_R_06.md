---
title: Running R Scripts
keywords: 
last_updated: Wed Apr 13 22:20:52 2016
---

## Possibilities for Executing R Scripts

### R console

{% highlight r %}
source("my_script.R")
{% endhighlight %}

### Command-line


{% highlight sh %}
Rscript my_script.R # or just ./myscript.R after making it executable
R CMD BATCH my_script.R # Alternative way 1 
R --slave < my_script.R # Alternative way 2
{% endhighlight %}
### Passing arguments from command-line to R

Create an R script named `test.R` with the following content:


{% highlight sh %}
myarg <- commandArgs()
print(iris[1:myarg[6], ])
{% endhighlight %}

Then run it from the command-line like this:

{% highlight sh %}
Rscript test.R 10
{% endhighlight %}

In the given example the number `10` is passed on from the command-line as an argument to the R script which is used to return to `STDOUT` the first 10 rows of the `iris` sample data. If several arguments are provided, they will be interpreted as one string and need to be split in R with the strsplit function. A more detailed example can be found [here](http://manuals.bioinformatics.ucr.edu/home/ht-seq\#TOC-Quality-Reports-of-FASTQ-Files-).


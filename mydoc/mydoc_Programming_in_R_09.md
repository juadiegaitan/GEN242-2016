---
title: Homework 5
keywords: 
last_updated: Wed Apr 13 17:35:30 2016
---

## Reverse and complement of DNA

__Task 1__: Write a `RevComp` function that returns the reverse and complement of a DNA sequence string. Include an argument that will allow to return only the reversed sequence, the complemented sequence or the reversed and complemented sequence. The following R functions will be useful for the implementation: 


{% highlight r %}
x <- c("ATGCATTGGACGTTAG")  
x <- substring(x, 1:nchar(x), 1:nchar(x)) 
x <- rev(x) 
x <- paste(x, collapse="")
chartr("ATGC", "TACG", x) 
{% endhighlight %}

__Task 2__: Write a function that applies the `RevComp` function to many sequences stored in a vector.

## Translate DNA into Protein

__Task 3__: Write a function that will translate one or many DNA sequences in all three reading frames into proteins. The following commands will simplify this task:


{% highlight r %}
AAdf <- read.table(file="http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/AA.txt", header=TRUE, sep="\t") 
AAv <- as.character(AAdf[,2]) 
names(AAv) <- AAdf[,1] 
y <- gsub("(...)", "\\1_", x) 
y <- unlist(strsplit(y, "_")) 
y <- y[grep("^...\$", y)] 
AAv[y] 
{% endhighlight %}
Submit the 3 functions in one well structured and annotated R script to the instructor. The script should include instructions on how to use the functions.



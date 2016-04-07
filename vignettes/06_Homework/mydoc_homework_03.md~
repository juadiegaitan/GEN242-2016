---
title: HW3 - Introduction to R
last_updated: 07-April-16
---

## A. Object Subsetting, Import and Export

- __Task 1__: Sort the rows of the `iris` data frame by its first column and sort its columns alphabetically by column names.
- __Task 2__: Subset the first 12 rows, export the result to a text file and view it in a spreadsheet program like Excel or Google Sheets. 
- __Task 3__: Change some column titles in your spreadsheet program and import the result into R.  

Before you start it can be helpful to evaluate the structure of the `iris` data set with the following commands:
{% highlight r %}
class(iris)
dim(iris)
colnames(iris)
{% endhighlight %}

<!---
Solution
{% highlight r %}
irismod <- iris[order(iris[,1]), order(colnames(iris))]
irismod <- irismod[1:12,]
write.table(irismod, file="irismod.xls", sep="\t", quote=FALSE, row.names=FALSE)
irisimport <- read.delim(file="irismod.xls", sep="\t")
{% endhighlight %}
-->

## Homework submission

Upload R script solving the homework assignemts to your private course GitHub repository under `Homework/HW3/HW3.R`.

## Due date

Most homeworks will be due one week after they are assigned. This one is due on Thu, April 14th at 6:00 PM.

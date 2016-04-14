---
title: Useful Utilities
keywords: 
last_updated: Wed Apr 13 16:38:00 2016
---

## Debugging Utilities

Several debugging utilities are available for R. They include:

* `traceback`
* `browser`
* `options(error=recover)`
* `options(error=NULL)`
* `debug`

The [Debugging in R page](http://www.stats.uwo.ca/faculty/murdoch/software/debuggingR/) provides an overview of the available resources.

## Regular Expressions

R's regular expression utilities work similar as in other languages. To learn how to use them in R, one can consult the main help page on this topic with `?regexp`.

The grep function can be used for finding patterns in strings, here letter `A` in vector `month.name`.

{% highlight r %}
month.name[grep("A", month.name)] 
{% endhighlight %}

{% highlight txt %}
## [1] "April"  "August"
{% endhighlight %}
Example for using regular expressions to substitute a pattern by another one using a back reference. Remember: single escapes `\` need to be double escaped `\\` in R.


{% highlight r %}
gsub('(i.*a)', 'xxx_\\1', "virginica", perl = TRUE) 
{% endhighlight %}

{% highlight txt %}
## [1] "vxxx_irginica"
{% endhighlight %}

## Interpreting a Character String as Expression

Some useful examples

Generates vector of object names in session

{% highlight r %}
mylist <- ls()
mylist[1] 
{% endhighlight %}

{% highlight txt %}
## [1] "i"
{% endhighlight %}

Executes 1st entry as expression


{% highlight r %}
get(mylist[1])
{% endhighlight %}

{% highlight txt %}
## [1] 150
{% endhighlight %}

Alternative approach 

{% highlight r %}
eval(parse(text=mylist[1])) 
{% endhighlight %}

{% highlight txt %}
## [1] 150
{% endhighlight %}

## Replacement, Split and Paste Functions for Strings

__Selected examples__

Substitution with back reference which inserts in this example `_` character

{% highlight r %}
x <- gsub("(a)","\\1_", month.name[1], perl=T) 
x
{% endhighlight %}

{% highlight txt %}
## [1] "Ja_nua_ry"
{% endhighlight %}

Split string on inserted character from above

{% highlight r %}
strsplit(x,"_")
{% endhighlight %}

{% highlight txt %}
## [[1]]
## [1] "Ja"  "nua" "ry"
{% endhighlight %}

Reverse a character string by splitting first all characters into vector fields


{% highlight r %}
paste(rev(unlist(strsplit(x, NULL))), collapse="") 
{% endhighlight %}

{% highlight txt %}
## [1] "yr_aun_aJ"
{% endhighlight %}

## Time, Date and Sleep

__Selected examples__

Return CPU (and other) times that an expression used (here ls)

{% highlight r %}
system.time(ls()) 
{% endhighlight %}

{% highlight txt %}
##    user  system elapsed 
##       0       0       0
{% endhighlight %}

Return the current system date and time

{% highlight r %}
date() 
{% endhighlight %}

{% highlight txt %}
## [1] "Wed Apr 13 16:06:06 2016"
{% endhighlight %}

Pause execution of R expressions for a given number of seconds (e.g. in loop)

{% highlight r %}
Sys.sleep(1) 
{% endhighlight %}

### Example

#### Import of Specific File Lines with Regular Expression

The following example demonstrates the retrieval of specific lines from an external file with a regular expression. First, an external file is created with the `cat` function, all lines of this file are imported into a vector with `readLines`, the specific elements (lines) are then retieved with the `grep` function, and the resulting lines are split into vector fields with `strsplit`.


{% highlight r %}
cat(month.name, file="zzz.txt", sep="\n")
x <- readLines("zzz.txt")
x[1:6] 
{% endhighlight %}

{% highlight txt %}
## [1] "January"  "February" "March"    "April"    "May"      "June"
{% endhighlight %}

{% highlight r %}
x <- x[c(grep("^J", as.character(x), perl = TRUE))]
t(as.data.frame(strsplit(x, "u")))
{% endhighlight %}

{% highlight txt %}
##                 [,1]  [,2] 
## c..Jan....ary.. "Jan" "ary"
## c..J....ne..    "J"   "ne" 
## c..J....ly..    "J"   "ly"
{% endhighlight %}
## Calling External Software

How to run External command-line software. Here example for running `blastall` from R

{% highlight r %}
system("blastall -p blastp -i seq.fasta -d uniprot -o seq.blastp")
{% endhighlight %}


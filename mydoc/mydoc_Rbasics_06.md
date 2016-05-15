---
title: Data Types 
keywords: 
last_updated: Sat May 14 21:39:23 2016
---

## Numeric data

Example: `1, 2, 3, ...`


{% highlight r %}
x <- c(1, 2, 3)
x
{% endhighlight %}

{% highlight txt %}
## [1] 1 2 3
{% endhighlight %}

{% highlight r %}
is.numeric(x)
{% endhighlight %}

{% highlight txt %}
## [1] TRUE
{% endhighlight %}

{% highlight r %}
as.character(x)
{% endhighlight %}

{% highlight txt %}
## [1] "1" "2" "3"
{% endhighlight %}

## Character data

Example: `"a", "b", "c", ...`


{% highlight r %}
x <- c("1", "2", "3")
x
{% endhighlight %}

{% highlight txt %}
## [1] "1" "2" "3"
{% endhighlight %}

{% highlight r %}
is.character(x)
{% endhighlight %}

{% highlight txt %}
## [1] TRUE
{% endhighlight %}

{% highlight r %}
as.numeric(x)
{% endhighlight %}

{% highlight txt %}
## [1] 1 2 3
{% endhighlight %}

## Complex data

Example: mix of both


{% highlight r %}
c(1, "b", 3)
{% endhighlight %}

{% highlight txt %}
## [1] "1" "b" "3"
{% endhighlight %}

## Logical data

Example: `TRUE` of `FALSE`


{% highlight r %}
x <- 1:10 < 5
x  
{% endhighlight %}

{% highlight txt %}
##  [1]  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
{% endhighlight %}

{% highlight r %}
!x
{% endhighlight %}

{% highlight txt %}
##  [1] FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
{% endhighlight %}

{% highlight r %}
which(x) # Returns index for the 'TRUE' values in logical vector
{% endhighlight %}

{% highlight txt %}
## [1] 1 2 3 4
{% endhighlight %}


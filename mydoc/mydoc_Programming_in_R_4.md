---
title: Data Types 
keywords: 
last_updated: Wed Apr 13 12:24:22 2016
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



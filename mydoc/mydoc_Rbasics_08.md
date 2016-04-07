---
title: Important Utilities
keywords: 
last_updated: Thu Apr  7 06:48:52 2016
---
	
## Combining Objects

The `c` function combines vectors and lists


{% highlight r %}
c(1, 2, 3)
{% endhighlight %}

{% highlight txt %}
## [1] 1 2 3
{% endhighlight %}

{% highlight r %}
x <- 1:3; y <- 101:103
c(x, y)
{% endhighlight %}

{% highlight txt %}
## [1]   1   2   3 101 102 103
{% endhighlight %}

{% highlight r %}
iris$Species[1:8]
{% endhighlight %}

{% highlight txt %}
## [1] setosa setosa setosa setosa setosa setosa setosa setosa
## Levels: setosa versicolor virginica
{% endhighlight %}

The `cbind` and `rbind` functions can be used to append columns and rows, respecively.

{% highlight r %}
ma <- cbind(x, y)
ma
{% endhighlight %}

{% highlight txt %}
##      x   y
## [1,] 1 101
## [2,] 2 102
## [3,] 3 103
{% endhighlight %}

{% highlight r %}
rbind(ma, ma)
{% endhighlight %}

{% highlight txt %}
##      x   y
## [1,] 1 101
## [2,] 2 102
## [3,] 3 103
## [4,] 1 101
## [5,] 2 102
## [6,] 3 103
{% endhighlight %}

## Accessing Dimensions of Objects

Length and dimension information of objects


{% highlight r %}
length(iris$Species)
{% endhighlight %}

{% highlight txt %}
## [1] 150
{% endhighlight %}

{% highlight r %}
dim(iris)
{% endhighlight %}

{% highlight txt %}
## [1] 150   5
{% endhighlight %}

## Accessing Name Slots of Objects

Accessing row and column names of 2D objects

{% highlight r %}
rownames(iris)[1:8]
{% endhighlight %}

{% highlight txt %}
## [1] "1" "2" "3" "4" "5" "6" "7" "8"
{% endhighlight %}

{% highlight r %}
colnames(iris)
{% endhighlight %}

{% highlight txt %}
## [1] "Sepal.Length" "Sepal.Width"  "Petal.Length" "Petal.Width"  "Species"
{% endhighlight %}

Return name field of vectors and lists

{% highlight r %}
names(myVec)
{% endhighlight %}

{% highlight txt %}
##  [1] "A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M" "N" "O" "P" "Q" "R" "S" "T" "U" "V" "W" "X"
## [25] "Y" "Z"
{% endhighlight %}

{% highlight r %}
names(myL)
{% endhighlight %}

{% highlight txt %}
## [1] "name"        "wife"        "no.children" "child.ages"
{% endhighlight %}





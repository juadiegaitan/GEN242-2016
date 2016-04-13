---
title: Loops
keywords: 
last_updated: Wed Apr 13 13:54:52 2016
---

## `for` loop

`for` loops iterate over elements of a looping vector.

__Syntax__

{% highlight r %}
for(variable in sequence) { 
	statements 
}
{% endhighlight %}
__Example__

{% highlight r %}
mydf <- iris
myve <- NULL
for(i in seq(along=mydf[,1])) {
	myve <- c(myve, mean(as.numeric(mydf[i,1:3])))
}
myve[1:8]
{% endhighlight %}

{% highlight txt %}
## [1] 3.333333 3.100000 3.066667 3.066667 3.333333 3.666667 3.133333 3.300000
{% endhighlight %}

__Note:__ Inject into objecs is much faster than append approach with `c`, `cbind`, etc.

__Example__

{% highlight r %}
myve <- numeric(length(mydf[,1]))
for(i in seq(along=myve)) {
	myve[i] <- mean(as.numeric(mydf[i,1:3]))
}
myve[1:8]
{% endhighlight %}

{% highlight txt %}
## [1] 3.333333 3.100000 3.066667 3.066667 3.333333 3.666667 3.133333 3.300000
{% endhighlight %}

### Conditional Stop of Loops

The `stop` function can be used to break out of a loop (or a function) when a condition becomes `TRUE`. In addition, an error message will be printed.

__Example__

{% highlight r %}
x <- 1:10
z <- NULL
for(i in seq(along=x)) { 
	if(x[i] < 5) { 
		z <- c(z, x[i]-1)  
	} else { 
		stop("values need to be < 5") 
	}
}
{% endhighlight %}

## `while` loop

Iterates as long as a condition is true.

__Syntax__

{% highlight r %}
while(condition) {
	statements
}
{% endhighlight %}

__Example__

{% highlight r %}
z <- 0
while(z<5) { 
	z <- z + 2
	print(z)  
}
{% endhighlight %}

{% highlight txt %}
## [1] 2
## [1] 4
## [1] 6
{% endhighlight %}

## The `apply` Function Family

### `apply`

__Syntax__

{% highlight r %}
apply(X, MARGIN, FUN, ARGs)
{% endhighlight %}

__Arguments__

* `X`: `array`, `matrix` or `data.frame`
* `MARGIN`: `1` for rows, `2` for columns
* `FUN`: one or more functions
* `ARGs`: possible arguments for functions

__Example__

{% highlight r %}
apply(iris[1:8,1:3], 1, mean)
{% endhighlight %}

{% highlight txt %}
##        1        2        3        4        5        6        7        8 
## 3.333333 3.100000 3.066667 3.066667 3.333333 3.666667 3.133333 3.300000
{% endhighlight %}

### `tapply`

Applies a function to vector components that are defined by a factor.

__Syntax__

{% highlight r %}
tapply(vector, factor, FUN)
{% endhighlight %}

__Example__

{% highlight r %}
iris[1:2,]
{% endhighlight %}

{% highlight txt %}
##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
## 1          5.1         3.5          1.4         0.2  setosa
## 2          4.9         3.0          1.4         0.2  setosa
{% endhighlight %}

{% highlight r %}
tapply(iris$Sepal.Length, iris$Species, mean)
{% endhighlight %}

{% highlight txt %}
##     setosa versicolor  virginica 
##      5.006      5.936      6.588
{% endhighlight %}

### `sapply` and `lapply`

Both apply a function to vector or list objects. The `lapply` function always returns a list object, while `sapply` returns `vector` or `matrix` objects when it is possible. 

__Examples__

{% highlight r %}
x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE,FALSE,FALSE,TRUE))
lapply(x, mean)
{% endhighlight %}

{% highlight txt %}
## $a
## [1] 5.5
## 
## $beta
## [1] 4.535125
## 
## $logic
## [1] 0.5
{% endhighlight %}

{% highlight r %}
sapply(x, mean)
{% endhighlight %}

{% highlight txt %}
##        a     beta    logic 
## 5.500000 4.535125 0.500000
{% endhighlight %}

Often used in combination with a function definition:

{% highlight r %}
lapply(names(x), function(x) mean(x))
{% endhighlight %}

{% highlight txt %}
## Warning in mean.default(x): argument is not numeric or logical: returning NA

## Warning in mean.default(x): argument is not numeric or logical: returning NA

## Warning in mean.default(x): argument is not numeric or logical: returning NA
{% endhighlight %}

{% highlight txt %}
## [[1]]
## [1] NA
## 
## [[2]]
## [1] NA
## 
## [[3]]
## [1] NA
{% endhighlight %}

{% highlight r %}
sapply(names(x), function(x) mean(x))
{% endhighlight %}

{% highlight txt %}
## Warning in mean.default(x): argument is not numeric or logical: returning NA

## Warning in mean.default(x): argument is not numeric or logical: returning NA

## Warning in mean.default(x): argument is not numeric or logical: returning NA
{% endhighlight %}

{% highlight txt %}
##     a  beta logic 
##    NA    NA    NA
{% endhighlight %}


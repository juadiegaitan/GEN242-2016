---
title: Control Structures
keywords: 
last_updated: Wed Apr 13 22:20:52 2016
---

## Important Operators

### Comparison operators

* `==` (equal)
* `!=` (not equal)
* `>` (greater than)
* `>=` (greater than or equal)
* `<` (less than)
* `<=` (less than or equal)

### Logical operators
		
* `&` (and)
* `|` (or) 
* `!` (not)

## Conditional Executions: `if` Statements

An `if` statement operates on length-one logical vectors.

__Syntax__

{% highlight r %}
if(TRUE) { 
	statements_1 
} else { 
	statements_2 
}
{% endhighlight %}

__Example__

{% highlight r %}
if(1==0) { 
	print(1) 
} else { 
	print(2) 
}
{% endhighlight %}

{% highlight txt %}
## [1] 2
{% endhighlight %}

## Conditional Executions: `ifelse` Statements

The `ifelse` statement operates on vectors.

__Syntax__

{% highlight r %}
ifelse(test, true_value, false_value)
{% endhighlight %}
__Example__

{% highlight r %}
x <- 1:10 
ifelse(x<5, x, 0)
{% endhighlight %}

{% highlight txt %}
##  [1] 1 2 3 4 0 0 0 0 0 0
{% endhighlight %}


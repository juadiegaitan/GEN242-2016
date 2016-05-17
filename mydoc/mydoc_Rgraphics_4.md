---
title: lattice Graphics
keywords: 
last_updated: Tue May 17 09:39:48 2016
---

- What is `lattice`?
    - High-level graphics system 
    - Developed by Deepayan Sarkar 
    - Implements Trellis graphics system from S-Plus
    - Simplifies high-level plotting tasks: arranging complex graphical features 
    - Syntax similar to R's base graphics

- Documentation and Help
    - [Manual](http://lmdvr.r-forge.r-project.org)
    - [Intro](http://www.his.sunderland.ac.uk/~cs0her/Statistics/UsingLatticeGraphicsInR.htm)
    - [Book](http://www.amazon.com/Lattice-Multivariate-Data-Visualization-Use/dp/0387759689)
		
Open a list of all functions available in the lattice package


{% highlight r %}
library(help=lattice) 
{% endhighlight %}

Accessing and changing global parameters:


{% highlight r %}
?lattice.options
?trellis.device
{% endhighlight %}

## Scatter Plot Sample


{% highlight r %}
library(lattice)
p1 <- xyplot(1:8 ~ 1:8 | rep(LETTERS[1:4], each=2), as.table=TRUE) 
plot(p1)
{% endhighlight %}

![](Rgraphics_files/scatter_plot_lattice-1.png)

## Line Plot Sample


{% highlight r %}
library(lattice)
p2 <- parallelplot(~iris[1:4] | Species, iris, horizontal.axis = FALSE, 
              layout = c(1, 3, 1))  
plot(p2)
{% endhighlight %}

![](Rgraphics_files/line_plot_lattice-1.png)


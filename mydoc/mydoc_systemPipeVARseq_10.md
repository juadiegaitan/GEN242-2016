---
title: Venn diagram of variants
keywords: 
last_updated: Thu May  5 11:34:28 2016
---

The venn diagram utilities defined by the `systemPipeR` package can be used to
identify common and unique variants reported for different samples
and/or variant callers. The below generates a 4-way venn diagram
comparing four sampes for each of the two variant callers.


{% highlight r %}
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_gatk_filtered.txt")
varlist <- sapply(names(outpaths(args))[1:4], function(x) as.character(read.delim(outpaths(args)[x])$VARID))
vennset_gatk <- overLapper(varlist, type="vennsets")
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_sambcf_filtered.txt")
varlist <- sapply(names(outpaths(args))[1:4], function(x) as.character(read.delim(outpaths(args)[x])$VARID))
vennset_bcf <- overLapper(varlist, type="vennsets")
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_vartools_filtered.txt")
varlist <- sapply(names(outpaths(args))[1:4], function(x) as.character(read.delim(outpaths(args)[x])$VARID))
vennset_vartools <- overLapper(varlist, type="vennsets")
pdf("./results/vennplot_var.pdf")
vennPlot(list(vennset_gatk, vennset_bcf, vennset_vartools), mymain="", mysub="GATK: red; BCFtools: blue; VariantTools: green", colmode=2, ccol=c("red", "blue", "green"))
dev.off()
{% endhighlight %}

![](systemPipeVARseq_files/vennplot_var.png)
<div align="center">Figure 2: Venn Diagram for 4 samples from GATK and BCFtools</div>



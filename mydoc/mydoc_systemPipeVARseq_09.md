---
title: Summary statistics of variants
keywords: 
last_updated: Thu May  5 11:34:28 2016
---

The `varSummary` function counts the number of variants for each feature type
included in the anntation reports.

## Summary for `GATK`


{% highlight r %}
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_gatk_filtered.txt")
varSummary(args)
write.table(varSummary(args), "./results/variantStats_gatk.xls", quote=FALSE, col.names = NA, sep="\t")
{% endhighlight %}

## Summary for `BCFtools`


{% highlight r %}
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_sambcf_filtered.txt")
varSummary(args)
write.table(varSummary(args), "./results/variantStats_sambcf.xls", quote=FALSE, col.names = NA, sep="\t")
{% endhighlight %}

## Summary for `VariantTools`  


{% highlight r %}
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_vartools_filtered.txt")
varSummary(args)
write.table(varSummary(args), "./results/variantStats_vartools.xls", quote=FALSE, col.names = NA, sep="\t")
{% endhighlight %}


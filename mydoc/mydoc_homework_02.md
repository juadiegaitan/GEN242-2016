---
title: Homework 2 - HW2
last_updated: 07-Mar-16
---

## Topic: Linux Basics

(1) Download Halobacterium proteome and inspect it
{% highlight sh %}
wget http://biocluster.ucr.edu/~tgirke/Linux.sh # downloads code from this slide
wget ftp://ftp.ncbi.nih.gov/genbank/genomes/Bacteria/Halobacterium_sp_uid217/AE004437.faa
less AE004437.faa # press q to quit
{% endhighlight %}

(2) How many protein sequences are stored in the downloaded file?
{% highlight sh %}
grep '>' AE004437.faa | wc
grep '^>' AE004437.faa --count
{% endhighlight %}

(3) How many proteins contain the pattern `WxHxxH` or `WxHxxHH`?
{% highlight sh %}
egrep 'W.H..H{1,2}' AE004437.faa --count
{% endhighlight %}

(4) Use `less` to find IDs for pattern matches or use `awk`
{% highlight sh %}
awk --posix -v RS='>' '/W.H..(H){1,2}/ { print ">" $0;}' AE004437.faa | less
awk --posix -v RS='>' '/W.H..(H){1,2}/ { print ">" $0;}' AE004437.faa | grep '^>' | cut -c 2- | cut -f 1 -d\ > myIDs
{% endhighlight %}

(5) Create a BLASTable database with `formatdb`
{% highlight sh %}
formatdb -i AE004437.faa -p T -o
{% endhighlight %}

(6) Query BLASTable database by IDs stored in a file (_e.g._ `myIDs`)
{% highlight sh %}
fastacmd -d AE004437.faa -i myIDs > myseq.fasta
{% endhighlight %}

(7) Run BLAST search for sequences stored in `myseq.fasta`
{% highlight sh %}
blastall -p blastp -i myseq.fasta -d AE004437.faa -o blastp.out -e 1e-6 -v 10 -b 10
blastall -p blastp -i myseq.fasta -d AE004437.faa -m 8 -e 1e-6 > blastp.tab
{% endhighlight %}

Additional exercise material in [Linux Manual](http://manuals.bioinformatics.ucr.edu/home/linux-basics#TOC-Exercises)

## Homework assignment

Perform above analysis on _Escherichia coli_ HS uid13959: ftp.ncbi.nih.gov/genbank/genomes/Bacteria/Escherichia_coli_HS_uid13959/CP000802.faa. 
Record result from final BLAST command (with `-m 8`) in text file.

## Homework submission

Upload result file to your private course GitHub repository under `Homework/HW2/HW2.txt`.

## Due date

Most homeworks will be due one week after they are assigned. This one is due on Thu, April 7th at 6:00 PM.

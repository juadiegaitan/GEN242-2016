---
title: Tutorial 1 - Introduction to Biocluster
last_updated: 29-Mar-16
---

## Hardware Infrastructure

### Compute nodes

- Over 3,000 CPU cores
- 50 batch compute nodes, each with 64 CPU cores and 512GB RAM
- 4 high-memory nodes, each 32 CPU cores and 1024GB RAM
- 4 GPU nodes, each with 10,000 cuda cores
    
### Interconnect 
- FDR IB @56Gbs 

### Storage

- Parallel GPFS storage system with 1.2PB usable space
- Backup of same architecture and similar amount

### User traffic

- Compute tasks need to be submitted via `qsub`
- Biocluster headnode only for login, not for compute tasks!
- Monitor cluster activity: `qstat` or `qstatMonitor`

### Manuals

- [Biocluster Manual](http://manuals.bioinformatics.ucr.edu/home/hpc)
- [Linux Manual](http://manuals.bioinformatics.ucr.edu/home/linux-basics#TOC-Exercises)


## Log into Biocluster

+ Login command on OS X or Linux 

{% highlight sh %}
ssh -X user@biocluster.ucr.edu
{% endhighlight %}
  
Type password

+ Windows: provide same information in [Putty](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html)
    
    + Host name: `biocluster.ucr.edu`
    + User name: ...
    + Password: ...

## Important Linux Commands

Finding help
{% highlight sh %}
man <program_name>
{% endhighlight %}
        
List content of current directory
{% highlight sh %}
ls
{% endhighlight %}

Print current working directory
{% highlight sh %}
pwd
{% endhighlight %}

Search in files and directories
{% highlight sh %}
grep
{% endhighlight %}

Word count
{% highlight sh %}
wc
{% endhighlight %}

Create directory
{% highlight sh %}
mkdir
{% endhighlight %}

Delete files and directories
{% highlight sh %}
rm
{% endhighlight %}

Move and rename files
{% highlight sh %}
mv
{% endhighlight %}

Copy files from internet to `pwd`
{% highlight sh %}
wget
{% endhighlight %}

Viewing files
{% highlight sh %}
less
{% endhighlight %}

Powerful text/code editors
{% highlight sh %}
vim or emacs
{% endhighlight %}

Less sophisticated text editors
{% highlight sh %}
pico or nano
{% endhighlight %}


## File Exchange

+ GUI applications
    + Windows: [WinSCP](http://winscp.net/eng/index.php)
	+ Mac OS X: [CyberDuck](http://cyberduck.en.softonic.com/mac)
	+ Win/OS X/Linux: [FileZilla](https://filezilla-project.org/)
        
+ SCP command line tool

{% highlight sh %}
scp file user@remotehost:/home/user/ # From local to remote 
scp user@remotehost:/home/user/file . # From remote to local 
{% endhighlight %}

	
## STD IN/OUT/ERR, Redirect & Wildcards
{% highlight sh %}
file.*                        # * wildcard to specify many files
ls > file                     # prints ls output into specified file
command < myfile              # uses file after '<' as STDIN
command >> myfile             # appends output of one command to file
command | tee myfile          # writes STDOUT to file and prints it to screen
command > myfile; cat myfile  # writes STDOUT to file and prints it to screen
command > /dev/null           # turns off progress info 
grep pattern file | wc        # Pipes (|) output of 'grep' into 'wc'
grep pattern nonexistingfile 2 > mystderr # prints STDERR to file
{% endhighlight %}


## Homework Assignment (HW2)

See HW2 page [here](http://girke.bioinformatics.ucr.edu/GEN242/mydoc/mydoc_homework_02.html).

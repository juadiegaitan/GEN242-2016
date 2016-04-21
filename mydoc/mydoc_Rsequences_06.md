---
title: Coerce `GRanges` to `RangedData` object and vice versa
keywords: 
last_updated: Wed Apr 20 19:57:09 2016
---

{% highlight r %}
gff_rd <- as(gff, "RangedData") 
gff_gr <- as(gff_rd, "GRanges") 
{% endhighlight %}

## Utilities for Range Containers

### Accessor and subsetting methods for GRanges objects


{% highlight r %}
gff[1:4]; gff[1:4, c("type", "ID")]; gff[2] <- gff[3] # Subsetting and replacement 
{% endhighlight %}

{% highlight txt %}
## GRanges object with 4 ranges and 10 metadata columns:
##     seqnames           ranges strand |   source       type     score     phase                  ID
##        <Rle>        <IRanges>  <Rle> | <factor>   <factor> <numeric> <integer>         <character>
##   1     Chr1 [   1, 30427671]      + |   TAIR10 chromosome      <NA>      <NA>                Chr1
##   2     Chr1 [3631,     5899]      + |   TAIR10       gene      <NA>      <NA>           AT1G01010
##   3     Chr1 [3631,     5899]      + |   TAIR10       mRNA      <NA>      <NA>         AT1G01010.1
##   4     Chr1 [3760,     5630]      + |   TAIR10    protein      <NA>      <NA> AT1G01010.1-Protein
##            Name                Note          Parent       Index Derives_from
##     <character>     <CharacterList> <CharacterList> <character>  <character>
##   1        Chr1                                            <NA>         <NA>
##   2   AT1G01010 protein_coding_gene                        <NA>         <NA>
##   3 AT1G01010.1                           AT1G01010           1         <NA>
##   4 AT1G01010.1                                            <NA>  AT1G01010.1
##   -------
##   seqinfo: 7 sequences from an unspecified genome
{% endhighlight %}

{% highlight txt %}
## GRanges object with 4 ranges and 2 metadata columns:
##     seqnames           ranges strand |       type                  ID
##        <Rle>        <IRanges>  <Rle> |   <factor>         <character>
##   1     Chr1 [   1, 30427671]      + | chromosome                Chr1
##   2     Chr1 [3631,     5899]      + |       gene           AT1G01010
##   3     Chr1 [3631,     5899]      + |       mRNA         AT1G01010.1
##   4     Chr1 [3760,     5630]      + |    protein AT1G01010.1-Protein
##   -------
##   seqinfo: 7 sequences from an unspecified genome
{% endhighlight %}

{% highlight r %}
c(gff[1:2], gff[401:402]) # GRanges objects can be concatenated with the c() function.
{% endhighlight %}

{% highlight txt %}
## GRanges object with 4 ranges and 10 metadata columns:
##       seqnames           ranges strand |   source           type     score     phase
##          <Rle>        <IRanges>  <Rle> | <factor>       <factor> <numeric> <integer>
##     1     Chr1 [   1, 30427671]      + |   TAIR10     chromosome      <NA>      <NA>
##     2     Chr1 [3631,     5899]      + |   TAIR10           mRNA      <NA>      <NA>
##   401     Chr5 [5516,     5769]      - |   TAIR10        protein      <NA>      <NA>
##   402     Chr5 [5770,     5801]      - |   TAIR10 five_prime_UTR      <NA>      <NA>
##                        ID        Name            Note          Parent       Index Derives_from
##               <character> <character> <CharacterList> <CharacterList> <character>  <character>
##     1                Chr1        Chr1                                        <NA>         <NA>
##     2         AT1G01010.1 AT1G01010.1                       AT1G01010           1         <NA>
##   401 AT5G01015.2-Protein AT5G01015.2                                        <NA>  AT5G01015.2
##   402                <NA>        <NA>                     AT5G01015.2        <NA>         <NA>
##   -------
##   seqinfo: 7 sequences from an unspecified genome
{% endhighlight %}

{% highlight r %}
seqnames(gff); ranges(gff); strand(gff); seqlengths(gff) # Accessor functions 
{% endhighlight %}

{% highlight txt %}
## factor-Rle of length 449 with 7 runs
##   Lengths:   72   22   38  118  172   13   14
##   Values : Chr1 Chr2 Chr3 Chr4 Chr5 ChrC ChrM
## Levels(7): Chr1 Chr2 Chr3 Chr4 Chr5 ChrC ChrM
{% endhighlight %}

{% highlight txt %}
## IRanges of length 449
##       start      end    width names
## [1]       1 30427671 30427671     1
## [2]    3631     5899     2269     2
## [3]    3631     5899     2269     3
## [4]    3760     5630     1871     4
## [5]    3631     3913      283     5
## ...     ...      ...      ...   ...
## [445] 11918    12241      324   445
## [446] 11918    12241      324   446
## [447] 11918    12241      324   447
## [448] 11918    12241      324   448
## [449] 11918    12241      324   449
{% endhighlight %}

{% highlight txt %}
## factor-Rle of length 449 with 13 runs
##   Lengths:  18  54  28  21  12 117   1 171   1  12   1   8   5
##   Values :   +   -   +   -   +   -   +   -   +   -   +   -   +
## Levels(3): + - *
{% endhighlight %}

{% highlight txt %}
##     Chr1     Chr2     Chr3     Chr4     Chr5     ChrC     ChrM 
## 30427671 19698289 23459830 18585056 26975502   154478   366924
{% endhighlight %}

{% highlight r %}
start(gff[1:4]); end(gff[1:4]); width(gff[1:4]) # Direct access to IRanges components
{% endhighlight %}

{% highlight txt %}
## [1]    1 3631 3631 3760
{% endhighlight %}

{% highlight txt %}
## [1] 30427671     5899     5899     5630
{% endhighlight %}

{% highlight txt %}
## [1] 30427671     2269     2269     1871
{% endhighlight %}

{% highlight r %}
values(gff); values(gff)[, "type"] # Accessing metadata component. 
{% endhighlight %}

{% highlight txt %}
## DataFrame with 449 rows and 10 columns
##       source       type     score     phase                  ID        Name                Note
##     <factor>   <factor> <numeric> <integer>         <character> <character>     <CharacterList>
## 1     TAIR10 chromosome        NA        NA                Chr1        Chr1                    
## 2     TAIR10       mRNA        NA        NA         AT1G01010.1 AT1G01010.1                    
## 3     TAIR10       mRNA        NA        NA         AT1G01010.1 AT1G01010.1                    
## 4     TAIR10    protein        NA        NA AT1G01010.1-Protein AT1G01010.1                    
## 5     TAIR10       exon        NA        NA                  NA          NA                    
## ...      ...        ...       ...       ...                 ...         ...                 ...
## 445   TAIR10       gene        NA        NA           ATMG00030   ATMG00030 protein_coding_gene
## 446   TAIR10       mRNA        NA        NA         ATMG00030.1 ATMG00030.1                    
## 447   TAIR10    protein        NA        NA ATMG00030.1-Protein ATMG00030.1                    
## 448   TAIR10       exon        NA        NA                  NA          NA                    
## 449   TAIR10        CDS        NA         0                  NA          NA                    
##                              Parent       Index Derives_from
##                     <CharacterList> <character>  <character>
## 1                                            NA           NA
## 2                         AT1G01010           1           NA
## 3                         AT1G01010           1           NA
## 4                                            NA  AT1G01010.1
## 5                       AT1G01010.1          NA           NA
## ...                             ...         ...          ...
## 445                                          NA           NA
## 446                       ATMG00030           1           NA
## 447                                          NA  ATMG00030.1
## 448                     ATMG00030.1          NA           NA
## 449 ATMG00030.1,ATMG00030.1-Protein          NA           NA
{% endhighlight %}

{% highlight txt %}
##   [1] chromosome      mRNA            mRNA            protein         exon           
##   [6] five_prime_UTR  CDS             exon            CDS             exon           
##  [11] CDS             exon            CDS             exon            CDS            
##  [16] exon            CDS             three_prime_UTR gene            mRNA           
##  [21] protein         five_prime_UTR  CDS             exon            CDS            
##  [26] exon            CDS             exon            CDS             exon           
##  [31] CDS             exon            CDS             exon            CDS            
##  [36] exon            CDS             exon            CDS             three_prime_UTR
##  [41] exon            three_prime_UTR exon            mRNA            protein        
##  [46] five_prime_UTR  CDS             exon            CDS             exon           
##  [51] CDS             exon            CDS             exon            CDS            
##  [56] exon            CDS             exon            CDS             three_prime_UTR
##  [61] exon            three_prime_UTR exon            gene            mRNA           
##  [66] protein         five_prime_UTR  exon            five_prime_UTR  CDS            
##  [71] three_prime_UTR exon            chromosome      gene            mRNA           
##  [76] protein         exon            CDS             exon            CDS            
##  [81] exon            CDS             three_prime_UTR gene            rRNA           
##  [86] exon            gene            rRNA            exon            gene           
##  [91] mRNA            protein         exon            CDS             chromosome     
##  [96] gene            mRNA            protein         exon            CDS            
## [101] gene            mRNA            protein         five_prime_UTR  CDS            
## [106] exon            CDS             exon            CDS             exon           
## [111] CDS             exon            CDS             exon            CDS            
## [116] exon            CDS             exon            CDS             three_prime_UTR
## [121] exon            gene            mRNA            protein         exon           
## [126] five_prime_UTR  CDS             exon            CDS             exon           
## [131] CDS             three_prime_UTR chromosome      gene            mRNA           
## [136] protein         CDS             exon            gene            mRNA           
## [141] protein         five_prime_UTR  CDS             exon            CDS            
## [146] exon            CDS             exon            CDS             exon           
## [151] CDS             exon            CDS             exon            CDS            
## [156] exon            CDS             exon            CDS             exon           
## [161] CDS             exon            CDS             exon            CDS            
## [166] exon            CDS             exon            CDS             exon           
## [171] CDS             exon            CDS             exon            CDS            
## [176] exon            CDS             exon            CDS             exon           
## [181] CDS             three_prime_UTR exon            mRNA            protein        
## [186] CDS             exon            CDS             exon            CDS            
## [191] exon            CDS             exon            CDS             exon           
## [196] CDS             exon            CDS             exon            CDS            
## [201] exon            CDS             exon            CDS             exon           
## [206] CDS             exon            CDS             exon            CDS            
## [211] exon            CDS             exon            CDS             exon           
## [216] CDS             exon            CDS             exon            CDS            
## [221] exon            CDS             exon            CDS             exon           
## [226] CDS             exon            CDS             exon            gene           
## [231] mRNA            protein         five_prime_UTR  CDS             exon           
## [236] CDS             exon            CDS             exon            CDS            
## [241] exon            CDS             exon            CDS             exon           
## [246] CDS             exon            CDS             three_prime_UTR exon           
## [251] chromosome      gene            mRNA            protein         five_prime_UTR 
## [256] CDS             exon            CDS             exon            CDS            
## [261] exon            CDS             exon            CDS             exon           
## [266] CDS             exon            CDS             exon            CDS            
## [271] exon            CDS             exon            CDS             exon           
## [276] CDS             exon            CDS             exon            CDS            
## [281] exon            CDS             exon            CDS             three_prime_UTR
## [286] exon            mRNA            protein         five_prime_UTR  CDS            
## [291] exon            CDS             exon            CDS             exon           
## [296] CDS             exon            CDS             exon            CDS            
## [301] exon            CDS             exon            CDS             exon           
## [306] CDS             exon            CDS             exon            CDS            
## [311] exon            CDS             exon            CDS             exon           
## [316] CDS             exon            CDS             exon            CDS            
## [321] three_prime_UTR exon            mRNA            protein         five_prime_UTR 
## [326] CDS             exon            CDS             exon            CDS            
## [331] exon            CDS             exon            CDS             exon           
## [336] CDS             exon            CDS             exon            CDS            
## [341] exon            CDS             exon            CDS             exon           
## [346] CDS             exon            CDS             exon            CDS            
## [351] exon            CDS             three_prime_UTR exon            mRNA           
## [356] protein         five_prime_UTR  CDS             exon            CDS            
## [361] exon            CDS             exon            CDS             exon           
## [366] CDS             exon            CDS             exon            CDS            
## [371] exon            CDS             exon            CDS             exon           
## [376] CDS             exon            CDS             exon            CDS            
## [381] exon            CDS             exon            CDS             exon           
## [386] CDS             exon            CDS             three_prime_UTR exon           
## [391] gene            mRNA            protein         five_prime_UTR  CDS            
## [396] exon            CDS             three_prime_UTR exon            mRNA           
## [401] protein         five_prime_UTR  CDS             exon            CDS            
## [406] three_prime_UTR exon            gene            mRNA            protein        
## [411] five_prime_UTR  CDS             exon            CDS             exon           
## [416] CDS             exon            CDS             exon            CDS            
## [421] three_prime_UTR exon            chromosome      gene            tRNA           
## [426] exon            gene            mRNA            protein         CDS            
## [431] exon            gene            tRNA            exon            exon           
## [436] chromosome      gene            mRNA            protein         CDS            
## [441] exon            gene            rRNA            exon            gene           
## [446] mRNA            protein         exon            CDS            
## Levels: chromosome gene mRNA protein exon five_prime_UTR CDS three_prime_UTR rRNA tRNA
{% endhighlight %}

{% highlight r %}
gff[elementMetadata(gff)[ ,"type"] == "gene"] # Returns only gene ranges.
{% endhighlight %}

{% highlight txt %}
## GRanges object with 21 ranges and 10 metadata columns:
##       seqnames         ranges strand   |   source     type     score     phase          ID
##          <Rle>      <IRanges>  <Rle>   | <factor> <factor> <numeric> <integer> <character>
##    19     Chr1 [ 5928,  8737]      -   |   TAIR10     gene      <NA>      <NA>   AT1G01020
##    64     Chr1 [11649, 13714]      -   |   TAIR10     gene      <NA>      <NA>   AT1G01030
##    74     Chr2 [ 1025,  2810]      +   |   TAIR10     gene      <NA>      <NA>   AT2G01008
##    84     Chr2 [ 3706,  5513]      +   |   TAIR10     gene      <NA>      <NA>   AT2G01010
##    87     Chr2 [ 5782,  5945]      +   |   TAIR10     gene      <NA>      <NA>   AT2G01020
##   ...      ...            ...    ... ...      ...      ...       ...       ...         ...
##   427     ChrC [  383,  1444]      -   |   TAIR10     gene      <NA>      <NA>   ATCG00020
##   432     ChrC [ 1717,  4347]      -   |   TAIR10     gene      <NA>      <NA>   ATCG00030
##   437     ChrM [  273,   734]      -   |   TAIR10     gene      <NA>      <NA>   ATMG00010
##   442     ChrM [ 8848, 11415]      -   |   TAIR10     gene      <NA>      <NA>   ATMG00020
##   445     ChrM [11918, 12241]      +   |   TAIR10     gene      <NA>      <NA>   ATMG00030
##              Name                Note          Parent       Index Derives_from
##       <character>     <CharacterList> <CharacterList> <character>  <character>
##    19   AT1G01020 protein_coding_gene                        <NA>         <NA>
##    64   AT1G01030 protein_coding_gene                        <NA>         <NA>
##    74   AT2G01008 protein_coding_gene                        <NA>         <NA>
##    84   AT2G01010                rRNA                        <NA>         <NA>
##    87   AT2G01020                rRNA                        <NA>         <NA>
##   ...         ...                 ...             ...         ...          ...
##   427   ATCG00020 protein_coding_gene                        <NA>         <NA>
##   432   ATCG00030                tRNA                        <NA>         <NA>
##   437   ATMG00010 protein_coding_gene                        <NA>         <NA>
##   442   ATMG00020                rRNA                        <NA>         <NA>
##   445   ATMG00030 protein_coding_gene                        <NA>         <NA>
##   -------
##   seqinfo: 7 sequences from an unspecified genome
{% endhighlight %}

### Useful utilities for GRanges objects


{% highlight r %}
gff <- gff[values(gff)$type != "chromosome"] # Remove chromosome ranges
strand(gff) <- "*" # Erases the strand information
reduce(gff) # Collapses overlapping ranges to continuous ranges.
{% endhighlight %}

{% highlight txt %}
## GRanges object with 22 ranges and 0 metadata columns:
##        seqnames         ranges strand
##           <Rle>      <IRanges>  <Rle>
##    [1]     Chr1 [ 3631,  5899]      *
##    [2]     Chr1 [ 5928,  8737]      *
##    [3]     Chr1 [11649, 13714]      *
##    [4]     Chr2 [ 1025,  2810]      *
##    [5]     Chr2 [ 3706,  5513]      *
##    ...      ...            ...    ...
##   [18]     ChrC [  383,  1444]      *
##   [19]     ChrC [ 1717,  4347]      *
##   [20]     ChrM [  273,   734]      *
##   [21]     ChrM [ 8848, 11415]      *
##   [22]     ChrM [11918, 12241]      *
##   -------
##   seqinfo: 7 sequences from an unspecified genome
{% endhighlight %}

{% highlight r %}
gaps(gff) # Returns uncovered regions.
{% endhighlight %}

{% highlight txt %}
## GRanges object with 43 ranges and 0 metadata columns:
##        seqnames           ranges strand
##           <Rle>        <IRanges>  <Rle>
##    [1]     Chr1 [   1, 30427671]      +
##    [2]     Chr1 [   1, 30427671]      -
##    [3]     Chr1 [   1,     3630]      *
##    [4]     Chr1 [5900,     5927]      *
##    [5]     Chr1 [8738,    11648]      *
##    ...      ...              ...    ...
##   [39]     ChrM  [    1, 366924]      -
##   [40]     ChrM  [    1,    272]      *
##   [41]     ChrM  [  735,   8847]      *
##   [42]     ChrM  [11416,  11917]      *
##   [43]     ChrM  [12242, 366924]      *
##   -------
##   seqinfo: 7 sequences from an unspecified genome
{% endhighlight %}

{% highlight r %}
disjoin(gff) # Returns disjoint ranges.
{% endhighlight %}

{% highlight txt %}
## GRanges object with 211 ranges and 0 metadata columns:
##         seqnames         ranges strand
##            <Rle>      <IRanges>  <Rle>
##     [1]     Chr1   [3631, 3759]      *
##     [2]     Chr1   [3760, 3913]      *
##     [3]     Chr1   [3914, 3995]      *
##     [4]     Chr1   [3996, 4276]      *
##     [5]     Chr1   [4277, 4485]      *
##     ...      ...            ...    ...
##   [207]     ChrC [ 1752,  4310]      *
##   [208]     ChrC [ 4311,  4347]      *
##   [209]     ChrM [  273,   734]      *
##   [210]     ChrM [ 8848, 11415]      *
##   [211]     ChrM [11918, 12241]      *
##   -------
##   seqinfo: 7 sequences from an unspecified genome
{% endhighlight %}

{% highlight r %}
coverage(gff) # Returns coverage of ranges.
{% endhighlight %}

{% highlight txt %}
## RleList of length 7
## $Chr1
## integer-Rle of length 30427671 with 45 runs
##   Lengths:     3630      129      154       82      281 ...      233      161      380 30413957
##   Values :        0        4        5        3        5 ...        4        2        4        0
## 
## $Chr2
## integer-Rle of length 19698289 with 14 runs
##   Lengths:     1024      248      185       53      362 ...      164      625      102 19691617
##   Values :        0        5        3        5        3 ...        3        0        5        0
## 
## $Chr3
## integer-Rle of length 23459830 with 29 runs
##   Lengths:     1652      145      139      111       95 ...      155      148      156 23453781
##   Values :        0        4        5        3        5 ...        3        5        4        0
## 
## $Chr4
## integer-Rle of length 18585056 with 72 runs
##   Lengths:     1179      357     1358      128      872 ...      212      114       74 18571697
##   Values :        0        5        0        5        3 ...        3        5        4        0
## 
## $Chr5
## integer-Rle of length 26975502 with 64 runs
##   Lengths:     1222       28       28      109       72 ...       76       55      174 26967058
##   Values :        0        4        7       13       16 ...        3        5        4        0
## 
## ...
## <2 more elements>
{% endhighlight %}

{% highlight r %}
findOverlaps(gff, gff[1:4]) # Returns the index pairings for the overlapping ranges. 
{% endhighlight %}

{% highlight txt %}
## Hits object with 55 hits and 0 metadata columns:
##        queryHits subjectHits
##        <integer>   <integer>
##    [1]         1           1
##    [2]         1           2
##    [3]         1           4
##    [4]         1           3
##    [5]         2           1
##    ...       ...         ...
##   [51]        16           1
##   [52]        16           2
##   [53]        16           3
##   [54]        17           1
##   [55]        17           2
##   -------
##   queryLength: 442
##   subjectLength: 4
{% endhighlight %}

{% highlight r %}
countOverlaps(gff, gff[1:4]) # Counts overlapping ranges 
{% endhighlight %}

{% highlight txt %}
##   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26 
##   4   4   4   4   3   4   3   3   3   3   3   3   3   3   3   3   2   0   0   0   0   0   0   0   0 
##  27  28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
##  52  53  54  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72  74  75  76  77 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
##  78  79  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  96  97  98  99 100 101 102 103 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
## 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
## 129 130 131 132 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
## 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
## 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
## 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
## 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 252 253 254 255 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
## 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275 276 277 278 279 280 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
## 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
## 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329 330 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
## 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346 347 348 349 350 351 352 353 354 355 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
## 356 357 358 359 360 361 362 363 364 365 366 367 368 369 370 371 372 373 374 375 376 377 378 379 380 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
## 381 382 383 384 385 386 387 388 389 390 391 392 393 394 395 396 397 398 399 400 401 402 403 404 405 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
## 406 407 408 409 410 411 412 413 414 415 416 417 418 419 420 421 422 424 425 426 427 428 429 430 431 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
## 432 433 434 435 437 438 439 440 441 442 443 444 445 446 447 448 449 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
{% endhighlight %}

{% highlight r %}
subsetByOverlaps(gff, gff[1:4]) # Returns only overlapping ranges 
{% endhighlight %}

{% highlight txt %}
## GRanges object with 17 ranges and 10 metadata columns:
##       seqnames       ranges strand   |   source            type     score     phase
##          <Rle>    <IRanges>  <Rle>   | <factor>        <factor> <numeric> <integer>
##     2     Chr1 [3631, 5899]      *   |   TAIR10            mRNA      <NA>      <NA>
##     3     Chr1 [3631, 5899]      *   |   TAIR10            mRNA      <NA>      <NA>
##     4     Chr1 [3760, 5630]      *   |   TAIR10         protein      <NA>      <NA>
##     5     Chr1 [3631, 3913]      *   |   TAIR10            exon      <NA>      <NA>
##     6     Chr1 [3631, 3759]      *   |   TAIR10  five_prime_UTR      <NA>      <NA>
##   ...      ...          ...    ... ...      ...             ...       ...       ...
##    14     Chr1 [5174, 5326]      *   |   TAIR10            exon      <NA>      <NA>
##    15     Chr1 [5174, 5326]      *   |   TAIR10             CDS      <NA>         0
##    16     Chr1 [5439, 5899]      *   |   TAIR10            exon      <NA>      <NA>
##    17     Chr1 [5439, 5630]      *   |   TAIR10             CDS      <NA>         0
##    18     Chr1 [5631, 5899]      *   |   TAIR10 three_prime_UTR      <NA>      <NA>
##                        ID        Name            Note                          Parent       Index
##               <character> <character> <CharacterList>                 <CharacterList> <character>
##     2         AT1G01010.1 AT1G01010.1                                       AT1G01010           1
##     3         AT1G01010.1 AT1G01010.1                                       AT1G01010           1
##     4 AT1G01010.1-Protein AT1G01010.1                                                        <NA>
##     5                <NA>        <NA>                                     AT1G01010.1        <NA>
##     6                <NA>        <NA>                                     AT1G01010.1        <NA>
##   ...                 ...         ...             ...                             ...         ...
##    14                <NA>        <NA>                                     AT1G01010.1        <NA>
##    15                <NA>        <NA>                 AT1G01010.1,AT1G01010.1-Protein        <NA>
##    16                <NA>        <NA>                                     AT1G01010.1        <NA>
##    17                <NA>        <NA>                 AT1G01010.1,AT1G01010.1-Protein        <NA>
##    18                <NA>        <NA>                                     AT1G01010.1        <NA>
##       Derives_from
##        <character>
##     2         <NA>
##     3         <NA>
##     4  AT1G01010.1
##     5         <NA>
##     6         <NA>
##   ...          ...
##    14         <NA>
##    15         <NA>
##    16         <NA>
##    17         <NA>
##    18         <NA>
##   -------
##   seqinfo: 7 sequences from an unspecified genome
{% endhighlight %}

## GRangesList Objects


{% highlight r %}
sp <- split(gff, seq(along=gff)) # Stores every range in separate component of a GRangesList object
split(gff, seqnames(gff)) # Stores ranges of each chromosome in separate component.
{% endhighlight %}

{% highlight txt %}
## GRangesList object of length 7:
## $Chr1 
## GRanges object with 71 ranges and 10 metadata columns:
##       seqnames         ranges strand   |   source            type     score     phase
##          <Rle>      <IRanges>  <Rle>   | <factor>        <factor> <numeric> <integer>
##     2     Chr1   [3631, 5899]      *   |   TAIR10            mRNA      <NA>      <NA>
##     3     Chr1   [3631, 5899]      *   |   TAIR10            mRNA      <NA>      <NA>
##     4     Chr1   [3760, 5630]      *   |   TAIR10         protein      <NA>      <NA>
##     5     Chr1   [3631, 3913]      *   |   TAIR10            exon      <NA>      <NA>
##     6     Chr1   [3631, 3759]      *   |   TAIR10  five_prime_UTR      <NA>      <NA>
##   ...      ...            ...    ... ...      ...             ...       ...       ...
##    68     Chr1 [13335, 13714]      *   |   TAIR10            exon      <NA>      <NA>
##    69     Chr1 [12941, 13173]      *   |   TAIR10  five_prime_UTR      <NA>      <NA>
##    70     Chr1 [11864, 12940]      *   |   TAIR10             CDS      <NA>         0
##    71     Chr1 [11649, 11863]      *   |   TAIR10 three_prime_UTR      <NA>      <NA>
##    72     Chr1 [11649, 13173]      *   |   TAIR10            exon      <NA>      <NA>
##                        ID        Name            Note                          Parent       Index
##               <character> <character> <CharacterList>                 <CharacterList> <character>
##     2         AT1G01010.1 AT1G01010.1                                       AT1G01010           1
##     3         AT1G01010.1 AT1G01010.1                                       AT1G01010           1
##     4 AT1G01010.1-Protein AT1G01010.1                                                        <NA>
##     5                <NA>        <NA>                                     AT1G01010.1        <NA>
##     6                <NA>        <NA>                                     AT1G01010.1        <NA>
##   ...                 ...         ...             ...                             ...         ...
##    68                <NA>        <NA>                                     AT1G01030.1        <NA>
##    69                <NA>        <NA>                                     AT1G01030.1        <NA>
##    70                <NA>        <NA>                 AT1G01030.1,AT1G01030.1-Protein        <NA>
##    71                <NA>        <NA>                                     AT1G01030.1        <NA>
##    72                <NA>        <NA>                                     AT1G01030.1        <NA>
##       Derives_from
##        <character>
##     2         <NA>
##     3         <NA>
##     4  AT1G01010.1
##     5         <NA>
##     6         <NA>
##   ...          ...
##    68         <NA>
##    69         <NA>
##    70         <NA>
##    71         <NA>
##    72         <NA>
## 
## ...
## <6 more elements>
## -------
## seqinfo: 7 sequences from an unspecified genome
{% endhighlight %}

{% highlight r %}
unlist(sp) # Returns data as GRanges object
{% endhighlight %}

{% highlight txt %}
## GRanges object with 442 ranges and 10 metadata columns:
##           seqnames         ranges strand   |   source           type     score     phase
##              <Rle>      <IRanges>  <Rle>   | <factor>       <factor> <numeric> <integer>
##       1.2     Chr1   [3631, 5899]      *   |   TAIR10           mRNA      <NA>      <NA>
##       2.3     Chr1   [3631, 5899]      *   |   TAIR10           mRNA      <NA>      <NA>
##       3.4     Chr1   [3760, 5630]      *   |   TAIR10        protein      <NA>      <NA>
##       4.5     Chr1   [3631, 3913]      *   |   TAIR10           exon      <NA>      <NA>
##       5.6     Chr1   [3631, 3759]      *   |   TAIR10 five_prime_UTR      <NA>      <NA>
##       ...      ...            ...    ... ...      ...            ...       ...       ...
##   438.445     ChrM [11918, 12241]      *   |   TAIR10           gene      <NA>      <NA>
##   439.446     ChrM [11918, 12241]      *   |   TAIR10           mRNA      <NA>      <NA>
##   440.447     ChrM [11918, 12241]      *   |   TAIR10        protein      <NA>      <NA>
##   441.448     ChrM [11918, 12241]      *   |   TAIR10           exon      <NA>      <NA>
##   442.449     ChrM [11918, 12241]      *   |   TAIR10            CDS      <NA>         0
##                            ID        Name                Note                          Parent
##                   <character> <character>     <CharacterList>                 <CharacterList>
##       1.2         AT1G01010.1 AT1G01010.1                                           AT1G01010
##       2.3         AT1G01010.1 AT1G01010.1                                           AT1G01010
##       3.4 AT1G01010.1-Protein AT1G01010.1                                                    
##       4.5                <NA>        <NA>                                         AT1G01010.1
##       5.6                <NA>        <NA>                                         AT1G01010.1
##       ...                 ...         ...                 ...                             ...
##   438.445           ATMG00030   ATMG00030 protein_coding_gene                                
##   439.446         ATMG00030.1 ATMG00030.1                                           ATMG00030
##   440.447 ATMG00030.1-Protein ATMG00030.1                                                    
##   441.448                <NA>        <NA>                                         ATMG00030.1
##   442.449                <NA>        <NA>                     ATMG00030.1,ATMG00030.1-Protein
##                 Index Derives_from
##           <character>  <character>
##       1.2           1         <NA>
##       2.3           1         <NA>
##       3.4        <NA>  AT1G01010.1
##       4.5        <NA>         <NA>
##       5.6        <NA>         <NA>
##       ...         ...          ...
##   438.445        <NA>         <NA>
##   439.446           1         <NA>
##   440.447        <NA>  ATMG00030.1
##   441.448        <NA>         <NA>
##   442.449        <NA>         <NA>
##   -------
##   seqinfo: 7 sequences from an unspecified genome
{% endhighlight %}

{% highlight r %}
sp[1:4, "type"] # Subsetting of GRangesList objects is similar to GRanges objects.
{% endhighlight %}

{% highlight txt %}
## GRangesList object of length 4:
## $1 
## GRanges object with 1 range and 1 metadata column:
##     seqnames       ranges strand |     type
##        <Rle>    <IRanges>  <Rle> | <factor>
##   2     Chr1 [3631, 5899]      * |     mRNA
## 
## $2 
## GRanges object with 1 range and 1 metadata column:
##     seqnames       ranges strand | type
##   3     Chr1 [3631, 5899]      * | mRNA
## 
## $3 
## GRanges object with 1 range and 1 metadata column:
##     seqnames       ranges strand |    type
##   4     Chr1 [3760, 5630]      * | protein
## 
## ...
## <1 more element>
## -------
## seqinfo: 7 sequences from an unspecified genome
{% endhighlight %}

{% highlight r %}
lapply(sp[1:4], length); sapply(sp[1:4], length) # Looping over GRangesList objects similar to lists
{% endhighlight %}

{% highlight txt %}
## $`1`
## [1] 1
## 
## $`2`
## [1] 1
## 
## $`3`
## [1] 1
## 
## $`4`
## [1] 1
{% endhighlight %}

{% highlight txt %}
## 1 2 3 4 
## 1 1 1 1
{% endhighlight %}


---
title: VAR-Seq1 - Comparison of Variant Callers
last_updated: 29-Apr-16
---

## VAR-Seq Workflow  

+ Read preprocessing
    + Quality filtering (trimming)
    + FASTQ quality report
+ Alignments
+ Alignment statistics
+ Variant calling
+ Variant filtering
+ Variant annotation
+ Combine results from many samples
+ Summary statistics of samples

## Challenge Project: Comparison of variant calling methods

+ Run workflow from start to finish (steps 1-8) on data set from Lu et al (2012)
+ Challenge project tasks
    + Compare performance of at least 2 variant callers, e.g. GATK, BCFtools and VariantTools. Include in your comparisons the following analysis/visualization steps:
        + Report unique and common variants identified by tested variant callers.
        + Compare the results from (a) with the variants identified by Lu et al, 2012
        + Plot results from a-b as venn diagrams
        + Plot the performance of the variant callers in form of ROC plots. As true result set one can use the union of the variants identified by all methods. 

## References

+ Lu P, Han X, Qi J, Yang J, Wijeratne AJ, Li T, Ma H (2012) Analysis of Arabidopsis genome-wide variations before and after meiosis and meiotic recombination by resequencing Landsberg erecta and all four products of a single meiosis. Genome Res 22: 508–518 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/22106370)
+ Kumar RD, Chang L-W, Ellis MJ, Bose R (2013) Prioritizing Potentially Druggable Mutations with dGene: An Annotation Tool for Cancer Genome Sequencing Data. PLoS One 8: e67980 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/23826350)
+ DePristo MA, Banks E, Poplin R, Garimella KV, Maguire JR, Hartl C, Philippakis AA, del Angel G, Rivas MA, Hanna M, et al (2011) A framework for variation discovery and genotyping using next-generation DNA sequencing data. Nat Genet 43: 491–498 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/21478889)
+ Li H (2011) A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics 27: 2987–2993 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/21903627)
+ McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, et al (2010) The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res 20: 1297–1303 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/20644199)




---
title: "Practicing R - ICGC Data Portal"
author: "Victoria Rios"
date: "3/7/2021"
output: html_document
---



## R Markdown
Data set source: [https://dcc.icgc.org/releases/release_28/Projects/BRCA-US]


```r
x <- read.csv('exp_seq.BRCA-US.tsv.gz',sep = '\t', nrows = 2000, na.strings = '?')
str(x)
```

```
## 'data.frame':	2000 obs. of  22 variables:
##  $ icgc_donor_id           : chr  "DO1808" "DO1808" "DO1808" "DO1808" ...
##  $ project_code            : chr  "BRCA-US" "BRCA-US" "BRCA-US" "BRCA-US" ...
##  $ icgc_specimen_id        : chr  "SP3930" "SP3930" "SP3930" "SP3930" ...
##  $ icgc_sample_id          : chr  "SA42163" "SA42163" "SA42163" "SA42163" ...
##  $ submitted_sample_id     : chr  "TCGA-E2-A15C-01A-31R-A12D-07" "TCGA-E2-A15C-01A-31R-A12D-07" "TCGA-E2-A15C-01A-31R-A12D-07" "TCGA-E2-A15C-01A-31R-A12D-07" ...
##  $ analysis_id             : int  7768 7768 7768 7768 7768 7768 7768 7768 7768 7768 ...
##  $ gene_model              : chr  "GAF" "GAF" "GAF" "GAF" ...
##  $ gene_id                 : chr  "ALPL" "ALPK1" "ALPK2" "ALPK3" ...
##  $ normalized_read_count   : num  4.49e-06 2.03e-05 7.66e-07 3.50e-06 0.00 ...
##  $ raw_read_count          : int  278 2189 140 607 0 0 347 20 17 1291 ...
##  $ fold_change             : logi  NA NA NA NA NA NA ...
##  $ assembly_version        : chr  "GRCh37" "GRCh37" "GRCh37" "GRCh37" ...
##  $ platform                : chr  "Illumina HiSeq" "Illumina HiSeq" "Illumina HiSeq" "Illumina HiSeq" ...
##  $ total_read_count        : logi  NA NA NA NA NA NA ...
##  $ experimental_protocol   : chr  "RNASeqV2_RSEM_genes https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor" "RNASeqV2_RSEM_genes https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor" "RNASeqV2_RSEM_genes https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor" "RNASeqV2_RSEM_genes https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor" ...
##  $ alignment_algorithm     : logi  NA NA NA NA NA NA ...
##  $ normalization_algorithm : logi  NA NA NA NA NA NA ...
##  $ other_analysis_algorithm: logi  NA NA NA NA NA NA ...
##  $ sequencing_strategy     : chr  "RNA-Seq" "RNA-Seq" "RNA-Seq" "RNA-Seq" ...
##  $ raw_data_repository     : chr  "CGHub" "CGHub" "CGHub" "CGHub" ...
##  $ raw_data_accession      : chr  "1a77082d-f130-436b-9028-23c8dd9dc565" "1a77082d-f130-436b-9028-23c8dd9dc565" "1a77082d-f130-436b-9028-23c8dd9dc565" "1a77082d-f130-436b-9028-23c8dd9dc565" ...
##  $ reference_sample_type   : logi  NA NA NA NA NA NA ...
```

Dropping the columns with NA values or not relevant information

```r
x <- subset(x, select = -c(fold_change,total_read_count,platform,raw_data_repository,experimental_protocol,alignment_algorithm,normalization_algorithm,other_analysis_algorithm,reference_sample_type))

x[is.na(x)] <- 0 # replacing NA values with a 0

str(x)
```

```
## 'data.frame':	2000 obs. of  13 variables:
##  $ icgc_donor_id        : chr  "DO1808" "DO1808" "DO1808" "DO1808" ...
##  $ project_code         : chr  "BRCA-US" "BRCA-US" "BRCA-US" "BRCA-US" ...
##  $ icgc_specimen_id     : chr  "SP3930" "SP3930" "SP3930" "SP3930" ...
##  $ icgc_sample_id       : chr  "SA42163" "SA42163" "SA42163" "SA42163" ...
##  $ submitted_sample_id  : chr  "TCGA-E2-A15C-01A-31R-A12D-07" "TCGA-E2-A15C-01A-31R-A12D-07" "TCGA-E2-A15C-01A-31R-A12D-07" "TCGA-E2-A15C-01A-31R-A12D-07" ...
##  $ analysis_id          : int  7768 7768 7768 7768 7768 7768 7768 7768 7768 7768 ...
##  $ gene_model           : chr  "GAF" "GAF" "GAF" "GAF" ...
##  $ gene_id              : chr  "ALPL" "ALPK1" "ALPK2" "ALPK3" ...
##  $ normalized_read_count: num  4.49e-06 2.03e-05 7.66e-07 3.50e-06 0.00 ...
##  $ raw_read_count       : int  278 2189 140 607 0 0 347 20 17 1291 ...
##  $ assembly_version     : chr  "GRCh37" "GRCh37" "GRCh37" "GRCh37" ...
##  $ sequencing_strategy  : chr  "RNA-Seq" "RNA-Seq" "RNA-Seq" "RNA-Seq" ...
##  $ raw_data_accession   : chr  "1a77082d-f130-436b-9028-23c8dd9dc565" "1a77082d-f130-436b-9028-23c8dd9dc565" "1a77082d-f130-436b-9028-23c8dd9dc565" "1a77082d-f130-436b-9028-23c8dd9dc565" ...
```
## Clustering
It's time to build the distance matrix. We use the euclidean distance method.

```r
d <- dist(x, method = 'euclidean')
```

```
## Warning in dist(x, method = "euclidean"): NAs introduced by coercion
```

## Hierarchial clustering
We proceed with average linkage method.

```r
h <- hclust(d, method = 'average')
```


Looking at the dendogram

```r
plot(h)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

Next, we can cut the dendrogram in order to create the desired number of clusters.

```r
cut_avg <- cutree(h, k=10)
```


```r
library(dplyr)

plot(h)
rect.hclust(h , k = 10, border = 2:6)
abline(h = 10, col = 'red')
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)
Now we can see the three clusters enclosed in three different colored boxes. 

```r
suppressPackageStartupMessages(library(dendextend))
avg_dend_obj <- as.dendrogram(h)
avg_col_dend <- color_branches(h, h = 10)
plot(avg_col_dend)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)
Now we  append the cluster results obtained back in the original dataframe under column name the cluster with mutate(), from the dplyr package and count how many observations were assigned to each cluster with the count() function.

```r
suppressPackageStartupMessages(library(dplyr))
x_clusters <- mutate(x, cluster = cut_avg)
count(x_clusters,cluster)
```

```
##    cluster    n
## 1        1 1893
## 2        2   82
## 3        3   11
## 4        4    7
## 5        5    2
## 6        6    1
## 7        7    1
## 8        8    1
## 9        9    1
## 10      10    1
```
Visualizing some relationships of the data

```r
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))

ggplotly(ggplot(x_clusters, aes(x=gene_id, y = normalized_read_count, color = factor(cluster))) + geom_point())
```

```
## PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.
```

```
## Error in path.expand(path): invalid 'path' argument
```




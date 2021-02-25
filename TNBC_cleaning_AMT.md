---
title: "RNASeq"
author:
  - "Victoria Ríos"
  - "Antonio Martínez"
date: "February 24, 2021"
output:
  html_document:
    keep_md: true
---



## Import RNASeq data frame

```r
data <- read.csv("TNBC__maseqv2.csv",stringsAsFactors=FALSE,header=FALSE,row.names=1)
colnames(data) <- data[2,] #Define column names
data <- data[-c(1:3),] #Remove trash
df.RNAseq <- data.frame(lapply(data,as.numeric),row.names=rownames(data)) #Make all columns numeric vectors

class(df.RNAseq) # Verify its class
```

```
## [1] "data.frame"
```

```r
str(df.RNAseq,list.len=10) # Check its structure
```

```
## 'data.frame':	20531 obs. of  1212 variables:
##  $ TCGA.3C.AAAU.01A.11R.A41B.07: num  0 16.4 12.9 52.2 408.1 ...
##  $ TCGA.3C.AALI.01A.11R.A41B.07: num  0 9.27 17.38 69.76 563.89 ...
##  $ TCGA.3C.AALJ.01A.31R.A41B.07: num  0.907 11.623 9.229 154.297 1360.834 ...
##  $ TCGA.3C.AALK.01A.11R.A41B.07: num  0 12.1 11.1 143.9 865.5 ...
##  $ TCGA.4H.AAAK.01A.12R.A41B.07: num  0 6.85 14.43 84.21 766.38 ...
##  $ TCGA.5L.AAT0.01A.12R.A41B.07: num  0 3.99 13.61 114.26 807.74 ...
##  $ TCGA.5L.AAT1.01A.12R.A41B.07: num  0 0 10.6 116 1108.4 ...
##  $ TCGA.5T.A9QA.01A.11R.A41B.07: num  0 1.46 9 107.56 1420.5 ...
##  $ TCGA.A1.A0SB.01A.11R.A144.07: num  0 15.3 14.4 116.4 657.3 ...
##  $ TCGA.A1.A0SD.01A.11R.A115.07: num  0 9.52 11.32 60.26 977.92 ...
##   [list output truncated]
```

Create the transposed RNASeq data frame  (columns <-> rows)

```r
library(data.table)
df.RNAseq.t <- transpose(df.RNAseq) # transposed data frame
colnames(df.RNAseq.t)<-rownames(df.RNAseq)
rownames(df.RNAseq.t)<-colnames(df.RNAseq)

class(df.RNAseq.t) # Verify its class
```

```
## [1] "data.frame"
```

```r
str(df.RNAseq.t, list.len=10) # Check its structure
```

```
## 'data.frame':	1212 obs. of  20531 variables:
##  $ ?|100130426              : num  0 0 0.907 0 0 ...
##  $ ?|100133144              : num  16.36 9.27 11.62 12.09 6.85 ...
##  $ ?|100134869              : num  12.93 17.38 9.23 11.08 14.43 ...
##  $ ?|10357                  : num  52.2 69.8 154.3 143.9 84.2 ...
##  $ ?|10431                  : num  408 564 1361 866 766 ...
##  $ ?|136542                 : num  0 0 0 0 0 0 0 0 0 0 ...
##  $ ?|155060                 : num  1187 516 592 553 261 ...
##  $ ?|26823                  : num  0 1.087 0 0.414 0.425 ...
##  $ ?|280660                 : num  0 0.544 0 0 0 ...
##  $ ?|317712                 : num  0 0 0 0 0 0 0 0 0 0 ...
##   [list output truncated]
```

Import the data of the patients which ID matches with TNBC cases.

```r
df.TNBC.id <- read.csv("TNBC_merged_id.csv") # Data frame with TNBC patients' IDs
df.id <- gsub("-",".",df.TNBC.id$patient.bcr_patient_barcode_upper)
class(df.id) # Verify its class
```

```
## [1] "character"
```

```r
str(df.id) # Check its structure
```

```
##  chr [1:1097] "TCGA.5L.AAT0" "TCGA.5L.AAT1" "TCGA.A1.A0SP" "TCGA.A2.A04V" ...
```
Filter the values of the transposed data frame (df.RNASeq.t) to keep only the rows of information that match the IDs of the TNBC patients (df.TNBC.id)

```r
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:data.table':
## 
##     between, first, last
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
match.df <- paste("^(", paste0(df.id, collapse = "|"), ")", sep = "")

tnbc.df <- filter(df.RNAseq.t, Reduce("|", lapply(rownames(df.RNAseq.t), grepl, pattern = match.df)))
```

Check the structure of the resulting filtered data frame

```r
str(tnbc.df, list.len = 10) # Check its structure
```

```
## 'data.frame':	1212 obs. of  20531 variables:
##  $ ?|100130426              : num  0 0 0.907 0 0 ...
##  $ ?|100133144              : num  16.36 9.27 11.62 12.09 6.85 ...
##  $ ?|100134869              : num  12.93 17.38 9.23 11.08 14.43 ...
##  $ ?|10357                  : num  52.2 69.8 154.3 143.9 84.2 ...
##  $ ?|10431                  : num  408 564 1361 866 766 ...
##  $ ?|136542                 : num  0 0 0 0 0 0 0 0 0 0 ...
##  $ ?|155060                 : num  1187 516 592 553 261 ...
##  $ ?|26823                  : num  0 1.087 0 0.414 0.425 ...
##  $ ?|280660                 : num  0 0.544 0 0 0 ...
##  $ ?|317712                 : num  0 0 0 0 0 0 0 0 0 0 ...
##   [list output truncated]
```
## Normalizing the data

Applying log2 to all the data set

```r
tnbc.df.2 <- log2(tnbc.df + 1)

cbind(tnbc.df[40:100,2], tnbc.df.2[40:100,2]) # comparing the original values against the normalized ones
```

```
##          [,1]      [,2]
##  [1,]  1.4905 1.3164354
##  [2,]  0.0000 0.0000000
##  [3,]  7.3845 3.0677248
##  [4,]  3.5615 2.1895083
##  [5,]  6.0066 2.8087145
##  [6,]  1.0049 1.0035303
##  [7,]  9.3125 3.3663222
##  [8,]  7.1857 3.0331058
##  [9,]  2.5659 1.8342662
## [10,]  5.3773 2.6729457
## [11,]  9.8281 3.4367082
## [12,] 12.6197 3.7676230
## [13,]  7.4762 3.0834176
## [14,] 10.1456 3.4784024
## [15,]  4.0290 2.3302716
## [16,]  3.4684 2.1597583
## [17,]  2.1993 1.6777563
## [18,]  9.4756 3.3889610
## [19,]  7.7660 3.1319187
## [20,]  9.8895 3.4448658
## [21,]  4.3960 2.4318903
## [22,]  4.1760 2.3718376
## [23,]  8.0526 3.1783322
## [24,]  6.2239 2.8527779
## [25,]  0.6524 0.7245630
## [26,]  6.7911 2.9618270
## [27,]  2.6265 1.8585778
## [28,]  1.6315 1.3958854
## [29,]  5.9806 2.8033510
## [30,]  0.4387 0.5247658
## [31,] 12.3539 3.7391892
## [32,]  6.4181 2.8910497
## [33,]  8.6466 3.2700205
## [34,]  3.1723 2.0608429
## [35,] 20.3388 4.4154071
## [36,]  7.5193 3.0907349
## [37,]  4.9341 2.5690292
## [38,]  7.2508 3.0445340
## [39,]  4.5768 2.4794375
## [40,]  6.9474 2.9904830
## [41,]  1.5860 1.3707223
## [42,]  2.6025 1.8489984
## [43,] 23.8942 4.6377377
## [44,]  1.8358 1.5037558
## [45,]  8.5185 3.2507342
## [46,]  3.0155 2.0055796
## [47,] 17.4747 4.2074790
## [48,]  5.7893 2.7632628
## [49,] 12.1946 3.7218757
## [50,]  7.6688 3.1158323
## [51,] 28.8549 4.8998958
## [52,] 15.0630 4.0056695
## [53,] 14.5507 3.9589076
## [54,] 32.0154 5.0450672
## [55,]  2.2206 1.6873295
## [56,]  4.0091 2.3245514
## [57,] 11.3741 3.6292517
## [58,] 10.5037 3.5240261
## [59,]  3.9466 2.3064372
## [60,]  7.8319 3.1427238
## [61,] 13.5805 3.8659683
```


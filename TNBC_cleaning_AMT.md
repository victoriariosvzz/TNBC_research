---
title: "RNASeq"
author:
  - "Victoria Ríos"
  - "Antonio Martínez"
date: "febrero 08, 2021"
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

Importamos el data frame de la información con los id's de los pacientes que son Triple Negativo

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


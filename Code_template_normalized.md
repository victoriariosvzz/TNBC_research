---
title: "Code_template_normalized"
author: "Victoria Rios"
date: "16/3/2021"
output: html_document
---



#Análisis de datos de BRCA TN TCGA

Primero definimos la opción para no se lean los strings como factores.


```r
options(stringsAsFactors = F)
```

Primero leemos los datos clínicos.

```r
brca_clin <- read.delim(file = "BRCA.clin.merged.txt", sep = "\t", header = F,
row.names = 1, as.is = T)
```

Obtenemos los datos para los 3 receptores en los renglones 24, 29 y 1160.

```r
print(rownames(brca_clin)[c(24,29,1160)])
```

```
## [1] "patient.breast_carcinoma_estrogen_receptor_status"             
## [2] "patient.breast_carcinoma_progesterone_receptor_status"         
## [3] "patient.lab_proc_her2_neu_immunohistochemistry_receptor_status"
```

Índices de muestras triple negativo (TN)

```r
index_TN <- which(as.character(brca_clin[24,]) == "negative" &
as.character(brca_clin[29,]) == "negative" &
as.character(brca_clin[1160,]) == "negative")
```

Filtramos para tener datos sólo de los TN y obtenemos las claves de las muestras


```r
brca_clin <- brca_clin[,index_TN]
tn_samples <- toupper(as.character(brca_clin[21,]))
```

Leer datos de RNA-seq

```r
f <- "BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"
brca_rnaseq <- read.delim(file = f, sep="\t", as.is=T, row.names = 1)
brca_rnaseq <- data.matrix(brca_rnaseq[-1,])
```

Extraer gene symbols

```r
gene_symbols <- as.character(sapply(rownames(brca_rnaseq), function(x)
unlist(strsplit(x = x, split = "\\|"))[1]))
```

Eliminar renglones con gene symbol igual a “?”

```r
brca_rnaseq <- brca_rnaseq[which(gene_symbols != "?"),]
gene_symbols <- gene_symbols[which(gene_symbols != "?")]
```

Ver el tipo de muestra de acuerdo al TCGA barcode y quitar los normales (tipo “11”) https://docs.gdc.
cancer.gov/Encyclopedia/pages/TCGA_Barcode/

```r
sample_type <- substr(x = colnames(brca_rnaseq), start = 14, stop = 15)
brca_rnaseq <- brca_rnaseq[,-which(sample_type == "11")]
sample_type <- sample_type[-which(sample_type == "11")]
```

Hacer match entre las muestras TN y los de RNAseq

```r
tn_match <- match(gsub(pattern = "-", replacement = "\\.", x = tn_samples),
substr(x = colnames(brca_rnaseq), start = 1, stop = 12))
tn_match_na <- which(is.na(tn_match))
brca_clin <- brca_clin[,-tn_match_na]
tn_samples <- tn_samples[-tn_match_na]
tn_match <- tn_match[-tn_match_na]
brca_rnaseq <- brca_rnaseq[,tn_match]

brca_rnaseq_t <- t(brca_rnaseq) # Transpose of Rnaseq 
```

# Normalization of the dataset

```r
brca_rnaseq_n <- log2(brca_rnaseq_t + 1)
```


It's time to build the distance matrix. We use the euclidean distance method.

```r
temp <- brca_rnaseq_n
row.names(temp) <- seq(1, nrow(temp))

d <- dist(temp, method = 'euclidean')
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

![plot of chunk unnamed-chunk-146](figure/unnamed-chunk-146-1.png)

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

![plot of chunk unnamed-chunk-148](figure/unnamed-chunk-148-1.png)

Now we can see the three clusters enclosed in three different colored boxes. 

```r
suppressPackageStartupMessages(library(dendextend))
avg_dend_obj <- as.dendrogram(h)
avg_col_dend <- color_branches(h, h = 10)
plot(avg_col_dend)
```

![plot of chunk unnamed-chunk-149](figure/unnamed-chunk-149-1.png)

# Clustering the original data set (normalized)

```r
k <- kmeans(brca_rnaseq_n, centers = 3)
```

The cluster assignments are in the cluster component:

```r
groups <- k$cluster
```

Note that because the first center is chosen at random, the final clusters are random. We impose some stability by repeating the entire function several times and averaging the results. The number of random starting values to use can be assigned through the nstart argument.

```r
k <- kmeans(x = brca_rnaseq_n, centers = 3, nstart = 25, iter.max = 800)
```


## Heatmaps
A powerful visualization tool for discovering clusters or patterns in your data is the heatmap. The idea is simple: plot an image of your data matrix with colors used as the visual cue and both the columns and rows ordered according to the results of a clustering algorithm. We will demonstrate this with the tissue_gene_expression dataset. We will scale the rows of the gene expression matrix.

The first step is compute:

```r
brca_rnaseq_df <- as.data.frame(brca_rnaseq_n)

x <- sweep(brca_rnaseq_df, 2, colMeans(brca_rnaseq_df))
```

### Filtering features
If the information about clusters in included in just a few features, including all the features can add enough noise that detecting clusters becomes challenging. One simple approach to try to remove features with no information is to only include those with high variance. In the movie example, a user with low variance in their ratings is not really informative: all the movies seem about the same to them. Here is an example of how we can include only the features with high variance.

```r
library(matrixStats)
x <- as.matrix(x)
sds <- colSds(x, na.rm = TRUE)
o <- order(sds, decreasing = TRUE)[1:25]
heatmap(x[1:100,o], col = RColorBrewer::brewer.pal(11, "Spectral"))
```

![plot of chunk unnamed-chunk-154](figure/unnamed-chunk-154-1.png)


# Clustering PCA results (normalized)

## PCA

```r
# Principal components
x <- prcomp(brca_rnaseq_n, scale = FALSE)
```


```r
library(factoextra)

#Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
fviz_eig(x)
```

![plot of chunk unnamed-chunk-156](figure/unnamed-chunk-156-1.png)
## Graph of individuals. 
Individuals with a similar profile are grouped together.

```r
fviz_pca_ind(x,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = T,     # Avoid text overlapping,
             max.overlaps = Inf
             )
```

```
## Warning: ggrepel: 77 unlabeled data points (too many overlaps). Consider increasing max.overlaps
```

![plot of chunk unnamed-chunk-157](figure/unnamed-chunk-157-1.png)

## Graph of variables. 
Positive correlated variables point to the same side of the plot. Negative correlated variables point to opposite sides of the graph.

```r
fviz_pca_var(x,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
```

```
## Warning: ggrepel: 20498 unlabeled data points (too many overlaps). Consider increasing max.overlaps
```

![plot of chunk unnamed-chunk-158](figure/unnamed-chunk-158-1.png)

## Biplot of individuals and variables

```r
fviz_pca_biplot(x, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
                )
```

```
## Warning: ggrepel: 99 unlabeled data points (too many overlaps). Consider increasing max.overlaps
```

```
## Warning: ggrepel: 20498 unlabeled data points (too many overlaps). Consider increasing max.overlaps
```

![plot of chunk unnamed-chunk-159](figure/unnamed-chunk-159-1.png)
 

## K-means with PCA results

```r
# Coordinates of individuals
ind.coord <- as.data.frame(get_pca_ind(x)$coord)

head(ind.coord)
```

```
##                                  Dim.1     Dim.2      Dim.3    Dim.4      Dim.5    Dim.6        Dim.7
## TCGA.A1.A0SP.01A.11R.A084.07 -34.77630 12.851299  -19.74839 18.66707 -53.942113 21.15655  -33.2286425
## TCGA.AR.A1AR.01A.31R.A137.07 -12.48833 -5.735631   84.85662 39.09367  -9.645538 10.33063   -0.5714474
## TCGA.AR.A2LR.01A.12R.A18M.07 -71.44650 12.818727  -81.06502 87.20191 -28.064658 21.16727   -1.6607350
## TCGA.E2.A158.01A.11R.A12D.07  29.69381 65.406327  -36.24245 21.87599  84.784312  5.07634  -25.7076986
## TCGA.A1.A0SK.01A.12R.A084.07 -28.63871 25.416391 -119.88952 93.92394 -13.944131 11.98572 -122.6576637
## TCGA.A2.A04U.01A.11R.A115.07  80.07158 26.603180  -79.45512 17.69526  53.746263 27.90487  -19.2740806
##                                    Dim.8      Dim.9     Dim.10    Dim.11     Dim.12     Dim.13      Dim.14
## TCGA.A1.A0SP.01A.11R.A084.07  32.4255610 -26.015145  21.572180 -17.39858  13.857970 -10.217833  35.6784286
## TCGA.AR.A1AR.01A.31R.A137.07  39.0788690   8.555215 -33.890396 -36.53837 -43.722354 -35.173244   1.9950241
## TCGA.AR.A2LR.01A.12R.A18M.07   5.0074069 -24.549684   5.313364 -23.95346   4.399015  -8.152651  -9.0925930
## TCGA.E2.A158.01A.11R.A12D.07 -45.2112843   5.468948 -38.422196  20.64253   2.284226 -77.590106 -38.1949562
## TCGA.A1.A0SK.01A.12R.A084.07 -19.4385870 -33.928770  47.551929   3.47101  19.248899   9.476779   0.1937712
## TCGA.A2.A04U.01A.11R.A115.07   0.5597773  11.351354  28.290620  60.66396  32.250652 -71.003925  39.7496354
##                                  Dim.15      Dim.16     Dim.17     Dim.18    Dim.19     Dim.20    Dim.21
## TCGA.A1.A0SP.01A.11R.A084.07 -72.559014    8.777630 -31.298231  53.794384 -14.05044   1.710274 -6.790262
## TCGA.AR.A1AR.01A.31R.A137.07  34.390751  -10.801192  -3.441544 -40.386603 -16.22348 -34.981763 -5.767731
## TCGA.AR.A2LR.01A.12R.A18M.07  -3.237657  -35.713997   4.705444  14.778406 -13.80474 -19.643196 -5.460084
## TCGA.E2.A158.01A.11R.A12D.07 -14.813178  -39.698332 -12.650943  -5.350899  16.27310 -25.524496 -4.593593
## TCGA.A1.A0SK.01A.12R.A084.07 -41.310145    3.933072  29.874107  57.743373 -83.94849  13.043463 25.111826
## TCGA.A2.A04U.01A.11R.A115.07 -12.630417 -101.354521 -18.978948 -12.721911 -35.86549  91.113813 70.012119
##                                  Dim.22      Dim.23     Dim.24     Dim.25     Dim.26     Dim.27     Dim.28
## TCGA.A1.A0SP.01A.11R.A084.07  17.881731  -5.3165428 20.1068284  -0.691597   6.654926  -3.329901 -0.2039229
## TCGA.AR.A1AR.01A.31R.A137.07  -2.893656 -30.2995134 -9.7539986   3.688905   8.373834 -10.154687 27.6593718
## TCGA.AR.A2LR.01A.12R.A18M.07  -6.156003  10.8991984 -0.8064582  28.502472 -10.411215  47.717288 12.2349584
## TCGA.E2.A158.01A.11R.A12D.07 105.886788   3.1130133 -7.9966973 -83.912428 -57.481160  40.365055 42.4704993
## TCGA.A1.A0SK.01A.12R.A084.07  -1.430018 -52.2563652 -3.1198878  -6.348291  59.173615  45.625197 48.8838791
## TCGA.A2.A04U.01A.11R.A115.07 -59.054222  -0.1305914  1.7272367  94.722042 -79.962720 -74.407984 39.3697850
##                                  Dim.29     Dim.30     Dim.31     Dim.32     Dim.33    Dim.34     Dim.35
## TCGA.A1.A0SP.01A.11R.A084.07  -4.486138  25.886705  50.884530 -19.618041   5.869925  3.700963   3.282115
## TCGA.AR.A1AR.01A.31R.A137.07 -20.925819 -23.247879  29.940532  13.823289  20.075964 13.442860  11.189542
## TCGA.AR.A2LR.01A.12R.A18M.07 -28.345566   2.220038 -14.645862  -8.713508 -22.003731 -6.946353 -59.591977
## TCGA.E2.A158.01A.11R.A12D.07 -25.086178  13.039794  31.038508  84.522416   8.172003 20.723622 -34.603525
## TCGA.A1.A0SK.01A.12R.A084.07 -25.581714 -67.744024  -8.794554  18.105411 -45.290980 35.069676  -2.342779
## TCGA.A2.A04U.01A.11R.A115.07  11.211352  20.337096  10.352724 -12.419858  12.309184 39.552753  44.420504
##                                  Dim.36     Dim.37     Dim.38    Dim.39     Dim.40     Dim.41    Dim.42
## TCGA.A1.A0SP.01A.11R.A084.07   2.405007   9.549094 -26.711301 -16.97839 44.8076054 -29.623469 -14.26171
## TCGA.AR.A1AR.01A.31R.A137.07  17.295222  16.646672   4.963408  19.40930 26.4253581  -6.871067 -16.24350
## TCGA.AR.A2LR.01A.12R.A18M.07  -5.001119  19.341095  30.468511  10.22646 25.1251592  14.148309  41.70219
## TCGA.E2.A158.01A.11R.A12D.07 -32.088742  33.795619  31.561484 -39.09153  2.2453030 -28.676301 -67.40118
## TCGA.A1.A0SK.01A.12R.A084.07  31.297599   1.935832 -15.904347   7.18604 -5.4511911  11.404106 -18.99131
## TCGA.A2.A04U.01A.11R.A115.07 -17.614530 -14.705793  45.976012  32.04638 -0.1706915  -5.478884 -14.46473
##                                   Dim.43     Dim.44     Dim.45     Dim.46     Dim.47     Dim.48     Dim.49
## TCGA.A1.A0SP.01A.11R.A084.07 -13.3734927 -27.375004  46.618161  13.121367  14.589641  13.150492 -24.141448
## TCGA.AR.A1AR.01A.31R.A137.07   0.2104237  -2.331240   2.037285   3.110355  17.203560  -5.874369  -9.160856
## TCGA.AR.A2LR.01A.12R.A18M.07   1.4560774  -2.712796 -42.958334 -50.739931  31.732351 -45.115153 -14.999787
## TCGA.E2.A158.01A.11R.A12D.07   2.0446854  38.306722 -40.594512  27.365021 -32.805437  -8.864394  -3.350651
## TCGA.A1.A0SK.01A.12R.A084.07 -17.6431569  -4.612007  10.652018  12.079040  63.803726  10.394155  24.465588
## TCGA.A2.A04U.01A.11R.A115.07 -23.9610426  24.125021 -34.018182   1.007010  -9.409442  33.723508  33.418386
##                                  Dim.50     Dim.51     Dim.52     Dim.53     Dim.54      Dim.55      Dim.56
## TCGA.A1.A0SP.01A.11R.A084.07  -1.834895 -23.527313  22.960193  25.186690  19.404866  -3.5483306  -2.2734115
## TCGA.AR.A1AR.01A.31R.A137.07   3.897542  23.777579  30.602769  29.100707   3.754367  19.3084408  22.5895862
## TCGA.AR.A2LR.01A.12R.A18M.07  12.826717   1.838866  12.544074  -2.959855 -18.075030  -5.7647437 -36.2128413
## TCGA.E2.A158.01A.11R.A12D.07 -30.996763 -41.299983 -32.829775 -39.014851   6.574026   0.1069002   2.9544631
## TCGA.A1.A0SK.01A.12R.A084.07  10.530565 -36.013504   3.578827  -9.964229 -37.885429  32.0526889 -17.3864777
## TCGA.A2.A04U.01A.11R.A115.07  16.267969  13.658379   9.305501  12.183873  -4.411118 -24.8341588   0.4494479
##                                   Dim.57      Dim.58     Dim.59     Dim.60     Dim.61     Dim.62     Dim.63
## TCGA.A1.A0SP.01A.11R.A084.07 -20.6815533  -3.3430677  24.405171   7.857186   4.906818  23.767414  15.969967
## TCGA.AR.A1AR.01A.31R.A137.07  -9.5212468 -44.7757193  41.017295 -14.421466  -8.562837 -19.494718 -12.093253
## TCGA.AR.A2LR.01A.12R.A18M.07  -0.5129237   2.6447380 -10.555549   1.595838 -52.347472  69.828737 -75.346388
## TCGA.E2.A158.01A.11R.A12D.07  13.3546505   2.0872150  21.773961 -16.124500 -10.124329  28.470956   2.393812
## TCGA.A1.A0SK.01A.12R.A084.07  36.5630223 -14.4618609 -31.519608 -13.968547  33.143522 -33.694567  54.707240
## TCGA.A2.A04U.01A.11R.A115.07  -1.7475333  -0.5462722  -5.230344  -5.857892  -4.988451   5.746244   4.250811
##                                  Dim.64     Dim.65     Dim.66     Dim.67     Dim.68     Dim.69     Dim.70
## TCGA.A1.A0SP.01A.11R.A084.07 -51.252036 -31.589024  -4.453504 -20.170967  20.100984  14.801319  36.965025
## TCGA.AR.A1AR.01A.31R.A137.07  40.512654  -7.128026 -10.668931   2.558195 -16.931709   4.493070   3.995374
## TCGA.AR.A2LR.01A.12R.A18M.07   7.564820  33.364674   1.315507  15.635552  56.541599 -10.644127  -7.915220
## TCGA.E2.A158.01A.11R.A12D.07   3.881929 -11.961267 -16.563567 -13.295722  -2.008618  16.345333 -12.043054
## TCGA.A1.A0SK.01A.12R.A084.07  -5.568843 -15.837847 -20.561271   4.065798 -25.772274  -5.189910 -36.046294
## TCGA.A2.A04U.01A.11R.A115.07  14.841924   2.193495 -32.649997  -7.339086  12.786632   4.413811  -8.193602
##                                  Dim.71     Dim.72     Dim.73    Dim.74     Dim.75     Dim.76     Dim.77
## TCGA.A1.A0SP.01A.11R.A084.07   7.039967  57.652925 -68.512775 -8.376189 -15.168022 -43.215322  12.425105
## TCGA.AR.A1AR.01A.31R.A137.07   5.588818  -8.002833 -30.939828 13.897478 -32.966748  14.370289  46.104140
## TCGA.AR.A2LR.01A.12R.A18M.07 -73.101640  37.216637  10.280272 30.807354  -1.685891  26.080439  -2.610625
## TCGA.E2.A158.01A.11R.A12D.07 -17.530399  -6.451361  -3.248381  8.169418   8.539531   0.853271 -17.354199
## TCGA.A1.A0SK.01A.12R.A084.07 -10.361159 -12.409765  -5.021559  2.478735 -14.424965  -3.349382 -16.162053
## TCGA.A2.A04U.01A.11R.A115.07  17.837718  -9.339983  -9.984512 -5.299915  -8.530021  -1.660178  23.290832
##                                 Dim.78     Dim.79     Dim.80    Dim.81     Dim.82     Dim.83     Dim.84
## TCGA.A1.A0SP.01A.11R.A084.07  19.71306  35.920693   9.472217 10.740176 -29.158516 -15.488296 -14.064879
## TCGA.AR.A1AR.01A.31R.A137.07  13.16975 -39.346785 -26.596783 32.792911  14.902895  -4.442085   4.569741
## TCGA.AR.A2LR.01A.12R.A18M.07  25.01312  55.132164   8.187548 -5.951364  24.402877   5.996857 -25.138483
## TCGA.E2.A158.01A.11R.A12D.07 -14.84326 -10.432895  -5.168031  5.828410  -5.674991 -12.382179   4.573035
## TCGA.A1.A0SK.01A.12R.A084.07  11.38127  -4.852673  -6.098550 -6.487277  10.887323   9.419676  20.879947
## TCGA.A2.A04U.01A.11R.A115.07  17.87854   6.360144  -3.592488 -6.896955  -2.293776   2.242428 -18.184612
##                                  Dim.85     Dim.86     Dim.87     Dim.88      Dim.89     Dim.90     Dim.91
## TCGA.A1.A0SP.01A.11R.A084.07 -19.712177   5.702813  34.851696   2.403405 -50.1620834 -38.740456 -21.366452
## TCGA.AR.A1AR.01A.31R.A137.07 -10.993212   1.939056  -1.363839  14.072878  -0.7266219  20.028133  -5.691862
## TCGA.AR.A2LR.01A.12R.A18M.07 -14.038698  22.892472  -6.057535   4.984753  -8.9876572  17.543872 -32.058130
## TCGA.E2.A158.01A.11R.A12D.07   3.556850 -12.592846  22.742886  22.189960   0.1405586 -12.740312  11.807307
## TCGA.A1.A0SK.01A.12R.A084.07  -7.434808  -8.750010 -30.011077 -11.013506   7.3997392   9.405639  26.474008
## TCGA.A2.A04U.01A.11R.A115.07   4.050987  -4.285420  15.662490  -2.515705  18.3487369   4.985036 -13.293763
##                                  Dim.92     Dim.93     Dim.94     Dim.95      Dim.96     Dim.97      Dim.98
## TCGA.A1.A0SP.01A.11R.A084.07  8.6563918  23.362425 31.0633257 -20.944230  -0.7617463  -6.308798 -33.8253192
## TCGA.AR.A1AR.01A.31R.A137.07  3.9321331  27.406989 15.7774341  19.855809 -55.3975045 -24.048802 -23.7001138
## TCGA.AR.A2LR.01A.12R.A18M.07 -0.7337504   1.687396 -0.3097997  29.328828 -17.3431687   1.566800  23.0450325
## TCGA.E2.A158.01A.11R.A12D.07 15.5272774   5.135998  3.6296475  11.803053   9.7146222 -11.071121 -14.9298832
## TCGA.A1.A0SK.01A.12R.A084.07 -9.2665013 -12.207698 -4.9137327   6.342720   0.8435729  -6.337024  -0.4530305
## TCGA.A2.A04U.01A.11R.A115.07  8.7508190   4.889520 -7.8994459  -3.437078   5.5568058   5.460420  -8.8599582
##                                 Dim.99    Dim.100      Dim.101    Dim.102    Dim.103    Dim.104    Dim.105
## TCGA.A1.A0SP.01A.11R.A084.07 55.475294   1.150328  22.77722442 -19.004764 -15.771729  25.602762  18.720141
## TCGA.AR.A1AR.01A.31R.A137.07 18.542254  -5.916490 -63.55146757 -29.681254 -24.107308 -32.091961  36.779852
## TCGA.AR.A2LR.01A.12R.A18M.07  8.080402 -16.796554   5.78003204 -11.054034  -1.637420  -0.144525  -7.494172
## TCGA.E2.A158.01A.11R.A12D.07  1.141194   3.389336  -8.64375824   8.091188  -1.398370  13.510710 -10.432846
## TCGA.A1.A0SK.01A.12R.A084.07 -9.414820  -2.493263 -13.14975791   6.490876  -5.353528  -7.351628   5.990928
## TCGA.A2.A04U.01A.11R.A115.07  6.326589   5.256607  -0.09513982   4.046043   2.074080  -2.177230 -11.569351
##                                 Dim.106    Dim.107   Dim.108      Dim.109    Dim.110    Dim.111    Dim.112
## TCGA.A1.A0SP.01A.11R.A084.07 -28.032201 -9.6884617 15.235802  18.12287768  -1.998397   2.862327  -6.936057
## TCGA.AR.A1AR.01A.31R.A137.07  43.384914 48.2772238 45.007975 -47.23173279 -34.880838 -13.563958  -4.942301
## TCGA.AR.A2LR.01A.12R.A18M.07   7.013450  1.0452628  2.297760  -1.06507639 -20.860566  -6.846620 -10.547611
## TCGA.E2.A158.01A.11R.A12D.07   3.730789 -6.7819445  5.261060   6.22709955  -6.050519   6.907259   1.577559
## TCGA.A1.A0SK.01A.12R.A084.07  -7.904604 -0.2339443  4.347976  -3.54975054  -1.029316  -6.971576   2.868685
## TCGA.A2.A04U.01A.11R.A115.07 -11.561800  4.8768069  5.103851   0.03944796  -2.453370   3.885556  -4.780379
##                                 Dim.113     Dim.114       Dim.115
## TCGA.A1.A0SP.01A.11R.A084.07  15.727734  -8.3014872  1.392608e-13
## TCGA.AR.A1AR.01A.31R.A137.07   2.313104  18.1036912  1.223915e-13
## TCGA.AR.A2LR.01A.12R.A18M.07 -12.186824   7.4571921 -6.368575e-14
## TCGA.E2.A158.01A.11R.A12D.07   4.015586 -11.7068084  7.283547e-14
## TCGA.A1.A0SK.01A.12R.A084.07  -8.309488  -0.5300006  6.061478e-14
## TCGA.A2.A04U.01A.11R.A115.07   4.411809   3.4555030  4.064160e-14
```



```r
kc <- kmeans(ind.coord, 3)

plot(ind.coord[,1:2],col=factor(kc$cluster))
```

![plot of chunk unnamed-chunk-161](figure/unnamed-chunk-161-1.png)

```r
# Add clusters obtained using the K-means algorithm
ind.coord$cluster <- factor(k$cluster)
```


```r
# Percentage of variance explained by dimensions
eigenvalue <- round(get_eigenvalue(x), 1)
variance.percent <- eigenvalue$variance.percent

head(eigenvalue)
```

```
##       eigenvalue variance.percent cumulative.variance.percent
## Dim.1     3760.7              4.5                         4.5
## Dim.2     2447.8              3.0                         7.5
## Dim.3     2124.6              2.6                        10.1
## Dim.4     1999.4              2.4                        12.5
## Dim.5     1603.8              1.9                        14.4
## Dim.6     1571.8              1.9                        16.3
```



```r
library(ggpubr)

ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex", size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )" )
) +
  stat_mean(aes(color = cluster), size = 4)
```

![plot of chunk unnamed-chunk-164](figure/unnamed-chunk-164-1.png)
# K-means with PCA results (not normalized)

## Clustering PCA results

## PCA

```r
# Principal components
x <- prcomp(brca_rnaseq_t, scale = FALSE)
```


```r
#Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
fviz_eig(x)
```

![plot of chunk unnamed-chunk-166](figure/unnamed-chunk-166-1.png)
## Graph of individuals. 
Individuals with a similar profile are grouped together.

```r
fviz_pca_ind(x,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = T,     # Avoid text overlapping,
             max.overlaps = Inf
             )
```

```
## Warning: ggrepel: 90 unlabeled data points (too many overlaps). Consider increasing max.overlaps
```

![plot of chunk unnamed-chunk-167](figure/unnamed-chunk-167-1.png)

## Graph of variables. 
Positive correlated variables point to the same side of the plot. Negative correlated variables point to opposite sides of the graph.

```r
fviz_pca_var(x,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
```

```
## Warning: ggrepel: 20502 unlabeled data points (too many overlaps). Consider increasing max.overlaps
```

![plot of chunk unnamed-chunk-168](figure/unnamed-chunk-168-1.png)

## Biplot of individuals and variables

```r
fviz_pca_biplot(x, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
                )
```

```
## Warning: ggrepel: 115 unlabeled data points (too many overlaps). Consider increasing max.overlaps
```

```
## Warning: ggrepel: 20502 unlabeled data points (too many overlaps). Consider increasing max.overlaps
```

![plot of chunk unnamed-chunk-169](figure/unnamed-chunk-169-1.png)
 

## Getting the coordinates of the PCA results as df

```r
# Coordinates of individuals
ind.coord <- as.data.frame(get_pca_ind(x)$coord)

head(ind.coord)
```

```
##                                   Dim.1     Dim.2      Dim.3        Dim.4     Dim.5      Dim.6      Dim.7
## TCGA.A1.A0SP.01A.11R.A084.07   8243.172 -15969.63  33717.902    -402.6725  12981.35 -27518.472  21986.561
## TCGA.AR.A1AR.01A.31R.A137.07  91012.882  83816.78   6000.662  -25101.9226 -73520.62  27929.323 -82940.885
## TCGA.AR.A2LR.01A.12R.A18M.07 -13156.537  28459.46 -16775.945  -29660.9140  74438.20 -20986.197   4709.628
## TCGA.E2.A158.01A.11R.A12D.07   7182.063  27692.11  20756.545  -45730.8206  37216.32  18404.675 -33246.020
## TCGA.A1.A0SK.01A.12R.A084.07  28343.501  28540.03  43244.433 -111898.9067  62659.17 -34468.081  17901.025
## TCGA.A2.A04U.01A.11R.A115.07 -96842.508  20495.49  13498.949  -24198.7547  18928.30   1101.741 -13315.074
##                                   Dim.8      Dim.9      Dim.10     Dim.11    Dim.12     Dim.13     Dim.14
## TCGA.A1.A0SP.01A.11R.A084.07 -26008.330  19288.453 -14395.0856  9291.5593 -60466.18  58903.423 -22652.801
## TCGA.AR.A1AR.01A.31R.A137.07 -45238.068  22709.930   1328.6710  8112.9212  73555.93 -14027.892 -62834.591
## TCGA.AR.A2LR.01A.12R.A18M.07  -9421.359  14902.018   8693.4560  1292.2461 -10256.00  11304.119   6740.231
## TCGA.E2.A158.01A.11R.A12D.07  33987.214 -16358.724    846.1996   523.6296  16201.62   3814.474  42372.316
## TCGA.A1.A0SK.01A.12R.A084.07 -15342.980  91048.538  20052.0784 13848.4274 -48701.07  18128.141  78777.386
## TCGA.A2.A04U.01A.11R.A115.07  10639.277   3458.297    275.2243  7097.8033  14030.33   8518.960  13460.006
##                                  Dim.15     Dim.16     Dim.17     Dim.18     Dim.19     Dim.20     Dim.21
## TCGA.A1.A0SP.01A.11R.A084.07  9585.7341 -10350.130  14749.016 -25898.049  21007.350 -16946.426 -31328.767
## TCGA.AR.A1AR.01A.31R.A137.07 -2246.8968  18238.348  17793.229  31162.723  24823.845   9280.965 -15980.845
## TCGA.AR.A2LR.01A.12R.A18M.07 16937.6764   7703.133  24631.474   2286.152  31608.857  -6096.921  -7807.341
## TCGA.E2.A158.01A.11R.A12D.07 19398.5981   6948.197  -5049.507  20531.016 -61339.537 -21743.007 -33864.064
## TCGA.A1.A0SK.01A.12R.A084.07 95880.8898  38794.113 -27883.000 -11059.916  43592.310 -23473.400  40473.612
## TCGA.A2.A04U.01A.11R.A115.07  -321.0825  -3219.611 -10382.134  -1619.060  -2377.485  -4095.632  -4772.906
##                                  Dim.22      Dim.23      Dim.24    Dim.25     Dim.26     Dim.27      Dim.28
## TCGA.A1.A0SP.01A.11R.A084.07  41960.270   -212.2412  -19907.718 -11697.07  20699.405  31856.266 -12552.3179
## TCGA.AR.A1AR.01A.31R.A137.07  23044.969 -84745.9262 -119571.718 -11460.38 -69242.803  18222.449 -12257.5036
## TCGA.AR.A2LR.01A.12R.A18M.07  10008.620   7954.0137    2492.191 -13273.42 -14770.074 -15112.108   7197.4997
## TCGA.E2.A158.01A.11R.A12D.07 -13902.776 -66587.1998   13987.494 -69164.64 -14700.907 -17333.096  15724.9721
## TCGA.A1.A0SK.01A.12R.A084.07   2624.729  35788.2725  -63114.075 -60612.62  -6465.644  74518.223 -27187.9414
## TCGA.A2.A04U.01A.11R.A115.07  -4342.375    349.4517    3958.493  -5353.45   5869.246  -3934.622   -597.2816
##                                  Dim.29    Dim.30     Dim.31      Dim.32     Dim.33     Dim.34     Dim.35
## TCGA.A1.A0SP.01A.11R.A084.07  -6450.234  12622.83  16185.986 13030.59450 -11704.318  -2420.545 -26338.697
## TCGA.AR.A1AR.01A.31R.A137.07 -70882.819 -39051.33  19901.216 60078.39264  54867.383  10216.399 -63101.244
## TCGA.AR.A2LR.01A.12R.A18M.07  -4609.517  16829.86  -8926.392   -59.53103  -3610.100  -2435.830  -4643.678
## TCGA.E2.A158.01A.11R.A12D.07 -57231.075 -26260.01  56557.306  4166.93024 -65612.125  -2157.567  77477.855
## TCGA.A1.A0SK.01A.12R.A084.07  12071.578  16476.22 -37657.119  6205.58897 -19269.443 -42487.909  -5401.662
## TCGA.A2.A04U.01A.11R.A115.07  -7209.950 -11909.70   1649.229  2758.98170  -9542.712  -8657.017  -5409.698
##                                  Dim.36    Dim.37    Dim.38    Dim.39      Dim.40     Dim.41     Dim.42
## TCGA.A1.A0SP.01A.11R.A084.07 -14817.080  17923.55 21400.201  5918.192 -41321.8859  12055.288  30110.070
## TCGA.AR.A1AR.01A.31R.A137.07 -14909.067 -24590.21  7095.390 15032.655 -11865.1908   3873.684 -27515.790
## TCGA.AR.A2LR.01A.12R.A18M.07  16541.681 -18519.05 -6426.847 12880.734  11022.2998 -24223.926  20631.415
## TCGA.E2.A158.01A.11R.A12D.07 -50993.075  53160.65 -2402.772 -5189.771 -54187.3475  -7952.787  80212.569
## TCGA.A1.A0SK.01A.12R.A084.07  41154.226  45261.10 -3206.423 -3434.718 -22403.0312 -20901.987 -35532.948
## TCGA.A2.A04U.01A.11R.A115.07   4235.311  13207.92 -3838.296 -1981.392    396.1441   1549.460   9199.472
##                                  Dim.43     Dim.44     Dim.45     Dim.46      Dim.47     Dim.48      Dim.49
## TCGA.A1.A0SP.01A.11R.A084.07   7494.938   5835.882  41809.675 -24264.650  -1554.3913 -18111.040 -37188.3496
## TCGA.AR.A1AR.01A.31R.A137.07  12062.412  26715.361  -4987.337  24787.324 -13823.9675   6132.200  17665.8275
## TCGA.AR.A2LR.01A.12R.A18M.07  -6617.034  18101.837  -3231.145 -22276.723 -21613.0382   8499.319 -28708.9426
## TCGA.E2.A158.01A.11R.A12D.07 -21511.696 -54177.042  20430.597  30834.859  37563.4146  51287.092   3150.2570
## TCGA.A1.A0SK.01A.12R.A084.07  40543.405 -20463.886 -26618.249  40989.240  -1128.4172 -19888.504  27912.1212
## TCGA.A2.A04U.01A.11R.A115.07  -8765.489   4295.924 -24881.918   2968.838   -872.0069  -9114.042    424.7231
##                                  Dim.50     Dim.51     Dim.52     Dim.53     Dim.54     Dim.55     Dim.56
## TCGA.A1.A0SP.01A.11R.A084.07  17832.607  10990.349  16454.508   3821.439  20820.992 -35136.257 -31424.010
## TCGA.AR.A1AR.01A.31R.A137.07  -1346.762  -7625.104 -18749.643  -1283.205   9501.687 -22389.328  33318.566
## TCGA.AR.A2LR.01A.12R.A18M.07  28283.069 -25399.025   3510.859  -9491.287 -13591.448 -23368.498  -4419.298
## TCGA.E2.A158.01A.11R.A12D.07   3706.100 -44583.185  18435.360 -14524.837  21347.082 -38259.984  12001.712
## TCGA.A1.A0SK.01A.12R.A084.07 -20338.200  25664.014 -25859.851   8733.758  14342.438  12465.014   2355.039
## TCGA.A2.A04U.01A.11R.A115.07  -4453.003  15232.737  -3329.077  -8449.795   9880.429   3402.538   5372.977
##                                   Dim.57     Dim.58      Dim.59     Dim.60     Dim.61     Dim.62     Dim.63
## TCGA.A1.A0SP.01A.11R.A084.07 -46387.0618  30577.314  35266.4027  30492.752 -47566.802 -32339.977 -23906.042
## TCGA.AR.A1AR.01A.31R.A137.07 -16288.8491  13761.839    280.4328  13687.506  30333.460  11111.492  -4059.153
## TCGA.AR.A2LR.01A.12R.A18M.07    686.0594  53415.195   7750.5376  -9454.612  -6540.784  21815.268 -19930.481
## TCGA.E2.A158.01A.11R.A12D.07  14558.9351  25651.937 -25005.1517  23297.612 -23927.789 -15306.915 -13678.218
## TCGA.A1.A0SK.01A.12R.A084.07   7638.6245 -18313.638    465.3323  23183.205  18017.736 -12566.262   8496.608
## TCGA.A2.A04U.01A.11R.A115.07  11566.0060   1439.242  10316.1969 -22879.587  17615.166  -6995.882 -37276.661
##                                  Dim.64      Dim.65     Dim.66     Dim.67      Dim.68     Dim.69     Dim.70
## TCGA.A1.A0SP.01A.11R.A084.07 -30689.123 -84261.0641  34251.000  43119.198 -24565.5811  -4371.419 -15955.756
## TCGA.AR.A1AR.01A.31R.A137.07  -9981.906  -9891.3456 -16167.295 -25243.868  -6197.6214 -11667.889   5098.419
## TCGA.AR.A2LR.01A.12R.A18M.07  16174.093  16073.1703   2162.015 -22685.413 -20014.8458 -56480.973  30436.990
## TCGA.E2.A158.01A.11R.A12D.07   2050.899  30429.8420 -28085.076  -5803.351  12974.5344  -3651.377 -26395.017
## TCGA.A1.A0SK.01A.12R.A084.07   4230.510  28601.4880 -18419.097 -19795.778 -14703.5788  17738.797    188.319
## TCGA.A2.A04U.01A.11R.A115.07 -21685.052   -429.5298  29282.635  17241.253    501.3464  -8972.609   8023.526
##                                  Dim.71     Dim.72     Dim.73     Dim.74     Dim.75    Dim.76     Dim.77
## TCGA.A1.A0SP.01A.11R.A084.07 -26987.386  -3863.954  29104.685   5767.681  16238.919  8319.971  -3297.700
## TCGA.AR.A1AR.01A.31R.A137.07  -2998.682  13382.280  -5240.716   2951.226  -1833.702 -1024.271   8940.078
## TCGA.AR.A2LR.01A.12R.A18M.07  36960.584 -65856.382 -33430.572 -34948.484 -35510.221 33231.973 -32181.914
## TCGA.E2.A158.01A.11R.A12D.07  18062.240  18109.985   2163.178   3294.067  -1631.649  3503.519  11854.773
## TCGA.A1.A0SK.01A.12R.A084.07  11390.441  -4718.846  -8430.031  -4383.649   8150.779 -7124.410   8006.371
## TCGA.A2.A04U.01A.11R.A115.07   5088.762  13297.746  28780.922  54549.153  -8684.933 -1914.905 -16320.463
##                                  Dim.78      Dim.79     Dim.80     Dim.81      Dim.82     Dim.83     Dim.84
## TCGA.A1.A0SP.01A.11R.A084.07  -4648.332 -111593.434  32419.663   9206.510    439.5604  21696.908 -4870.5311
## TCGA.AR.A1AR.01A.31R.A137.07   2666.909  -15140.401   8657.713   2793.587   6730.4435  -7133.213  2588.6789
## TCGA.AR.A2LR.01A.12R.A18M.07 -26166.537   10708.827  93803.626 -66205.969 -18364.8399   8448.097 14261.4357
## TCGA.E2.A158.01A.11R.A12D.07 -11866.264   13572.947  16951.459   5619.053   5109.4510  13420.388 -4288.1043
## TCGA.A1.A0SK.01A.12R.A084.07   8953.903    9572.416 -19975.062  10876.831   2760.1149   5305.364  -370.1069
## TCGA.A2.A04U.01A.11R.A115.07  22735.239  -23018.112 -12604.118  15264.070  -6018.8075 -15852.028 65809.8291
##                                   Dim.85     Dim.86       Dim.87     Dim.88       Dim.89      Dim.90
## TCGA.A1.A0SP.01A.11R.A084.07 -5076.39382 -31794.192   2467.56802   15796.86  13219.21918 -17492.1513
## TCGA.AR.A1AR.01A.31R.A137.07   -95.12124   2657.741   3873.75184   -1647.88  -6363.00759    270.0358
## TCGA.AR.A2LR.01A.12R.A18M.07  6337.23015  27077.172  29153.56067  -29627.67 -40839.60533 -19677.9784
## TCGA.E2.A158.01A.11R.A12D.07   864.26899  -2636.719 -17653.11663   15404.75   4449.10520   1007.6573
## TCGA.A1.A0SK.01A.12R.A084.07 10454.78213 -13935.512    -16.20552   15279.47   5558.83043 -10907.9100
## TCGA.A2.A04U.01A.11R.A115.07  8797.05044  74726.137 -32983.61634 -112230.19     77.91193 -95328.7994
##                                  Dim.91     Dim.92    Dim.93     Dim.94     Dim.95      Dim.96     Dim.97
## TCGA.A1.A0SP.01A.11R.A084.07  -6166.300 -10577.769  6581.899   8298.433 -14546.790  -7657.7231   9673.936
## TCGA.AR.A1AR.01A.31R.A137.07  -1114.205   8595.836  4667.719  -2690.472  -2992.440   -296.6319  -8770.410
## TCGA.AR.A2LR.01A.12R.A18M.07 -25546.455  50611.067 -2761.288 -28009.783  34994.978 -27758.4585  16288.691
## TCGA.E2.A158.01A.11R.A12D.07  10373.814   4829.874  6216.838   2122.443 -17216.601    496.2396   3165.126
## TCGA.A1.A0SK.01A.12R.A084.07  -5240.346   9550.446  7125.998  -5751.493  -6585.141  -2312.4466  -1366.415
## TCGA.A2.A04U.01A.11R.A115.07  53359.114  15558.397 -2807.333  24919.804   6383.277  36121.6009 -29747.686
##                                   Dim.98    Dim.99    Dim.100    Dim.101     Dim.102  Dim.103    Dim.104
## TCGA.A1.A0SP.01A.11R.A084.07   2802.9094  4670.257 -14590.345 -10878.481   2642.2256 7667.250  -5718.151
## TCGA.AR.A1AR.01A.31R.A137.07   7221.1930  1930.103  -8734.230   1404.901   1349.5253 1036.978  -4359.356
## TCGA.AR.A2LR.01A.12R.A18M.07   -478.3200 14360.378  18244.123   6508.452  17571.4041 1134.662   5481.299
## TCGA.E2.A158.01A.11R.A12D.07   1298.7540 -1678.524   4951.051  -7835.457   -441.3185 3825.917 -10949.119
## TCGA.A1.A0SK.01A.12R.A084.07    563.1625 -9278.129  -3378.593  -8003.718   4868.1567 1565.925  -1675.445
## TCGA.A2.A04U.01A.11R.A115.07 -31608.5767 -7466.445 -14411.676 -10281.240 -12192.8556 1157.241  12952.560
##                                 Dim.105   Dim.106   Dim.107   Dim.108   Dim.109    Dim.110    Dim.111
## TCGA.A1.A0SP.01A.11R.A084.07  -171.5965  8932.300 -3788.093  3840.577 -1755.323  2259.0378 -4003.1128
## TCGA.AR.A1AR.01A.31R.A137.07 -3408.3588 -1491.566  4503.434 -1820.119  -972.269 -1147.9286  3684.9508
## TCGA.AR.A2LR.01A.12R.A18M.07 -3673.6029  4143.622 11250.168 -8691.531 -4073.611 10760.2630  2333.9381
## TCGA.E2.A158.01A.11R.A12D.07   296.4453 -1992.548 -4241.317 -2028.937  5056.914   533.0822  -190.3658
## TCGA.A1.A0SK.01A.12R.A084.07  1686.0503  1917.306  3737.778 -3139.718 -4387.189  1077.5490   817.8213
## TCGA.A2.A04U.01A.11R.A115.07 -8250.8067 23631.373 -7107.954 10550.252 12607.482 -6685.6494 16033.6751
##                                 Dim.112    Dim.113   Dim.114       Dim.115
## TCGA.A1.A0SP.01A.11R.A084.07  2603.7074  -1324.357 1512.1789 -1.056761e-10
## TCGA.AR.A1AR.01A.31R.A137.07  3295.9094   1289.173 2355.6513 -1.584514e-10
## TCGA.AR.A2LR.01A.12R.A18M.07   179.8578 -11731.219 1224.5489 -1.076907e-10
## TCGA.E2.A158.01A.11R.A12D.07  4250.2081  -2425.956  923.0512 -4.688130e-11
## TCGA.A1.A0SK.01A.12R.A084.07 -1235.8662  -3412.250  823.8293  2.865349e-10
## TCGA.A2.A04U.01A.11R.A115.07   588.6379   2019.171  520.6298 -6.990211e-11
```

## K means

```r
kc <- kmeans(ind.coord, 3)

plot(ind.coord[,1:2],col=factor(kc$cluster))
```

![plot of chunk unnamed-chunk-171](figure/unnamed-chunk-171-1.png)


```r
# Add clusters obtained using the K-means algorithm
ind.coord$cluster <- factor(k$cluster)
```


```r
# Percentage of variance explained by dimensions
eigenvalue <- round(get_eigenvalue(x), 1)
variance.percent <- eigenvalue$variance.percent

head(eigenvalue)
```

```
##       eigenvalue variance.percent cumulative.variance.percent
## Dim.1 4081548676              5.3                         5.3
## Dim.2 1754939455              2.3                         7.6
## Dim.3 1589841941              2.1                         9.7
## Dim.4 1231535632              1.6                        11.3
## Dim.5 1163653345              1.5                        12.8
## Dim.6 1025456073              1.3                        14.2
```



```r
library(ggpubr)

ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex", size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )" )
) +
  stat_mean(aes(color = cluster), size = 4)
```

![plot of chunk unnamed-chunk-174](figure/unnamed-chunk-174-1.png)



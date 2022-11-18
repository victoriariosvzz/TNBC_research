#### Clustering pipeline functions ####


# Function to prepare the clinical data to be consistent and homogeneous
feature_engineering_clinical_custom = function(clinical_dataframe, features_to_engineer) {
  # Define a loop to go over each feature and make the classes compliant across
  # data source denomination
  
  for (feature in features_to_engineer) {
    # Do a specific transformation depending on the clinical parameter
    if (feature == "histological_type") {
      clinical_dataframe$histological_type <-
        mapvalues(
          clinical_dataframe$histological_type,
          c("Mixed (not specified)", "Not specified"),
          c("Mixed", NA)
        )
    } else if (feature == "histological_grade") {
      clinical_dataframe$histological_grade <-
        mapvalues(
          clinical_dataframe$histological_grade,
          c(
            "UNKNOWN",
            "N/A",
            "Unknown",
            "II-III",
            "I LOW GRADE (WELL DIFFERENTIATED)",
            "3",
            "2",
            "1"
          ),
          c(NA, NA, NA, 2, 1, 3, 2, 1)
        )
      
      # Change to factor
      clinical_dataframe$histological_grade <-
        as.factor(clinical_dataframe$histological_grade)
      
    } else if (feature == "chemotherapy") {
      clinical_dataframe$chemotherapy <-
        mapvalues(
          clinical_dataframe$chemotherapy,
          c(
            "Doxorubicin",
            "Paclitaxel",
            "Adriamicin",
            "N/A",
            "Epirubicin",
            "NO",
            "YES"
          ),
          c(1, 1, 1, NA, 1, 0, 1)
        )
      
      # Change to factor
      clinical_dataframe$chemotherapy <-
        as.factor(clinical_dataframe$chemotherapy)
      
    } else if (feature == "histological_grade") {
      clinical_dataframe$histological_grade <-
        mapvalues(clinical_dataframe$histological_grade,
                  c("1",
                    "2",
                    "3"),
                  c(1, 2, 3))
      
      # Change to factor
      clinical_dataframe$histological_grade <-
        as.factor(clinical_dataframe$histological_grade)
      
    } else if (feature == "tabaquism") {
      clinical_dataframe$tabaquism <-
        mapvalues(clinical_dataframe$tabaquism,
                  c("Neg",
                    "Pos"),
                  c(0, 1))
      
      # Change to factor
      clinical_dataframe$tabaquism <- as.factor(clinical_dataframe$tabaquism)
      
    } else if (feature == "hispanic") {
      clinical_dataframe$hispanic <-
        mapvalues(clinical_dataframe$hispanic,
                  c("NO",
                    "YES"),
                  c(0, 1))
      
      # Change to factor
      clinical_dataframe$hispanic <- as.factor(clinical_dataframe$hispanic)
      
    } else if (feature == "Relapse") {
      clinical_dataframe$Relapse <-
        mapvalues(clinical_dataframe$Relapse,
                  c("TRUE",
                    "FALSE"),
                  c(1, 0))
      
      # Change to factor
      clinical_dataframe$Relapse <- as.factor(clinical_dataframe$Relapse)
      
    } else if (feature == "TNBC_status") {
      clinical_dataframe$TNBC_status <-
        mapvalues(clinical_dataframe$TNBC_status,
                  c(TRUE,
                    FALSE),
                  c(1, 0))
      
      # Change to factor
      clinical_dataframe$TNBC_status <-
        as.factor(clinical_dataframe$TNBC_status)
      
    } else if (feature == "pathological_history") {
      clinical_dataframe$pathological_history <-
        mapvalues(
          clinical_dataframe$pathological_history,
          c(
            "Neg",
            "Hypothyroidism",
            "Diabetes",
            "Ocular Surgery",
            "3 C sections",
            "Hysterechtomy",
            "Hypertension and Diabetes"
          ),
          c(0, 1, 1, 1, 1, 1, 1)
        )
      
      # Change to factor
      clinical_dataframe$pathological_history <-
        as.factor(clinical_dataframe$pathological_history)
      
    } else if (feature == "family_pathological_hsitory") {
      clinical_dataframe$family_pathological_hsitory <-
        mapvalues(
          clinical_dataframe$family_pathological_hsitory,
          c(
            "Neg",
            "Maternal Aunt Breast Cancer",
            "Diabetes",
            "Liver Cancer",
            "Father Pancreas Cancer. Maternal Aunt Breast Cancer. Father Aunts Ovarian Cancer. Cousin Breast Cancer",
            "Father Laringeal Cancer",
            "Fater Liver Cancer",
            "Mother and Sister Ovarian Cancer"
          ),
          c(0, 1, 1, 1, 1, 1, 1, 1)
        )
      
      # Change to factor
      clinical_dataframe$family_pathological_hsitory <-
        as.factor(clinical_dataframe$family_pathological_hsitory)
      
    } else if (feature == "radiotherapy") {
      clinical_dataframe$radiotherapy <-
        mapvalues(clinical_dataframe$radiotherapy,
                  c("Yes",
                    "No",
                    "N/A",
                    "NO",
                    "YES"),
                  c(1, 0, NA, 0, 1))
      # Change to factor
      clinical_dataframe$radiotherapy <-
        as.factor(clinical_dataframe$radiotherapy)
      
    } else if (feature == "hormone_therapy") {
      clinical_dataframe$hormone_therapy <-
        mapvalues(clinical_dataframe$hormone_therapy,
                  c("NO",
                    "YES"),
                  c(0, 1))
      # Change to factor
      clinical_dataframe$hormone_therapy <-
        as.factor(clinical_dataframe$hormone_therapy)
      
    } else if (feature == "methastasis") {
      clinical_dataframe$methastasis <-
        mapvalues(clinical_dataframe$methastasis,
                  c("no",
                    "yes",
                    "0",
                    "1"),
                  c(0, 1, 0, 1))
      # Change to factor
      clinical_dataframe$methastasis <-
        as.factor(clinical_dataframe$methastasis)
      
    } else if (feature == "relapse_free_status") {
      clinical_dataframe$relapse_free_status <-
        mapvalues(clinical_dataframe$relapse_free_status,
                  c("0:Not Recurred",
                    "1:Recurred"),
                  c(0, 1))
      # Change to numeric
      clinical_dataframe$relapse_free_status <-
        as.numeric(clinical_dataframe$relapse_free_status)
      
    } else if (feature == "tumor_stage" || feature == "Stage") {
      clinical_dataframe[,feature] <-
        mapvalues(
          clinical_dataframe[,feature],
          c("N/A",
            "UNKNOWN",
            "3B",
            "2B",
            "1",
            "2",
            "3",
            "4",
            "IIIB", 
            "IIB",
            "IIA",
            "IIIC",
            "IIIA",
            "II",
            "III"),
          c(NA, NA, 3, 2, 1, 2, 3, 4, 3, 2, 2, 3, 3, 2, 3)
        )
      # Change to factor
      clinical_dataframe[,feature] <-
        as.factor(clinical_dataframe[,feature])
      
    } else if (feature == "race" || feature == "Race") {
      clinical_dataframe[, feature] <-
        mapvalues(
          clinical_dataframe[, feature],
          c(
            "W",
            "B",
            "H",
            "UNKNOWN",
            "WHITE",
            "BLACK OR AFRICAN AMERICAN",
            "OTHER",
            "African American",
            "Caucasian"       
          ),
          c(
            "white",
            "african american",
            "hispanic or latino",
            NA,
            "white",
            "african american",
            NA,
            "african american",
            "white"
          )
        )
    } else if (feature == "age" || feature == "Age") {
      # Change age to numeric
      clinical_dataframe[,feature] <- as.numeric(clinical_dataframe[,feature])
      
      # Add a new variant of the age as ranges
      clinical_dataframe$age_label <- clinical_dataframe[,feature]
      
      for (age_row in 1:nrow(clinical_dataframe)) {
        if (is.na(clinical_dataframe[age_row, feature]) == FALSE) {
          if (clinical_dataframe[age_row, feature] < 40) {
            clinical_dataframe[age_row, "age_label"] <- "30s"
          }
          if (clinical_dataframe[age_row, feature] >= 40 &
              clinical_dataframe[age_row, feature] < 50) {
            clinical_dataframe[age_row, "age_label"] <- "40s"
          }
          if (clinical_dataframe[age_row, feature] >= 50 &
              clinical_dataframe[age_row, feature] < 60) {
            clinical_dataframe[age_row, "age_label"] <- "50s"
          }
          if (clinical_dataframe[age_row, feature] >= 60 &
              clinical_dataframe[age_row, feature] < 70) {
            clinical_dataframe[age_row, "age_label"] <- "60s"
          }
          if (clinical_dataframe[age_row, feature] >= 70 &
              clinical_dataframe[age_row, feature] < 80) {
            clinical_dataframe[age_row, "age_label"] <- "70s"
          }
          if (clinical_dataframe[age_row, feature] >= 80 &
              clinical_dataframe[age_row, feature] < 90) {
            clinical_dataframe[age_row, "age_label"] <- "80s"
          }
          if (clinical_dataframe[age_row, feature] >= 90) {
            clinical_dataframe[age_row, "age_label"] <- "90s"
          }
        }
      }
    } else if (feature == "tumor_size_mm") {
      # Change tumor size to numeric
      clinical_dataframe$tumor_size_mm <-
        as.numeric(clinical_dataframe$tumor_size_mm)
      # Add a new variant of the tumor size as ranges
      clinical_dataframe$tumor_size <- clinical_dataframe$tumor_size_mm
      
      for (age_row in 1:nrow(clinical_dataframe)) {
        if (is.na(clinical_dataframe[age_row, "tumor_size_mm"]) == FALSE) {
          if (clinical_dataframe[age_row, "tumor_size_mm"] < 10) {
            clinical_dataframe[age_row, "tumor_size"] <- "Below 10 mm"
          }
          if (clinical_dataframe[age_row, "tumor_size_mm"] >= 10 &
              clinical_dataframe[age_row, "tumor_size_mm"] < 20) {
            clinical_dataframe[age_row, "tumor_size"] <- "Between 10 and 19 mm"
          }
          if (clinical_dataframe[age_row, "tumor_size_mm"] >= 20 &
              clinical_dataframe[age_row, "tumor_size_mm"] < 30) {
            clinical_dataframe[age_row, "tumor_size"] <- "Between 20 and 29 mm"
          }
          if (clinical_dataframe[age_row, "tumor_size_mm"] >= 30 &
              clinical_dataframe[age_row, "tumor_size_mm"] < 40) {
            clinical_dataframe[age_row, "tumor_size"] <- "Between 30 and 39 mm"
          }
          if (clinical_dataframe[age_row, "tumor_size_mm"] >= 40 &
              clinical_dataframe[age_row, "tumor_size_mm"] < 50) {
            clinical_dataframe[age_row, "tumor_size"] <- "Between 40 and 49 mm"
          }
          if (clinical_dataframe[age_row, "tumor_size_mm"] >= 50 &
              clinical_dataframe[age_row, "tumor_size_mm"] < 60) {
            clinical_dataframe[age_row, "tumor_size"] <- "Between 50 and 59 mm"
          }
          if (clinical_dataframe[age_row, "tumor_size_mm"] >= 60 &
              clinical_dataframe[age_row, "tumor_size_mm"] < 70) {
            clinical_dataframe[age_row, "tumor_size"] <- "Between 60 and 69 mm"
          }
          if (clinical_dataframe[age_row, "tumor_size_mm"] >= 70 &
              clinical_dataframe[age_row, "tumor_size_mm"] < 80) {
            clinical_dataframe[age_row, "tumor_size"] <- "Between 70 and 79 mm"
          }
          if (clinical_dataframe[age_row, "tumor_size_mm"] >= 80 &
              clinical_dataframe[age_row, "tumor_size_mm"] < 90) {
            clinical_dataframe[age_row, "tumor_size"] <- "Between 80 and 89 mm"
          }
          if (clinical_dataframe[age_row, "tumor_size_mm"] >= 90) {
            clinical_dataframe[age_row, "tumor_size"] <- "Greater than 89 mm"
          }
        }
      }
      
    } else if (feature == "days_to_last_followup" ||
               feature == "days_to_death") {
      clinical_dataframe$Surv_event <- NA
      clinical_dataframe$Surv_days <- NA
      
      for (i in 1:nrow(clinical_dataframe)) {
        if (is.na(clinical_dataframe$days_to_last_followup[i]) != T) {
          clinical_dataframe$Surv_days[i] <-
            clinical_dataframe$days_to_last_followup[i]
          clinical_dataframe$Surv_event[i] <- 0
        }
        if (is.na(clinical_dataframe$days_to_death[i]) != T) {
          clinical_dataframe$Surv_days[i] <-
            clinical_dataframe$days_to_death[i]
          clinical_dataframe$Surv_event[i] <- 1
        }
        if (is.na(clinical_dataframe$days_to_last_followup[i]) == T &&
            is.na(clinical_dataframe$days_to_death[i]) == T &&
            is.na(clinical_dataframe$vital_status[i]) != T &&
            is.na(clinical_dataframe$survival_months[i]) != T) {
          if (is.na(clinical_dataframe$survival_months[i]) != T) {
            clinical_dataframe$Surv_days[i] <-
              clinical_dataframe$survival_months[i] / 30
          }
          if (clinical_dataframe$vital_status[i] == "alive" ||
              clinical_dataframe$vital_status[i] == "Living") {
            clinical_dataframe$Surv_event[i] <- 0
          } else if (clinical_dataframe$vital_status[i] == "dead" ||
                     clinical_dataframe$vital_status[i] == "Died of Other Causes" ||
                     clinical_dataframe$vital_status[i] == "Died of Disease") {
            clinical_dataframe$Surv_event[i] <- 1
          }
        }
      }
    }
  }
  
  clinical_dataframe <-
    clinical_dataframe[, !(names(clinical_dataframe) %in% c("last_known_alive"))]
  if ("relapse_free_status" %in% features_to_engineer) {
    clinical_dataframe$relapse_free_status <-
      as.factor(clinical_dataframe$relapse_free_status)
  }
  
  return(clinical_dataframe)
}

# Function to create a consensus cluster heatmap plot per resulting clusters based
# on the consensus clustering results
consensus_complex_heatmap = function(res_train, k_clusters_train) {
  ccl <- list()
  ha_row <- list()
  x <- c(
    "skyblue",
    "gold",
    "violet",
    "darkorchid",
    "slateblue",
    "forestgreen",
    "violetred",
    "orange",
    "midnightblue",
    "grey31",
    "black"
  )
  
  names(x) <- as.character(seq(1, 11, by = 1))
  for (i in seq(2, k_clusters_train)) {
    # get cc matrix and labels
    ccmatrix <- res_train$realdataresults[[i]]$consensus_matrix
    db_source <- as.data.frame(RNA_clin_data[
      c(rownames(res_train$realdataresults[[i]]$ordered_annotation)),
      "source_db"])
    rownames(db_source) <-
      c(rownames(res_train$realdataresults[[i]]$ordered_annotation))
    colnames(db_source) <- c("db_source")
    annon <- res_train$realdataresults[[i]]$ordered_annotation
    
    # do heatmap
    n <- 10
    seq <- rev(seq(0, 255,
                   by = 255 / (n)))
    palRGB <- cbind(seq,
                    seq, 255)
    mypal <- rgb(palRGB,
                 maxColorValue = 255)
    ha = HeatmapAnnotation(
      df = data.frame(Cluster = as.character(annon[, 1])),
      Source = as.character(db_source[, 1]),
      col = list(
        Cluster = x,
        Source = c(
          "TCGA_RNAseq" = "green",
          "Ref2_RNAseq" = "red",
          "RefGSE20271_Array" = "blue",
          "RefGSE575678_Array" = "pink",
          "Ref4_RNAseq" = "yellow"
        )
      )
    )
    ccl[[i]] <- Heatmap(
      ccmatrix,
      name = "Consensus_index",
      top_annotation = ha,
      col = mypal,
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_column_names = FALSE
    )
  }
  return(ccl)
}

clinical_plots = function(RNA_clin_data, clinical_vars_names, count_or_percentage, k_clusters) {
  # We split the clinical data in cluster 1, cluster 2 and cluster 3
  clin_data1 <-
    RNA_clin_data[which(RNA_clin_data["consensuscluster"] == 1),]
  clin_data2 <-
    RNA_clin_data[which(RNA_clin_data["consensuscluster"] == 2),]
  if (k_clusters == 3) {
    clin_data3 <-
      RNA_clin_data[which(RNA_clin_data["consensuscluster"] == 3),]
  }
  
  # Plotting categorical and numeric variables per cluster
  if (count_or_percentage == "counts") {
    for (vars in clinical_vars_names) {
      var1 <- as.data.frame(table(clin_data1[vars]))
      cluster_info <- as.data.frame(rep("1",
                                        nrow(var1)))
      var1 <- cbind(cluster_info,
                    var1)
      colnames(var1) <- c("cluster_info",
                          "variable",
                          "count")
      
      
      var2 <- as.data.frame(table(clin_data2[vars]))
      cluster_info <- as.data.frame(rep("2",
                                        nrow(var2)))
      var2 <- cbind(cluster_info, var2)
      colnames(var2) <- c("cluster_info",
                          "variable",
                          "count")
      if (k_clusters == 3) {
        var3 <- as.data.frame(table(clin_data3[vars]))
        cluster_info <- as.data.frame(rep("3",
                                          nrow(var3)))
        var3 <- cbind(cluster_info, var3)
        colnames(var3) <- c("cluster_info",
                            "variable",
                            "count")
      }
      
      if (k_clusters == 3) {
        merged_df <- rbind(var1,
                           var2,
                           var3)
      }else if (k_clusters == 2) {
        merged_df <- rbind(var1,
                           var2)
      }
      
      ## PLOTS
      # Grouped barcharts
      plt <- ggplot(merged_df,
                    aes(x = cluster_info,
                        y = count,
                        fill = variable)) +
        geom_bar(position = "dodge",
                 stat = "identity")+ 
        stat_compare_means()+
        labs(title = as.character(vars))
      
      if (vars %in% c("age",
                      "tumor_size_mm",
                      "age_at_menarche",
                      "Surv_months",
                      "Age")) {
        
        if (k_clusters == 3) {
          merged_df <- rbind(clin_data1[c("consensuscluster", vars)],
                             clin_data2[c("consensuscluster", vars)],
                             clin_data3[c("consensuscluster", vars)])
        }else if (k_clusters == 2) {
          merged_df <- rbind(clin_data1[c("consensuscluster", vars)],
                             clin_data2[c("consensuscluster", vars)])
        }
        colnames(merged_df) <- c("consensuscluster",
                                 "vars")
        
        plt <- ggboxplot(merged_df, 
                         x="consensuscluster",
                         y="vars",
                         color = "consensuscluster",
                         ggtheme = theme_gray(),
                         add='jitter') +
          stat_compare_means() +
          labs(title = as.character(vars)) +
          xlab("Cluster") + ylab(vars)
        
      }
      print(plt)
    }
  }else if (count_or_percentage == "percentage") {
    for (vars in clinical_vars_names) {
      var1 <- as.data.frame(table(clin_data1[vars]))
      cluster_info <- as.data.frame(rep("1",
                                        nrow(var1)))
      var1 <- cbind(cluster_info,
                    var1)
      
      var1_perc <- var1
      var1_perc$Freq <-  (var1$Freq/sum(var1$Freq))
      
      colnames(var1) <- c("cluster_info",
                          "variable",
                          "count")
      colnames(var1_perc) <- c("cluster_info",
                               "variable",
                               "percentage")
      
      
      var2 <- as.data.frame(table(clin_data2[vars]))
      cluster_info <- as.data.frame(rep("2",
                                        nrow(var2)))
      var2 <- cbind(cluster_info, var2)
      
      var2_perc <- var2
      var2_perc$Freq <-  (var2$Freq/sum(var2$Freq))
      
      colnames(var2) <- c("cluster_info",
                          "variable",
                          "count")
      colnames(var2_perc) <- c("cluster_info",
                               "variable",
                               "percentage")
      if (k_clusters == 3) {
        var3 <- as.data.frame(table(clin_data3[vars]))
        cluster_info <- as.data.frame(rep("3",
                                          nrow(var3)))
        var3 <- cbind(cluster_info, var3)
        
        var3_perc <- var3
        var3_perc$Freq <-  (var3$Freq/sum(var3$Freq))
        
        colnames(var3) <- c("cluster_info",
                            "variable",
                            "count")
        colnames(var3_perc) <- c("cluster_info",
                                 "variable",
                                 "percentage")
        merged_df <- rbind(var1_perc,
                           var2_perc,
                           var3_perc)
      }else if (k_clusters == 2) {
        merged_df <- rbind(var1_perc,
                           var2_perc)
      }
      
      ## PLOTS
      # Grouped barcharts
      plt <- ggplot(merged_df,
                    aes(x = cluster_info,
                        y = percentage,
                        fill = variable)) +
        geom_bar(position = "dodge",
                 stat = "identity")+ 
        stat_compare_means()+
        labs(title = as.character(vars)) + scale_y_continuous(labels = scales::percent)
      
      if (vars %in% c("age",
                      "tumor_size_mm",
                      "age_at_menarche",
                      "Surv_months",
                      "Age")) {
        if (k_clusters == 3) {
          merged_df <- rbind(clin_data1[c("consensuscluster", vars)],
                             clin_data2[c("consensuscluster", vars)],
                             clin_data3[c("consensuscluster", vars)])
        }else if (k_clusters == 2) {
          merged_df <- rbind(clin_data1[c("consensuscluster", vars)],
                             clin_data2[c("consensuscluster", vars)])
        }
        colnames(merged_df) <- c("consensuscluster",
                                 "vars")
        
        plt <- ggboxplot(merged_df, 
                         x="consensuscluster",
                         y="vars",
                         color = "consensuscluster",
                         ggtheme = theme_gray(),
                         add='jitter') +
          stat_compare_means() +
          labs(title = as.character(vars)) +
          xlab("Cluster") + ylab(vars)
        
      }
      
      print(plt)
    }
  }
}

## Statistical tests functions

# Function to perform chi-squared test for the categorical variables of the 
# clusters and identify which labels are significant.
chi_squared_clusters <- function(variable_name, clinical_df, k_clusters) {
  ## Saving the possible outcomes of the feature
  mut_chisq_test_names <- c()
  variable_to_check <- variable_name
  values_to_check <- unique(clinical_df[,variable_name])
  values_to_check <- values_to_check[is.na(values_to_check) != TRUE]
  
  final_results <- c()
  
  for (value in c(values_to_check)) {
    cluster1_no <-
      length(which(
        clinical_df[, variable_to_check] == value &
          clinical_df[, "consensuscluster"] == 1 &
          is.na(clinical_df[, variable_to_check]) != TRUE
      ))
    cluster2_no <-
      length(which(
        clinical_df[, variable_to_check] == value &
          clinical_df[, "consensuscluster"] == 2 &
          is.na(clinical_df[, variable_to_check]) != TRUE
      ))
    
    if (k_clusters == 3) {
      cluster3_no <-
        length(which(
          clinical_df[, variable_to_check] == value &
            clinical_df[, "consensuscluster"] == 3 &
            is.na(clinical_df[, variable_to_check]) != TRUE
        ))
    }
    
    
    cluster1 <-
      length(which(
        clinical_df[, variable_to_check] != value &
          clinical_df[, "consensuscluster"] == 1 &
          is.na(clinical_df[, variable_to_check]) != TRUE
      ))
    cluster2 <-
      length(which(
        clinical_df[, variable_to_check] != value &
          clinical_df[, "consensuscluster"] == 2 &
          is.na(clinical_df[, variable_to_check]) != TRUE
      ))
    
    if (k_clusters == 3) {
      cluster3 <-
        length(which(
          clinical_df[, variable_to_check] != value &
            clinical_df[, "consensuscluster"] == 3 &
            is.na(clinical_df[, variable_to_check]) != TRUE
        ))
    }
    
    if (k_clusters == 3) {
      df <-
        data.frame(
          cluster1 = c(0, 0),
          cluster2 = c(0, 0),
          cluster3 = c(0, 0)
        )
      
      
      df$cluster1 <- c(cluster1_no, cluster1)
      df$cluster2 <- c(cluster2_no, cluster2)
      df$cluster3 <- c(cluster3_no, cluster3)
    }else if (k_clusters == 2) {
      df <-
        data.frame(
          cluster1 = c(0, 0),
          cluster2 = c(0, 0)
        )
      
      
      df$cluster1 <- c(cluster1_no, cluster1)
      df$cluster2 <- c(cluster2_no, cluster2)
    }
    
    rownames(df) <- c("0", "1")
    
    chisq_test_var <- chisq.test(df, simulate.p.value = T)
    mut_chisq_test <- c(chisq_test_var$p.value)
    
    ## CALCULATING EFFECT SIZE
    # Custom function
    Cramers_V <- function(chi, n, df_es)
      sqrt((chi) / (n * df_es))
    # Find degrees of freedom - min row or col - 1
    df_es <- min(dim(df)) - 1
    # Calculate
    effect_size <-
      Cramers_V(chi = chisq_test_var$statistic,
                n = sum(df),
                df = df_es)
    
    # we save the P values and FDR of each test
    mut_chisq_test <-
      cbind(mut_chisq_test,
            p.adjust(mut_chisq_test, method = 'fdr'),
            effect_size,
            df_es)
    colnames(mut_chisq_test) <-
      c("p_value", "fdr", "effect_size", "DOF")
    
    # we order the resulting matrix by P value
    if (mut_chisq_test[1]<0.05) {
      final_results[[value]] <- mut_chisq_test
    }
  }
  return(final_results)
}

#Function to perform the Kruskal-Wallis test on the categorical variables and get 
#a global value of significance. At the same time, the chi-square function is run 
#to understand which labels are the most significant.
categorical_kruskal_wallis_test <-
  function(variable_name, clinical_df, formula, k_clusters) {
    print(variable_name)
    print("---")
    print("Calculating the frequencies per cluster...")
    # Calculate the frequencies of each class into each cluster
    temp_table <-
      table(clinical_df[c(variable_name, "consensuscluster")])
    
    # counts
    print("counts")
    print(temp_table)
    
    # percentage
    print("percentages")
    print(colPercents(temp_table, 2))
    
    # --------------
    # Kruskal-Wallis test (is the variable significant?)
    kruskal_wallis_general <-
      kruskal.test(formula, data = clinical_df)
    print(kruskal_wallis_general)
    
    # Chi-squared test
    chisq_test_res <-
      chi_squared_clusters(variable_name, clinical_df, k_clusters)
    
    return(chisq_test_res)
  }

# Function to perform logistic regression on the numeric variables and return 
# the significance of each pair of clusters variables.
clinical_test <-
  function(variable_name,
           clinical_df,
           formula,
           predictor_type,
           k_clusters) {
    print(variable_name)
    print("---")
    print("Running logistic regression per each cluster combination...")
    
    # Removing rows with missing values
    clinical_df_nonull <-
      clinical_df[!is.na(clinical_df[, variable_name]),]
    
    # Logistic regression on each possible cluster combination
    significant_excluding <- c()
    for (cluster in 1:k_clusters) {
      model <-
        glm(formula = formula,
            data = clinical_df_nonull[clinical_df_nonull$consensuscluster != cluster, ],
            family = binomial)
      significant_excluding[[as.numeric(cluster)]] <-
        cbind(summary(model)$coeff[-1, 4], summary(model)$coeff[-1, 4] < 0.05)
    }
    print(significant_excluding)
    if (predictor_type == "categorical") {
      print("Running Kruskal Wallis test and Chi-square test for independence...")
      categorical_kruskal_wallis_test(
        variable_name = variable_name,
        clinical_df = clinical_df_nonull,
        formula = formula,
        k_clusters = k_clusters
      )
      
      # Table of frequencies
      temp_table <-
        table(clinical_df[!is.na(clinical_df[, variable_name]), c(variable_name, "consensuscluster")])
      
      # counts
      temp_table
      
      # percentage
      colPercents(temp_table, 2)
      
      # Plots
      tab_var <-
        table(clinical_df[, c("consensuscluster", variable_name)])
      
      if (k_clusters == 3) {
        plt <-
          barplot(
            tab_var,
            legend.text = c('cluster1', 'cluster2', 'cluster3'),
            xlab = variable_name,
            ylab = "counts"
          )
      }else if (k_clusters == 2) {
        plt <-
          barplot(
            tab_var,
            legend.text = c('cluster1', 'cluster2'),
            xlab = variable_name,
            ylab = "counts"
          )
      }
      
      print(plt)
      
      plt_prop <- ggplot(clinical_df,
                         aes_string(x = "consensuscluster",
                             fill = variable_name)) +
        geom_bar(position = "fill") +
        labs(y = "Proportion")
      print(plt_prop)
    }
    else if (predictor_type == "numerical") {
      print("Running ANOVA test...")
      
      aov_variable <-
        aov(clinical_df_nonull[, variable_name] ~ clinical_df_nonull$consensuscluster)
      print(variable_name)
      print(summary(aov_variable))
      
      
      # Plots
      group_by(clinical_df[!is.na(clinical_df[, variable_name]), ], clinical_df[!is.na(clinical_df[, variable_name]), "consensuscluster"]) %>%
        summarise(
          count = n(),
          mean = mean(variable_name, na.rm = TRUE),
          sd = sd(variable_name, na.rm = TRUE)
        )
      
      ggplot(clinical_df, aes_string(y = variable_name)) +
        geom_boxplot(
          aes_string(color = "consensuscluster", fill = "consensuscluster"),
          bins = 30,
          alpha = 0.3,
          outlier.size = 2,
          outlier.colour = "black"
        )
    }
  }


# Function for mutation analysis
mutation_analysis <- function(merged_mut, RNA_seq_merged_allgenes, k_clusters, annon) {
  # We drop mutations related to: UTR, introns, RNA and silent.
  mut_filters <-
    c("3'UTR",
      "5'Flank",
      "3'Flank",
      "5'UTR",
      "IGR",
      "Intron",
      "Silent",
      "RNA")
  mut_filter_index <-
    which(merged_mut$Variant_Classification %in% mut_filters)
  merged_mut <- merged_mut[-mut_filter_index, ]
  
  # We obtain the different samples and genes to make an array where it will store 0 and 1 depending on whether
  # mutation is present.
  u_samples <- sort(unique(merged_mut$Tumor_Sample_Barcode))
  u_genes <- sort(unique(merged_mut$Hugo_Symbol))
  mut_matrix <-
    matrix(
      0,
      nrow = length(u_genes),
      ncol = length(u_samples),
      
      dimnames = list(u_genes, u_samples)
    )
  mut_matrix <- mut_matrix[, -1]
  
  # We fill the array by doing one loop per sample
  for(i in colnames(mut_matrix)) {
    aux_genes <-
      merged_mut$Hugo_Symbol[which(merged_mut$Tumor_Sample_Barcode == i)]
    mut_matrix[aux_genes, i] <- 1
  }
  
  # Defining the cluster groups with the patient IDs
  cluster = list()
  
  for (i in seq(1, k_clusters)) {
    cluster[[i]] <- filter(annon, annon$consensuscluster == i)
  }
  
  #We observe the difference in the naming convention of mutations and rna-seq data, as well as the difference in dimensions. We can compare the first 16 characters.
  raw_brca_rnaseq <- RNA_seq_merged_allgenes
  
  mut_rnaseq_match <- match(gsub(
    pattern = "-",
    replacement = "\\.",
    x = substr(colnames(mut_matrix), 1, 16)
  ),
  substr(rownames(raw_brca_rnaseq), 1, 16))
  
  #Clean the mutation matrix and expression matrix from NA values
  raw_brca_rnaseq <- t(RNA_seq_merged_allgenes)
  
  # to store the non NA values postitions
  mut_rnaseq_match_nonNA <- which(!is.na(mut_rnaseq_match))
  
  # to remove the NA values from mut_matrix_rna and mut_rnaseq_match
  mut_matrix_rna <- mut_matrix[, mut_rnaseq_match_nonNA]
  mut_rnaseq_match <- mut_rnaseq_match[mut_rnaseq_match_nonNA]
  
  rna_mut <- raw_brca_rnaseq[, mut_rnaseq_match]
  
  # We filter the genes in RNA-seq that show no variation and we keep the most frequent mutations only.
  # We convert all the dataframe as numeric
  rna_mut <-
    mutate_all(data.frame(rna_mut), function(x)
      as.numeric(as.character(x)))
  
  # We apply the mad function
  rna_mut_mad <- apply(rna_mut, 1, mad)
  mad_zero <- which(rna_mut_mad == 0)
  rna_mut <- rna_mut[-mad_zero, ]
  
  # We count the times each mutation appears in the matrix and store the ones with an appearance of >9 times
  # We count the times each mutation appears in the matrix
  mut_matrix_rna_count <- apply(mut_matrix_rna, 1, sum)
  
  # if the mutation count is greater than 9, we add it to the mut_matrix_rna
  mut_matrix_rna <- mut_matrix_rna[which(mut_matrix_rna_count > 9), ]
  
  rna_results <- list("rna_mut" = rna_mut, "mut_matrix_rna" = mut_matrix_rna)
  
  return(rna_results)
}

mutation_t_test <- function(mut_matrix_rna, rna_mut) {
  mut_rna_test <- c()
  mut_rna_test_names <- c()
  for (i in 1:nrow(mut_matrix_rna)) {
    # we apply the function "function(x)" to each of the row values of rna_mut where each row value is x
    aux_mut_test <- t(sapply(1:nrow(rna_mut), function(x) {
      # mutation: the mutation row (presence or absence of mutation per gene sequence)
      # exp: the gene row of the rna_mut matrix
      temp_df <-  data.frame(mut = mut_matrix_rna[i,],
                             exp = as.matrix(rna_mut)[x,])
      temp_df <- na.omit(temp_df)
      
      if (length(temp_df) > 0 && length(unique(temp_df$mut)) == 2) {
        aux_test <-
          t.test(
            formula = exp ~ mut,
            data = temp_df
          )
      }
    }))
    
    mut_rna_test_names <- c(sapply(1:nrow(rna_mut), function(x) {
      # mutation: the mutation row (presence or absence of mutation per gene sequence)
      # exp: the gene row of the rna_mut matrix
      temp_df <-  data.frame(mut = mut_matrix_rna[i,],
                             exp = as.matrix(rna_mut)[x,])
      temp_df <- na.omit(temp_df)
      
      if (length(temp_df) > 0 && length(unique(temp_df$mut)) == 2) {
        mut_rna_test_names <- paste(rownames(mut_matrix_rna)[i], rownames(rna_mut),sep = "_")
      }
    }))
    
    # appending each test into a vector
    mut_rna_test <- c(mut_rna_test, aux_mut_test)
  }
  
  # we save the P values and FDR of each test
  mut_rna_test <- unlist(mut_rna_test)
  mut_rna_test <- cbind(mut_rna_test, p.adjust(mut_rna_test, method='fdr'))
  colnames(mut_rna_test) <- c("p_value", "fdr")
  rownames(mut_rna_test) <- mut_rna_test_names
  
  # we order the resulting matrix by P value
  mut_rna_test <- mut_rna_test[order(mut_rna_test[, "p_value"]), ]
  
  return(mut_rna_test)
}

# Function to create the cluster slices of the mutations and rnaseq matrices
create_cluster_rnaseq_matrix <- function(rna_mut, mut_matrix_rna, cluster_num) {
  mut_rnaseq_match_temp <-
    match(gsub(
      pattern = "-",
      replacement = "\\.",
      x = substr(colnames(rna_mut), 1, 16)
    ),
    substr(colnames(RNA_seq_merged[, which(RNA_seq_merged["consensuscluster",] == cluster_num)]), 1, 16))
  
  # to store the NA values postitions
  NA_cluster_temp <- which(is.na(mut_rnaseq_match_temp))
  
  # to remove the NA values from mut_matrix_rna and mut_rnaseq_match
  mut_matrix_rna_temp <- mut_matrix_rna[,-NA_cluster_temp]
  mut_rnaseq_match_temp <- mut_rnaseq_match_temp[-NA_cluster_temp]
  rna_mut_cluster_temp <- rna_mut[-NA_cluster_temp]
  
  results_cluster_rnaseq <- list("mut_matrix_rna" = mut_matrix_rna_temp, "mut_rnaseq_match" = mut_rnaseq_match_temp, "rna_mut_cluster" = rna_mut_cluster_temp)
  
  return(results_cluster_rnaseq)
}

chi2_test_clusters_pair <- function(mut_matrix_rna, mut_matrix_rna_cluster1, mut_matrix_rna_cluster2, clusters_num_list) {
  mut_rna_test <- c()
  mut_rna_test_names <- c()
  for (i in 1:nrow(mut_matrix_rna)) {
    # mutation: the mutation row (presence or absence of mutation per gene sequence)
    # exp: the gene row of the rna_mut matrix
    ones1 <- sum(mut_matrix_rna_cluster1[i,])
    zeros1 <- length(mut_matrix_rna_cluster1[i,]) - ones1
    
    ones2 <- sum(mut_matrix_rna_cluster2[i,])
    zeros2 <- length(mut_matrix_rna_cluster2[i,]) - ones2
    
    mut_clust <- as.table(rbind(c(zeros1, zeros2), c(ones1, ones2)))
    dimnames(mut_clust) <- list(status = c("No mutation", "Mutation"),
                                cluster = clusters_num_list)
    
    aux_mut_test <-
      chisq.test(mut_clust, simulate.p.value = T)$p.value
    
    # appending each test into a vector
    mut_rna_test <- c(mut_rna_test, aux_mut_test)
    # appending each name into a vector
    mut_rna_test_names <-
      c(mut_rna_test_names, paste(rownames(mut_matrix_rna)[i]))
  }
  
  # we save the P values and FDR of each test
  mut_rna_test <- cbind(mut_rna_test, p.adjust(mut_rna_test, method = 'fdr'))
  colnames(mut_rna_test) <- c("p_value", "fdr")
  rownames(mut_rna_test) <- mut_rna_test_names
  
  # we order the resulting matrix by P value
  mut_rna_test_12 <- as.data.frame(mut_rna_test[order(mut_rna_test[, "p_value"]),])
  
  print("Significant mutations")
  significant_mutations <-
    mut_rna_test_12[mut_rna_test_12$p_value < 0.05, ]
  significant_mutations_names12 <- rownames(significant_mutations)
  
  significant_mutations_values <-
    as.data.frame(cbind(
      rowSums(mut_matrix_rna_cluster1[significant_mutations_names12, ]),
      rowSums(mut_matrix_rna_cluster2[significant_mutations_names12, ])
    ))
  
  colnames(significant_mutations_values) <- clusters_num_list
  
  chi2_mutations <- list("significant_mutations_values" = significant_mutations_values, "mut_rna_test" = mut_rna_test)
  
  return(chi2_mutations)
}

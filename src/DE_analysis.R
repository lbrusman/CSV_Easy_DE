suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(janitor))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ROCR))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(ggcorrplot))
suppressPackageStartupMessages(library(RUVSeq))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(EDASeq))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(R.utils))

#get arguments from command line
args = commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1]

keys <- attachLocally(args)
cat("Command-line arguments attached to global environment:\n");
# print(keys);
str(mget(keys, envir=globalenv()))

set.seed(123)

meta <- read.table(metadata, sep = "\t", header = TRUE)
colnames(meta) <- gsub("_", "", colnames(meta))

files <- list.files(indir, pattern = file_pattern, full.names = TRUE)
og_celltypes <- list.files(indir, pattern = file_pattern, full.names = FALSE, recursive = FALSE) %>% str_split_i(file_pattern, 1)
celltypes_fix <- gsub("_", "", og_celltypes)

# print(meta$Sex)
# print(meta["Sex"])
# print("here is wtf is going on with subset_on1")

#let's get list of n_cells per cell type
meta[,celltypes_fix] <- sapply(meta[,celltypes_fix], as.numeric)
all_n_cells <- colSums(meta[,celltypes_fix], na.rm = TRUE)
names(all_n_cells) <- celltypes_fix

covariates_list <- str_split(covariates, pattern = "_") %>% unlist()

deseq_formula <- paste0("~", paste0(covariates_list, collapse = "+"), "+", contrast_var)
print(deseq_formula)

#get variables for correlation with latent vars
correlate_vars <- str_split(correlate_vars, pattern = "_") %>% unlist()

# Setup additional variables
max.k.global = 10 # Roughly 1/2 samples is ok

#get RLE plot thresholds for RNA or ATAC
if (assay == "RNA") {
  f_tresh = 10
}
if (assay == "ATAC") {
  f_tresh = 100
}

# Don't touch this one
latent_vars <- paste0("W_", 1:max.k.global)

# Write here the covariants you want to test against the Ws and EXCLUDE from the final formula (Usually just disease)
covariants <- contrast_var

#get subsets out
subset1 <- str_split(subset1, pattern = "_") %>% unlist()
print("here is split subset1")
print(subset1)
subset2 <- str_split(subset2, pattern = "_") %>% unlist()

#do RUVseq for T1D vs. non-diabetic
deseq_stats = data.frame()
lm_res_all = data.frame()

for (i in og_celltypes) {
    print(i)
    file.use <- paste0(indir, i, file_pattern)
    print(file.use)
    cell.use <- gsub("_", "", i)
    print(cell.use)
    print(paste0("Analyzing: ", cell.use))
    raw_counts <- read.table(file.use, header = TRUE, row.names = 1)

    print(paste0("Donors detectected with celltype: ", ncol(raw_counts)))

    meta_cell <- meta
    meta_cell[,donor_col] <- ifelse(grepl("^[[:digit:]]+", meta_cell[,donor_col])==TRUE, paste0("X", meta_cell[,donor_col]), meta_cell[,donor_col])
    meta_cell[,donor_col] <- gsub("-", "\\.", meta_cell[,donor_col])
    
    meta_cell <- meta_cell[meta_cell[,donor_col] %in% colnames(raw_counts),]
    rownames(meta_cell) = meta_cell[,donor_col]
    # Ensure that the column names of raw_counts are in the same order as the library identifiers in meta_cell
    meta_cell <- meta_cell[colnames(raw_counts),]

    if (!is.na(subset_on1) & subset_on1 != "NA") { #changed from subset_on1 != "NA"
      print("trying to subset1")
      meta_cell <- meta_cell[which(meta_cell[,subset_on1] %in% subset1),]
      print("here's what's left of meta_cell")
      print(head(meta_cell))
    }

    if (subset_on2 != "NA") {
      meta_cell <- meta_cell[which(meta_cell[,subset_on2] %in% subset2),]
    }
        
    Analysis = paste0(cell.use, "_", contrast_id)
    print(paste0("Analysis: ", Analysis))
    #first let's get rid of NA for that var
    meta_cell <- meta_cell[!is.na(meta_cell[,contrast_var]),]

    if (is.logical(meta_cell[,contrast_var])) {
        print("column is logical! changing to character")
        meta_cell[,contrast_var] <- as.character(meta_cell[,contrast_var])
        # control_grp <- as.character(control_grp)
        # experimental_grp <- as.character(experimental_grp)
        print(paste("new class is", class(meta_cell[,contrast_var])))
      }

    # Create vectors of samples for each condition
    if (is.character(meta_cell[,contrast_var])) {
      ctrl_samples <- meta_cell[,donor_col][meta_cell[,contrast_var] == control_grp]
      exp_samples <- meta_cell[,donor_col][meta_cell[,contrast_var] == experimental_grp]

      if (is.null(ctrl_samples)) {
        print("not enough ctrl samples")
        { next }
      }
      if (is.null(exp_samples)) {
        print("not enough exp samples")
        { next }
      }

      # Subset meta_cell for the current contrast
      meta_cell.use <- meta_cell[meta_cell[,donor_col] %in% c(ctrl_samples, exp_samples),]
    }
    else {
      meta_cell.use <- meta_cell[!is.na(meta_cell[,contrast_var]),]
    }

    #scale continuous var cols if we want
    if (scale == TRUE) {
      num_cols <- unlist(lapply(meta_cell.use[,covariates_list], is.numeric), use.names = TRUE)
      num_cols <- names(num_cols)[which(num_cols)]

      meta_cell.use <- meta_cell.use %>% mutate(across(all_of(num_cols), ~scale(.)))

    }

    rownames(meta_cell.use) <- meta_cell.use[,donor_col]

    tryCatch({
      # Subset raw_counts for the donors in meta_cell
      raw_counts_subset <- raw_counts[, rownames(meta_cell.use)]
                                  }, error = function(e) {
              message("Error in reading counts for -> ", cell.use, " ", contrast_id, ": ", e$message)
          })        
      
      # Ensure that the column names of raw_counts are in the same order as the library identifiers in meta_cell
      message("Donors detectected in this contrast: ", nrow(meta_cell.use))

      #first check if enough cells in total
      cells_this_type <- all_n_cells[cell.use]

      if (cells_this_type < 100) {
        print(paste("cell type", cell.use, "doesn't have enough total cells"))

        {next}
      }

      keep_donors <- list()
      for (d in meta_cell.use[,donor_col]) {
        n_cells_donor <- meta_cell.use[which(meta_cell.use[,donor_col] == d), cell.use]

        if (n_cells_donor >= 20) {
          keep_donors <- keep_donors %>% append(d)
        }
      }

      if (length(keep_donors) < 2) {
        print(paste0("n donors left is ", length(keep_donors), " which is not enough..."))
        print(paste0("not enough donors left for cell type ", cell.use))
        { next }
      }

      print(paste0("Number of donors to keep based on minimum number of cells: ", length(keep_donors)))
      # message("Here are the donors: ", keep_donors)

      raw_counts_subset <- raw_counts_subset[,unlist(keep_donors)]

      meta_cell.use <- meta_cell.use[meta_cell.use[,donor_col] %in% unlist(keep_donors),]

      print("here are the unique contrast levels")
      print(unique(meta_cell.use[,contrast_var]))
      print(meta_cell.use[,contrast_var])


      #then check if enough cells in donor
      
      # Check if enough samples are present for both conditions
      if (length(unique(meta_cell.use[,contrast_var])) >= 2) {
          # tryCatch({
              
              ####====Apply filters====####
              ## Rebecca's new filters
              # Filter genes being tested: Ensure that 50% of donors per condition have a count > 5
              if (is.character(meta_cell.use[,contrast_var])){
                  ctrl_counts <- raw_counts_subset[,colnames(raw_counts_subset) %in% ctrl_samples]
                  exp_counts <- raw_counts_subset[,colnames(raw_counts_subset) %in% exp_samples]


                  if (!is.data.frame(ctrl_counts)) {
                    print("not enough samples in control group")
                    { next }
                  }
                  if (!is.data.frame(exp_counts)) {
                    print("not enough samples in experimental group")
                    { next }
                  }

                  ctrl_cutoff <- floor(length(ctrl_samples)/2)
                  exp_cutoff <- floor(length(exp_samples)/2)

                  ctrl_genes_to_keep <- list()
                  # loop through each row in the data frame
                  for (i in 1:nrow(ctrl_counts)) {
                  
                    # check if there are at least n_ND values greater than 5 in the current row
                    if (sum(ctrl_counts[i, ] >= 5) >= ctrl_cutoff) {
                        
                        # add the row name to the result list
                        ctrl_genes_to_keep <- append(ctrl_genes_to_keep, rownames(ctrl_counts[i, ]))
                    }
                  }
                  
                  exp_genes_to_keep <- list()
                  # loop through each row in the data frame
                  for (i in 1:nrow(exp_counts)) {
                  
                    # check if there are at least n_ND values greater than 5 in the current row
                    if (sum(exp_counts[i, ] >= 5) >= exp_cutoff) {
                        
                        # add the row name to the result list
                        exp_genes_to_keep <- append(exp_genes_to_keep, rownames(exp_counts[i, ]))
                    }
                  }
                  # genes_to_keep <- intersect(ctrl_genes_to_keep, exp_genes_to_keep)
                  genes_to_keep <- unique(c(ctrl_genes_to_keep, exp_genes_to_keep))
                  print("here are the first genes to keep")
                  print(genes_to_keep[1:50])

                  raw_counts_filter <- raw_counts_subset[rownames(raw_counts_subset) %in% genes_to_keep,]
              } 
              if (is.numeric(meta_cell.use[,contrast_var])) { #if contrast_var is numeric, we need to keep all samps to subset genes
                all_samps_cutoff <- floor(nrow(meta_cell.use)/2)

                genes_to_keep <- list()
                for (i in 1:nrow(raw_counts_subset)) {
                  if (sum(raw_counts_subset[i, ] >= 5) >= all_samps_cutoff) {
                    genes_to_keep <- append(genes_to_keep, rownames(raw_counts_subset[i, ]))
                  }
                }
                raw_counts_filter <- raw_counts_subset[rownames(raw_counts_subset) %in% genes_to_keep,]
              }

              if (length(genes_to_keep) < 10) {
                print("not enough genes with enough counts")
                { next }
              }

              #now re-define max.k based on number of donors
              max.k <- min(max.k.global, ncol(raw_counts_filter)-3)
              latent_vars <- paste0("W_", 1:max.k)

              counts_filtered <- as.matrix(raw_counts_filter)

              #open PDF to save all plots
              pdf(paste0(outdir, Analysis, "_all_plots.pdf"))

              ####====FIRST RUN WITH BASE MODEL====####
              # Run DESeq2
              if (is.character(meta[,contrast_var])) {
                  meta_cell.use[,contrast_var] <- factor(meta_cell.use[,contrast_var], levels = c(control_grp, experimental_grp)) #first one will be used as reference

                  # #lets get a filter to say there need to be 3 samples in each group
                  # dots <- as.symbol(contrast_var)
                  # count_samps <- meta_cell.use %>% group_by(.dots = dots) %>% summarize(n=n())
                  # print(count_samps)

              }

              skip <- FALSE
              tryCatch({
                dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                                            colData = meta_cell.use, 
                                            design = as.formula(deseq_formula))
                dds <- DESeq(dds, quiet = TRUE)

              },
              error = function(e) {
                print(paste0("Initial DESeq failed with message: ", e$message))
                skip <<- TRUE
              })
              if (skip) {
                { next }
              }
              
              
              res <- results(dds)
              res <- na.omit(res[order(res$padj), ])
              res <- res %>% as.data.frame()
              res <- cbind(gene = rownames(res), res)

              de_features_PVAL <- nrow(res[res$pvalue < 0.05,])
              de_features_FDR <- nrow(res[res$padj < 0.05,])

              p <- EnhancedVolcano(res,
                           lab = rownames(res),
                           x = "log2FoldChange",
                           y = "padj",
                           pCutoff = 0.05,
                           FCcutoff = 0,
                           cutoffLineType = 'blank',
                           xlim = c(min(res$log2FoldChange)-0.1, max(res$log2FoldChange)+0.1),
                           title = paste0("Base DESeq ", Analysis),
                            subtitle = paste0("n_DEGs=", de_features_FDR))
              print(p)

              # ggsave(paste0(outdir, Analysis, ".DESeq_base_volcano.pdf"), p, bg = "white")

              #save base DESeq results
              # Write DESeq results to file
              write.table(res, 
                      file = paste0(outdir,
                                    Analysis, ".DESeq_base.dds.res.tsv"), 
                      row.names = FALSE,
                      sep = '\t', quote = FALSE)

              
              ####====RUV-SEQ on base model ====####                
              #get BG features
              bg_features <- rownames(res)[res$pvalue > 0.5]
              bg_features <- intersect(bg_features, rownames(counts_filtered))
              message("bg features: ", length(bg_features))
              message("baseline de-genes: ", de_features_FDR)
              
              message("  - Running latent variables against covariants")
              
              # Create ruv-seq compatible object
              set_0 <- newSeqExpressionSet(counts_filtered,
                                          phenoData = meta_cell.use)
              set_0 <- betweenLaneNormalization(set_0, which="upper")
              set_max <- RUVg(set_0, bg_features, k = max.k)

              ## make corr mat with all known covariates and latent vars
              #set up matrices to hold stats values and p-values
              t_mat <- matrix(nrow = length(correlate_vars), ncol = length(latent_vars))
              rownames(t_mat) <- correlate_vars
              colnames(t_mat) <- latent_vars
          
              p_mat <- matrix(nrow = length(correlate_vars), ncol = length(latent_vars))
              rownames(p_mat) <- correlate_vars
              colnames(p_mat) <- latent_vars
          
              for (v in correlate_vars) {
                  print(v)
                  for (x in latent_vars) {
                      this_form <- as.formula(paste0(x, " ~ ", v))
                      stats_out <- lm(this_form, data = pData(set_max), na.action=na.omit)
                      # stats_out <- summary(stats_out)
                      stats_p <- summary(stats_out)$coefficients[2,4]
                      stats_t <- summary(stats_out)$coefficients[2,3]
          
                      t_mat[v, x] <- stats_t
                      p_mat[v, x] <- stats_p
                      
                      }
              }

              #plot correlation matrix
              # Get p-value matrix - bonferroni corrected
              p.df <- as.data.frame(p_mat) #*nrow(p_mat)*ncol(p_mat)
              # Function to get asteriks
              labs.function <- function(x){
                case_when(x >= 0.05 ~ "",
                          x < 0.05 & x >= 0.01 ~ "*",
                          x < 0.01 & x >= 0.001 ~ "**",
                          x < 0.001 ~ "***")
              }
              
              # Get asteriks matrix based on p-values
              p.labs <- p.df %>% mutate_all(labs.function)
              
              # Reshaping asteriks matrix to match ggcorrplot data output
              p.labs$Var1 <- as.factor(rownames(p.labs))
              p.labs <- reshape2::melt(p.labs, id.vars = "Var1", variable.name = "Var2", value.name = "lab")
              
              # Initial ggcorrplot
              cor.plot <- ggcorrplot(t_mat, hc.order = FALSE,
                                    lab = TRUE,
                                    title = Analysis, 
                                    lab_size = 2.25) + #original plots used 2.75
                        scale_fill_gradient2(low = "#3A86FF", 
                                            mid = "white", 
                                            high = "#FD2244",
                                            limits = c(-15, 15),
                                            midpoint = 0,
                                            oob = scales::squish) +
                        labs(fill = "t-statistic",
                            caption = "Bonferroni corrected p-values: * p < 0.05, ** p < 0.01, *** p < 0.001")
              p.labs$in.df <- ifelse(is.na(match(paste0(p.labs$Var1, p.labs$Var2),
                                                    paste0(cor.plot[["data"]]$Var1, cor.plot[["data"]]$Var2))),
                                        "No", "Yes")
              
              p.labs <- select(filter(p.labs, in.df == "Yes"), -in.df)
              
              # Add asteriks to ggcorrplot
              cor.plot.labs <- cor.plot +
                geom_text(aes(x = p.labs$Var1,
                              y = p.labs$Var2),
                          label = p.labs$lab,
                          nudge_y = 0.25,
                          size = 5) +
                guides(fill = guide_legend(title = "t-statistic"))
              
              print(cor.plot.labs)

              # ggsave(filename=paste0(outdir, 
              #                         Analysis, 
              #                         "_covariates_latent_var_correlation_mat.pdf"),
              #        plot = cor.plot.labs, device = "pdf", bg = "white")

              ####==== Find disease associated Ws ====#### 

              lm_results <- data.frame()
              for (x in latent_vars) {
                  this_form <- as.formula(paste0(x, " ~ ", contrast_var))
                  stats_out <- lm(this_form, data = pData(set_max), na.action=na.omit)
                  stats_p <- summary(stats_out)$coefficients[2,4]
                  stats_t <- summary(stats_out)$coefficients[2,3]

                  temp_df <- data.frame(W = x,
                                        T_stat = stats_t,
                                        p_val = stats_p)

                  lm_results <- rbind(lm_results, temp_df)
                          
              }

              lm_results$p_adj <- p.adjust(lm_results$p_val, method = "bonferroni")

              # Store it to add it in
              lm_results_sig = filter(lm_results, p_adj < 0.05) #filter on nominal or adjusted?
              disease_associated_Ws = unique(lm_results_sig$W)
              message("  - association found with disease status: \n    ", paste0(disease_associated_Ws,
                                                        collapse = "|"))
              
              ####==== Run Ruvseq on each K and run ANOVA ====#### 
              anova_res_all <- data.frame()
              for (K.i in 1:max.k) {
                  message("  Running RUVg: ", K.i, "/", max.k)

                  # Set up SeqExpressionSet
                  set <- RUVg(set_0, cIdx = bg_features, k = K.i)

                  # Store ANOVA results
                  #Get RLE manually
                  x1 <- normCounts(set)
                  y1 <- log(x1+1)
                  medi <- apply(y1, 1, median)
                  rle_mat <- apply(y1, 2, function(x) x - medi)

                  rle_long <- t(rle_mat) %>% as.data.frame()
                  gene_names <- colnames(rle_long)
                  rle_long[,donor_col] <- rownames(rle_long)

                  rle_long <- rle_long %>% merge(pData(set), by = donor_col)
                  rle_long <- rle_long %>% pivot_longer(cols = gene_names)

                  # ANOVA: RLE ~ sample
                  aov_form <- as.formula(paste("value ~", donor_col))
                  fit <- aov(aov_form, data = rle_long)
                  anova_summary <- summary(fit)[[1]]
                  f_stat <- anova_summary[donor_col, "F value"]
                  res_var <- anova_summary["Residuals", "Mean Sq"]

                  # Store ANOVA results
                  Anove_stat_tmp <- data.frame(
                                              celltype = cell.use,
                                              contrast = contrast_id,
                                              donors = nrow(pData(set)),
                                              k = paste0("k_", K.i),
                                              k_num = K.i,
                                              f_statistic = f_stat,
                                              residual_variance = res_var
                  )
                  anova_res_all <- rbind(anova_res_all, Anove_stat_tmp)
              }

              ####==== Find Best K ====#### 

              df <- anova_res_all[order(anova_res_all$k_num),]
              print(head(df))

              # Find local minima
              mins_idx <- which(df$f_statistic < lag(df$f_statistic) & df$f_statistic < lead(df$f_statistic))
              min_k     <- df$k_num[mins_idx]
              min_f     <- df$f_statistic[mins_idx]

              # If just one local minimum, take it immediately
              if (length(mins_idx) == 1) {
                best_k <- min_k
              } 
              else if (length(mins_idx) > 1) {
                # Compute ≥10-point drops between successive minima
                drops <- c(FALSE, (min_f[-1] <= min_f[-length(min_f)] - f_tresh))

                # Build candidates and filter by drop
                candidates <- data.frame(k = min_k,
                                        f = min_f,
                                        drop = drops) %>% filter(drop)
                print("here is candidates")
                print(head(candidates))

                if (nrow(candidates) > 0) {
                  # If some meet the drop criterion, pick the one with smallest f
                  best_k <- candidates$k[which.min(candidates$f)]
                  # best_k <- candidates %>% slice_min(order_by = f, n = 1) %>% pull(k)
                } 
                else {
                  # Otherwise, fall back to the single deepest local minimum
                  #     (i.e. the one with the lowest f among all minima)
                  best_k <- df$k_num[which.min(df$f_statistic)]

                } 
              } 
              else {
                # No local minima at all?  Fall back to the global minimum of f
                best_k <- df$k_num[which.min(df$f_statistic)]
              }

              message("  - Found BestK -> ", best_k)                    

              #Plot RLE plot

              y_base <- round(anova_res_all$f_statistic[which(anova_res_all$k == "k_1")])
              y_min <- round(min(anova_res_all$f_statistic, na.rm = TRUE))

              # Plot with raw residual variance
              k_levels <- paste0("k_", 1:max.k)
              anova_res_all$k <- factor(anova_res_all$k, levels = k_levels)
              gg <- ggplot(anova_res_all, aes(x = k, y = f_statistic, group = 1)) +
                geom_line() +

                geom_point(
                  aes(size = residual_variance, fill = residual_variance),
                  shape = 21,
                  color = "black"
                ) +

                scale_fill_viridis(
                  name = "Residual variance",
                  option = "D"
                ) +
                scale_size_continuous(
                  name = "Residual variance"
                ) +
                scale_y_continuous(
                  breaks = c(y_base, y_min)
                ) +

                theme_classic() +
                labs(
                  x     = "Number of RUV factors (k)",
                  y     = "\n F-statistic\n(ANOVA: RLE ~ sample)",
                  title = paste0(Analysis, "\nBest-k -> ", best_k)
                ) +
                theme(
                  axis.title    = element_text(size = 14),
                  axis.text     = element_text(size = 12),
                  axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1)
                )

              print(gg)

              # ggsave(filename = paste0(outdir, 
              #                           Analysis, 
              #                           "_Anova_F-stat.pdf"),
              #                             plot = gg, device = "pdf", bg = "white")
              
              
              ####==== Generate best formula ====####                 
              all_ws <- 1:best_k            
              all_ws <- paste0("W_", all_ws)

              W_terms = all_ws[!all_ws %in% disease_associated_Ws]

              best_formula <- paste0("~", paste0(W_terms, collapse = "+"), "+", contrast_var)
              message("  - Best formula ->  ", best_formula)
              
              
              ####==== Final Deseq Run ====####  
              # Run Deseq again - to test number of de-features and more
              if (is.character(meta[,contrast_var])) {
                  pData(set_max)[,contrast_var] <- factor(pData(set_max)[,contrast_var], levels = c(control_grp, experimental_grp))
              }
              
              dds <- DESeqDataSetFromMatrix(countData = counts(set_max),
                                            colData = pData(set_max),
                                            design = as.formula(best_formula))
              dds <- DESeq(dds, quiet = TRUE)
              res <- results(dds)
              res <- na.omit(res[order(res$padj), ])
              res <- res %>% as.data.frame()
              res <- cbind(gene = rownames(res), res)
              print(head(res))

              RUV_de_features <- nrow(res[res$pvalue < 0.05,])
              RUV_fdr_features <- nrow(res[res$padj < 0.05,])

              p <- EnhancedVolcano(res,
                           lab = rownames(res),
                           x = "log2FoldChange",
                           y = "padj",
                           pCutoff = 0.05,
                           FCcutoff = 0,
                           cutoffLineType = 'blank',
                           xlim = c(min(res$log2FoldChange)-0.1, max(res$log2FoldChange)+0.1),
                           title = paste0("Best RUVSeq ", Analysis),
                            subtitle = paste0("n_DEGs=", RUV_fdr_features))

              print(p)

              dev.off()

              # ggsave(paste0(outdir, Analysis, ".RuVseq.best_volcano.pdf"), p, bg = "white")

              message("  - de-features: ", de_features_FDR)

              # Write DESeq results to file
              write.table(res, 
                      file = paste0(outdir,
                                    Analysis, ".RuVseq.dds.res.tsv"), 
                      row.names = FALSE,
                      sep = '\t', quote = FALSE)
              
              
              # Gather stats
              deseq.stats_tmp <- data.frame(
                  celltype = cell.use,
                  contrast = contrast_id,
                  Donors = nrow(meta_cell.use),
                  initial_formula = deseq_formula,
                  RUV_formula = best_formula,
                  initial_de_features = de_features_PVAL,
                  initial_fdr_features = de_features_FDR,
                  RUV_de_features = RUV_de_features,
                  RUV_fdr_features = RUV_fdr_features)

              deseq_stats <- rbind(deseq_stats, deseq.stats_tmp)

              write.table(deseq_stats, 
                      file = paste0(outdir,
                                    Analysis, ".RuVseq.deseq.stats.tsv"), 
                                    row.names = FALSE,                      
                                    sep = '\t', quote = FALSE)
          
  } 
}

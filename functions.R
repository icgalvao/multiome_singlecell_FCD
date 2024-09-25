library(RColorBrewer)

### Functions


# Visualize colors 
seePalette <- function(){
  # create data that references the palette
  colorKey = data.frame(colorName=names(colors))
  # plot with ggplot, referencing the palette
  ggplot(data=colorKey, aes(x=1, y = 1:nrow(colorKey), fill=colorName, label=colorName)) +
    geom_tile() +
    scale_fill_manual(values = colors) +
    theme_void()+
    theme(legend.position="none") + 
    geom_text()
}

Viz <- function(obj, red){
  
  print(DimPlot(obj,
                reduction = red, label = T))
  
  
  # Verify batch effect and other effects 
  
  # By donor
  print(DimPlot(obj,
                reduction = red, 
                group.by = "patient"))
  
  # By sample 
  print(DimPlot(obj,
                reduction = red, 
                group.by = "sample"))
  
  
  # By sequencing batch
  print(DimPlot(obj,
                reduction = red, 
                group.by = "seq_batch"))
  
  # By histological assessment 
  print(DimPlot(obj,
                reduction = red,
                group.by = "histological_assessment"))
  
  # By lobe
  print(DimPlot(obj,
                reduction = red, 
                group.by = "lobe"))
  
  # By sex
  print(DimPlot(obj,
                reduction = red, 
                group.by = "sex"))
  
  # By age
  print(DimPlot(obj,
                reduction = red, 
                group.by = "age"))
  
}



# Wrapper function to process GEX until determining dataset dimensionality 
# determining the dimensionality must be done separately for each sample 
processGEX1 <- function(obj, vars_to_regress = "percent.mt"){
  
  DefaultAssay(obj) <- "RNA"
  
  # Normalize data 
  obj <- NormalizeData(obj, normalization.method = "LogNormalize") 
  
  # Find variable features 
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  
  # Scale data 
  all.genes <- rownames(obj)
  obj <- ScaleData(obj, features = all.genes, vars.to.regress = vars_to_regress)
  
  # Run PCA - linear dimensional reduction
  obj <- RunPCA(obj, 
                npcs = 50, 
                features = VariableFeatures(obj),
                verbose = FALSE)
  
  return(obj)
  
}




# When using future parallelization, set the strategy "multiseession"  if using 
# RStudio or "multicore" if using R 



# Wrapper function to process GEX until determining dataset dimensionality 
# determining the dimensionality must be done separately for each sample 
processGEX1_RStudio <- function(obj, vars_to_regress = "percent.mt"){
  
  DefaultAssay(obj) <- "RNA"
  
  
  # Alow paralellization 
  plan("multisession", workers = 16)
  
  # Normalize data 
  obj <- NormalizeData(obj, normalization.method = "LogNormalize") 
  
  
  # Convert back to sequential 
  plan("sequential")
  
  # Find variable features 
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  
  
  # Enable parallelization 
  plan("multisession", workers = 16)
  
  # Scale data 
  all.genes <- rownames(obj)
  obj <- ScaleData(obj, features = all.genes, vars.to.regress = vars_to_regress)
  
  # convert back to sequential 
  plan("sequential")
  
  # Run PCA - linear dimensional reduction
  obj <- RunPCA(obj, 
                npcs = 50, 
                features = VariableFeatures(obj),
                verbose = FALSE)
  
  return(obj)
  
}



# Wrapper function to process GEX until determining dataset dimensionality 
# determining the dimensionality must be done separately for each sample 
processGEX1_R <- function(obj, vars_to_regress = "percent.mt"){
  
  DefaultAssay(obj) <- "RNA"
  
  
  # Alow paralellization 
  plan("multicore", workers = 16)
  
  # Normalize data 
  obj <- NormalizeData(obj, normalization.method = "LogNormalize") 
  
  
  # Convert back to sequential 
  plan("sequential")
  
  # Find variable features 
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  
  
  # Enable parallelization 
  plan("multicore", workers = 16)
  
  # Scale data 
  all.genes <- rownames(obj)
  obj <- ScaleData(obj, features = all.genes, vars.to.regress = vars_to_regress)
  
  # convert back to sequential 
  plan("sequential")
  
  # Run PCA - linear dimensional reduction
  obj <- RunPCA(obj, 
                npcs = 50, 
                features = VariableFeatures(obj),
                verbose = FALSE)
  
  return(obj)
  
}




# Wrapper function to process GEX using SCTtransform() until determining dataset 
#dimensionality 
# determining the dimensionality must be done separately for each sample 
processGEX1_SCT <- function(obj, vars_to_regress){
  
  DefaultAssay(obj) <- "RNA"
  
  
  # Process data with SCT 
  obj <- SCTransform(obj,
                     assay = "RNA",
                     new.assay.name = "SCT",
                     method = "glmGamPoi", 
                     vars.to.regress = vars_to_regress,
                     vst.flavor = "v2",
                     verbose = FALSE)
  
  
  # Run PCA - linear dimensional reduction
  obj <- RunPCA(obj, 
                npcs = 50, 
                features = VariableFeatures(obj),
                verbose = FALSE)
  
  return(obj)
  
}



# Wrapper function to process GEX after determining the dimensionality of the dataset 

processGEX2 <- function(obj, ndims = 35, res){
  
  # Run UMAP - non-linear dimensional reduction 
  obj <- RunUMAP(obj,
                 dims = 1:ndims,
                 reduction = "pca",
                 reduction.name = "umap.rna")
  
  # Find Neighbors 
  obj <- FindNeighbors(obj, 
                       dims = 1:ndims)
  
  # Find clusters 
  obj <- FindClusters(obj,
                      resolution = res, 
                      algorithm = 3)
  
  return(obj)
  
}





# Wrapper function to process scATAC-seq data 

processATAC <- function(obj, assay = "peaksMACS2", res = 0.5){
  
  # Set default assay to the one with new peaks 
  DefaultAssay(obj) <- assay
  
  # Process data through latent semantic indexing 
  obj <- FindTopFeatures(obj, min.cutoff = 5) #feature selection, for dimensional reduction (peaks in > 5 cells)
  obj <- RunTFIDF(obj)   #TF-IDF normalization in the matrix
  obj <- RunSVD(obj)     # linear dimensional reduction (analogous to PCA)
  
  
  # Run UMAP - non-linear dimensional reduction 
  obj <- RunUMAP(object = obj, reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac")
  
  # Find neighbors 
  obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = 2:30)
  
  # Find Clusters
  obj <- FindClusters(object = obj, resolution = res, algorithm = 3)
  
  
  return(obj)
  
}



# Add pretty labels 
labels <- function(plot, id){
  
  p <- LabelClusters(plot,
                      id,
                      repel = TRUE, 
                      box = TRUE,
                      size = 3,
                      fill = "white",  
                      segment.color = 'grey50',
                      fontface = 'bold',
                      alpha = 0.5,
                      segment.alpha = 0.8,
                      label.size = NA,
                      force = 2,
                      segment.size = 0.5,
                      arrow = arrow(length = unit(0.01, 'npc'))
  )


  return(p)

}


# Print paged dataframe on rmarkdown 
print_paged_df <- function(...) {
  cat(rmarkdown:::print.paged_df(rmarkdown::paged_table(...)))
}




# The following function was taken from the escape package github 
# As the new version of escape does not work with the R version we are using 
# It is used to perform statistical tests on metadata columns 
# To use it, you must pass an object with columns for the enrichment scores 
# and the column for the variable you want to group by as the enriched parameter

# Functions for escape v1.12 

#' Perform significance testing between groups and enrichement scores.
#' 
#' This functions takes the enrichment scores and performs statistical 
#' testing to evaluate the difference by group selected. The function 
#' can perform 5 tests: 1) Welch's T test (T.test), 2) Logistic 
#' Regression (LR), 3) Wilcoxon Rank Sum Test (Wilcoxon), 
#' 4) one-way ANOVA (ANOVA), and 5) Kruskal-Wallis (KW). The latter 
#' two output will include the individual comparisons between groups 
#' using TukeyHSD for ANOVA and pairwise Wilcoxon Rank Sum Test 
#' for KW. The output includes adjusted p-values based on the 
#' Benjamini Hochberg method. 
#' 
#'
#' @param enriched The output of \code{\link{enrichIt}}.
#' @param group The parameter to group for the comparison, should a column of 
#' the enriched input
#' @param gene.sets Names of gene sets to compare
#' @param fit The test used for significance, 2 group: Wilcoxon, LR, T.test.
#' Multigroup: ANOVA or KW.
#' @importFrom dplyr select_if
#' @importFrom broom tidy
#' @importFrom reshape2 melt
#' @importFrom stats TukeyHSD median glm wilcox.test pairwise.wilcox.test kruskal.test
#' 
#' @examples 
#' ES2 <- readRDS(url(
#' "https://ncborcherding.github.io/vignettes/escape_enrichment_results.rds"))
#' output <- getSignificance(ES2, group = "Type", fit = "T.test")
#' 
#' @export
#'
#' @seealso \code{\link{enrichIt}} for generating enrichment scores.
#' @return Data frame of test statistics
getSignificance <- function(enriched, 
                            group = NULL,
                            gene.sets = NULL,
                            fit = NULL) {
  fit <- match.arg(fit,  choices = c("T.test", "ANOVA", "Wilcoxon", "LR", "KW"))
  group2 <- enriched[,group]
  gr_names <- unique(group2)
  if (!is.null(gene.sets)) {
    input <- enriched[,colnames(enriched) %in% gene.sets]
  } else {
    input <- select_if(enriched, is.numeric)
  }
  medians <- get.medians(input, group2)
  output <- NULL
  if (fit == "T.test" || fit == "Wilcoxon" || fit == "LR") {
    if (length(unique(group2)) != 2) {
      message("Ensure the group selection has only two levels for T.test 
                fit") 
    } else {
      if (fit == "T.test") {
        out <- lapply(input, function(x) t.test(x ~ group2))
        stat <- "T"
      } else if (fit == "Wilcoxon") {
        out <- lapply(input, function(x) wilcox.test(x ~ group2))
        stat <- "W"
      }  else if (fit == "LR") {
        group2 <- ifelse(group2 == gr_names[1], 0,1)
        out <- lapply(input, function(x) glm(group3 ~ x, family = "binomial"))
        out <- lapply(out, function(x) tidy(x)[2,])
        stat <- "L"
      }
      for (i in seq_along(out)) {
        df <- out[[i]]
        mat <- c(df$statistic, df$p.value)
        output <- rbind(output,mat)
      }
      output <- as.data.frame(output)
      colnames(output) <- c(paste0(stat, ".statistic"), "p.value")
    }
  } else if (fit == "ANOVA") {
    if (length(unique(group2)) <= 2) {
      message("Ensure the group selection has more than two levels 
                for ANOVA fit") }
    out <- lapply(input, function(x) aov(x ~ group2))
    for (i in seq_along(out)) {
      df <- out[[i]]
      tukey <- TukeyHSD(df)
      ind.p.values <- melt(tukey$group2[,4])
      names.ind.p.values <- gsub("-", "v", rownames(ind.p.values))
      names.ind.p.values <- paste0(names.ind.p.values,".p.value")
      fval <- summary(df)[[1]]$'F value'[[1]]
      pval <- summary(df)[[1]]$'Pr(>F)'[[1]]
      output <- rbind(output, c(fval, pval, t(ind.p.values)))
    }
    output <- as.data.frame(output)
    colnames(output) <- c("f.value", "p.value", names.ind.p.values)
  } else if (fit == "KW") {
    if (length(unique(group2)) <= 2) {
      message("Ensure the group selection has more than two levels 
                for Kruskal-Wallis test")}
    out <- lapply(input, function(x) kruskal.test(x ~ group2))
    out.ind <- lapply(input, function(x) pairwise.wilcox.test(x, group2, p.adjust.method = "BH"))
    for (i in seq_along(out)) {
      ind.p.values <- na.omit(melt(out.ind[[i]]$p.value))
      names.ind.p.values <- paste0(ind.p.values$Var1, "v", ind.p.values$Var2)
      names.ind.p.values <- paste0(names.ind.p.values,".p.value")
      ind.p.values <- ind.p.values[,3]
      Chi.squared <- out[[i]]$statistic 
      pval <- out[[i]]$p.value
      output <- rbind(output, c(Chi.squared, pval, t(ind.p.values)))
    }
    output <- as.data.frame(output)
    colnames(output) <- c("Chi.square", "p.value", names.ind.p.values)
  }
  rownames(output) <- colnames(input)
  output$FDR <- p.adjust(output$p.value) 
  output <- cbind.data.frame(output, medians)
  return(output)
}

get.medians<- function(input, group2) {
  input <- cbind.data.frame(group2, input)
  num <- ncol(input)-1
  med <- input %>%
    group_by(group2) %>%
    summarise(across(seq_len(all_of(num)), median))
  new.grouping <- unlist(med[,1])
  med <- as.data.frame(med[,seq_len(num) + 1])
  rownames(med) <- paste0("median.", new.grouping)
  med <- as.data.frame(t(med))
  return(med)
}


# Violin plots

# Define a function to plot violin plots 


asterisks <- function(x){
  
  if (x < 0.0001){
    test <- "****"
  } else if (x < 0.001){
    test <- "***"
  } else if (x < 0.01){
    test <- "**"
  } else if (x < 0.05){
    test <- "*"
  } else {
    test <- ""
  }
  
  return(test)  
  
}



# The function below plots custom violin plots comparing FCD IIa, IIb, and Normal
# This function automatically adds the asterisks for significance when 
#all three groups are in the image 

# Parameters:
# @ so: seurat object 
# @ variable: the column in metadata to plot
# @ title: the plot title, defaults to the variable name 
# @ res: the table with the statistical tests results, must contain an FDR column
# and pairwise wilcox tests - the table generated from getSignificance works perfectly


plot_custom_ha <- function(so, variable, title = variable, res){
  
  
  res <- test_res %>% filter(gene.set == variable)
  
  Switch <- res$FDR < 0.1
  
  dittoPlot(so, 
            var = variable,
            group.by = "ha3", 
            plots = c("jitter", "vlnplot", "boxplot"), 
            jitter.size = 0.1,
            ylab = "Enrichment Score",
            main = title,
            theme = theme_classic() + theme(plot.title = element_text(size = 10, hjust = 0.5))) + 
    scale_fill_manual(values = ha_cols) + 
    labs(fill = "Histological\nassessment") +
    scale_y_continuous(expand = expansion(add = c(0.05, 0.3*max(so@meta.data[, variable])))) + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) + 
    annotate(geom = "segment", x = 1.05, xend = 1.95, y = max(so@meta.data[, variable]) + 0.05*max(so@meta.data[, variable]), yend = max(so@meta.data[, variable]) + 0.05*max(so@meta.data[, variable]), linewidth = 1.3) +
    annotate(geom = "segment", x = 2.05, xend = 2.95, y = max(so@meta.data[, variable]) + 0.05*max(so@meta.data[, variable]), yend = max(so@meta.data[, variable]) + 0.05*max(so@meta.data[, variable]), linewidth = 1.3) +
    annotate(geom = "segment", x = 1.05, xend = 2.95, y = max(so@meta.data[, variable]) + 0.15*max(so@meta.data[, variable]), yend = max(so@meta.data[, variable]) + 0.15*max(so@meta.data[, variable]), linewidth = 1.3) +
    {if(Switch)annotate(geom = "text", x = 1.5, y = max(so@meta.data[, variable]) + 0.06*max(so@meta.data[, variable]), label = asterisks(res$`FCD IIbvFCD IIa.p.value`), size = 10)} + 
    {if(Switch)annotate(geom = "text", x = 2.5, y = max(so@meta.data[, variable]) + 0.06*max(so@meta.data[, variable]), label = asterisks(res$`NormalvFCD IIb.p.value`), size = 10)} + 
    {if(Switch)annotate(geom = "text", x = 2, y = max(so@meta.data[, variable]) + 0.16*max(so@meta.data[, variable]), label = asterisks(res$`NormalvFCD IIa.p.value`), size = 10)}
  
  
}


# Set1 palette for more than 8 colors

set_qual_pal <- function(n){ 
  
  library(RColorBrewer)
  
  getPalette = colorRampPalette(brewer.pal(8, "Set1"))
  colors = getPalette(n)
  
  return(colors)
  
}


# Define the colors to be  used in the final plots 
# Defining colors for the histological assessment
ha_cols <- brewer.pal(3, "Dark2")
names(ha_cols) <- c("FCD IIa", "FCD IIb", "Normal")

ha4_palette <- brewer.pal(4, "Dark2")
names(ha4_palette) <- c("FCD IIa", "FCD IIb", "Internal Controls", 
                        "Autopsy Controls")

ha4_abbr_palette <- brewer.pal(4, "Dark2")
names(ha4_abbr_palette) <- c("FCD IIa", "FCD IIb", "Int. Controls", 
                             "Aut. Controls")

# Palette for cell types 
pal <- c("#E41A1C", "#00b300", "#3FCA2B", "#148210", "#33B6EF", "#896191", 
         "#CB6651", "#FF980A", "#FFF22D", "#C9992C", "#653B04", "#F781BF", 
         "#0d9e33", "#33B6EF")
names(pal) <- c("Astro", "Neuron", "EN Upper", "EN Deep", "Endo/VLMC", "IN", 
               "Microglia", "Oligo", "OPC", "OPC Diff.", "Unclassified",
               "VLMC", "EN", "Endo")

# Palette for grey/white matter 
colors_matter <- c(GM="#737272", WM="#faeebe")



# Cell type markers

EN_markers <- c("KCNIP4", "R3HDM1", "PPP3CA", "SATB2", "SLC17A7")
uEN_markers <- c("CUX2", "CUX1", "CALB1")
dEN_markers <- c("SEMA3E", "FOXP2", "RORB", "NR4A2")
IN_markers <- c("GAD1", "GAD2")
Lamp5_markers <- c("LAMP5", "FGF13", "ADARB2")
Pvalb_markers <- c("PVALB", "BTBD11", "ZNF804A", "SCN1A-AS1")
Sst_markers <- c("SST", "GRIK1", "SYNPR")
Scng_markers <- c("ADARB2", "PRELID2", "CXCL14")
astro_markers <- c("SLC1A2", "SLC1A3", "GFAP", "AQP4")
micro_markers <- c("CD68", "PTPRC", "DOCK8", "CD86")
endo_markers <- c("CLDN5", "PECAM1", "CD34", "FLT1")
oligo_markers <- c("PLP1", "MOBP", "MBP", "ST18")
OPC_markers <- c("PDGFRA", "VCAN", "TNR", "GPR17")






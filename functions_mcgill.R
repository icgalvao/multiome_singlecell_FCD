# Script containing functions and examples of function usage
# For faceted UMAPs and bar plots within joint Seurat objects

library(tidyverse)
library(data.table)
library(glue)
library(DT)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggrepel)
library(Seurat)
library(Signac)
library(circlize)
library(reshape2)
library(scDblFinder)
library(gridExtra)

# ---- Custom functions for formatting ggplots ----

#' Remove the legend in a ggplot
#'
#' @return A theme element to hide legend
#' @export
#'
#' @author Selin Jessa
#'
#' @examples
#' # gg <- mpg %>%
#' # filter(class %in% c("compact", "suv")) %>%
#' # ggplot(aes(x = displ, y = hwy)) +
#' # geom_point(aes(color = factor(year))) +
#' # facet_wrap(~ class, ncol = 2)
#'
#' # gg
#' # gg + noLegend()
hide_legend <- function() {
  
  theme(legend.position = "none")
  
}

#' Remove axis ticks and tick labels from a ggplot
#'
#' @return A theme element to remove ticks
#' @export
#'
#' @author Selin Jessa
hide_ticks <- function() {
  
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())
  
}

#' Apply a clean theme to a ggplot2 object
#'
#' @references https://github.com/sjessa/ggmin
#'
#' @importFrom ggplot2 theme_light theme
#' @author Selin Jessa
#' @export
theme_min <- function(base_size = 11, base_family = "",
                      border_color = "grey90",
                      border_size = 1) {
  
  theme_light(base_size = 11, base_family = "") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, color = border_color, size = border_size),
      axis.ticks = element_line(color = border_color),
      strip.background = element_rect(fill = NA, color = NA),
      strip.text.x = element_text(color = "black", size = rel(1.2)),
      strip.text.y = element_text(color = "black", size = rel(1.2)),
      title = element_text(size = rel(0.9)),
      axis.text = element_text(color = "black", size = rel(0.8)),
      axis.title = element_text(color = "black", size = rel(0.9)),
      legend.title = element_text(color = "black", size = rel(0.9)),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), color = "black"),
      legend.key = element_rect(color = NA, fill = NA),
      legend.background = element_rect(color = NA, fill = NA)
    )
}



# ---- Custom function to plot UMAP ----

# Plot using embeddings (UMAP coordinates) extracted from Seurat objects
# Create plots using ggplot for more customization
# Function modified from Samantha Worme by Rodrigo Lopez-Gutierrez

plot_umap <- function(embedding,
                      colour_by = NULL,
                      colours = NULL,
                      colour_by_type = "discrete",
                      label = FALSE,
                      point_size = 0.8,
                      alpha = 0.8,
                      legend = ifelse((is.null(colour_by)) && (label), FALSE, TRUE),
                      label_repel = TRUE,
                      label_size = 4,
                      cells = NULL,
                      order_by = NULL,
                      clusters_to_label = NULL,
                      hide_ticks = TRUE,
                      title = NULL,
                      label_short = FALSE,
                      na_color = "gray80",
                      limits = NULL,
                      hide_axes = FALSE) {
  
  if (!is.null(order_by)) {
    
    # Check the variable is usable for sorting
    if (!is.numeric(embedding[[order_by]]) && !is.factor(embedding[[order_by]])) {
      
      stop("The variable specified in 'order_by' is neither numeric ",
           "nor a factor. If the column is of type character, consider ",
           "converting it to a factor. Otherwise, pass the name of a numeric column.")
      
    }
    
  } else {
    
    order_by <- colour_by
    
  }
  
  if (!is.null(cells)) embedding <- embedding %>% filter(cell %in% cells)
  if (!is.null(order_by)) embedding <- embedding %>% arrange(desc(is.na(!!sym(order_by))), !!sym(order_by))
  
  gg <- ggplot(embedding, aes(x = UMAP_1, y = UMAP_2))
  
  if (is.null(limits)) lims <- c(NA, NA)
  else lims <- limits
  
  gg <- gg +
    geom_point(aes_string(colour = colour_by), size = point_size, alpha = alpha)
  
  if (!is.null(colours)) {
    
    if (colour_by_type == "discrete") gg <- gg + scale_color_manual(values = colours, na.value = na_color)
    
    else if (colour_by_type == "continuous") {
      
      gg <- gg + scale_color_gradientn(colours = colours,
                                       na.value = na_color,
                                       limits = lims)
    }
    
  } else {
    
    if (colour_by_type == "continuous") {
      
      gg <- gg + scale_color_gradientn(colours = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "OrRd"))(n = 100),
                                       na.value = na_color,
                                       limits = lims)
    }
  }
  
  # Label clusters
  if (label) {
    
    centers <- embedding %>%
      group_by(seurat_clusters) %>%
      summarise(mean_x = median(UMAP_1),
                mean_y = median(UMAP_2))
    gg <- gg + ggrepel::geom_label_repel(data = centers,
                                         aes(x = mean_x, y = mean_y),
                                         label = centers$seurat_clusters,
                                         size = label_size,
                                         segment.color = 'grey50',
                                         fontface = 'bold',
                                         alpha = 0.8,
                                         segment.alpha = 0.8,
                                         label.size = NA,
                                         force = 2,
                                         segment.size = 0.5,
                                         arrow = arrow(length = unit(0.01, 'npc')))
    
  }
  
  gg <- gg + theme_min()
  
  if (hide_ticks) gg <- gg + hide_ticks()
  
  # Set the right axes titles
  if (hide_axes) gg <- gg + xlab(NULL) + ylab(NULL)
  
  # More aesthetics
  if (!legend) gg <- gg + hide_legend()
  
  if (!is.null(title)) gg <- gg + ggtitle(title)
  
  return(gg)
  
}



# ---- Example of extracting embedding from Seurat object ----

# ** Example code, to be adapted for your Seurat object! **
# From Rodrigo Lopez-Gutierrez
# Here, the seurat object is stored in a variable "seurat"
# 
# seurat_embedding <- as.data.frame(seurat[["umap"]]@cell.embeddings[, c(1,2)]) %>% rownames_to_column(var = "Cell")
# seurat_embedding <- left_join(seurat@meta.data %>% rownames_to_column("Cell"), seurat_embedding, by = "Cell")
# 
# # ---- Function to plot faceted UMAPs ----
# 
# # Faceted UMAPs highlighting different sets of cells at a time
# # Function requires above function plot_umap
# # By Rodrigo Lopez-Gutierrez
# # Here, the metadata column used for faceting UMAPs is "orig.ident"
# 
# samples_grid_order <- unique(seurat@meta.data$orig.ident)
# 
# ggs <- lapply(samples_grid_order, function(n) { 
#   # output number of cells per cluster for in plot title
#   num_cells <- seurat_embedding %>% 
#     filter(orig.ident == n) %>% 
#     nrow(.)
#   
#   seurat_embedding$highlight <- ifelse(seurat_embedding$orig.ident == n, "yes", "no")
#   seurat_embedding$order <- ifelse(seurat_embedding$orig.ident == n, 2, 1)
#   
#   gg <- plot_umap(seurat_embedding, colour_by = "highlight", order_by = "order", colours = c(yes = "#313A88", no = "gray80"), point_size = 0.8, label = FALSE, legend = FALSE, hide_axes = TRUE) + ggtitle(paste0(n, " - ", num_cells, " cells"))
#   
#   return(gg)
# })



# ---- Function to plot bar plot proportions in Seurat objects ----

# Function to inspect proportions of specific factor at each cluster in the integrated Seurat objects
# By Rodrigo Lopez-Gutierrez

# @param --- df
# data frame corresponding to the seurat metadata

# @param --- col_palette
# color palette vector with colors corresponding to the selected lvl_analysis

# @param --- lvl_analysis
# Select analysis level, e.g.: orig.ident, scDblFinder.class
##### Level must be a column in data frame with discrete values
##### EXAMPLE: lvl_analysis = "orig.ident"

# @param --- x_axis_title
# Rename the title of the x-axis for the plot generated

clust_lvl_prop <- function(df, 
                           col_palette, 
                           lvl_analysis = "orig.ident", 
                           x_axis_title = "Proportion of cells from each sample"){
  
  # counts
  data <- df %>% 
    dplyr::group_by(seurat_clusters) %>% 
    dplyr::mutate(count_cluster = n()) %>% 
    dplyr::arrange(desc(!!sym(lvl_analysis))) %>%
    dplyr::group_by(!!sym(lvl_analysis), seurat_clusters) %>% 
    dplyr::mutate(count_level_cluster = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(seurat_clusters, !!sym(lvl_analysis),
           count_cluster, count_level_cluster) %>% 
    distinct()
  
  # proportions for readability
  data <- data %>% 
    dplyr::mutate(prop_level_cluster = count_level_cluster / count_cluster)
  
  # Now plot data for selected level proportions
  plot_data <- data %>% dplyr::select(seurat_clusters, count_cluster, prop_level_cluster, count_level_cluster, !!sym(lvl_analysis)) %>% distinct() 
  
  plot_labels <- plot_data %>% 
    dplyr::arrange(desc(!!sym(lvl_analysis))) %>%
    dplyr::group_by(seurat_clusters) %>% 
    dplyr::summarise(label_level = paste(count_level_cluster, collapse = ", "))
  
  prop_plot <- ggplot(plot_data, aes(x = prop_level_cluster, y = seurat_clusters)) + 
    geom_bar(aes(fill = !!sym(lvl_analysis)), stat = "identity") + 
    geom_text(data = plot_labels,
              aes(x  = 1.03, y = seurat_clusters, label = label_level), hjust = 0,
              size = 2.5) +
    coord_cartesian(clip = "off") + 
    scale_fill_manual(values = col_palette) +  
    theme(panel.background = element_blank(), 
          legend.position = "none") + 
    xlim(0, 1.5) +
    ylab("Cluster") +
    xlab(x_axis_title) +
    guides(fill = guide_legend(reverse = T))
  
  return(prop_plot)
  
}

# Functions obtained from fungenomics/cytobox - mcgill university Dr. Kleiman's 
#lab 

annot_lvl_prop <- function(df, 
                           col_palette, 
                           annot = "consensus.annot.simp",
                           lvl_analysis = "orig.ident", 
                           x_axis_title = "Proportion of cells from each sample"){
  
  # counts
  data <- df %>% 
    dplyr::group_by(!!sym(annot)) %>% 
    dplyr::mutate(count_cluster = n()) %>% 
    dplyr::arrange(desc(!!sym(lvl_analysis))) %>% 
    dplyr::group_by(!!sym(lvl_analysis), !!sym(annot)) %>% 
    dplyr::mutate(count_level_cluster = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(!!sym(annot), !!sym(lvl_analysis),
                  count_cluster, count_level_cluster) %>% 
    distinct()
  
  # proportions for readability
  data <- data %>% 
    dplyr::mutate(prop_level_cluster = count_level_cluster / count_cluster)
  
  # Now plot data for selected level proportions
  plot_data <- data %>% dplyr::select(!!sym(annot), count_cluster, prop_level_cluster, count_level_cluster, !!sym(lvl_analysis)) %>% distinct() 
  
  plot_labels <- plot_data %>% 
    dplyr::arrange(desc(!!sym(lvl_analysis))) %>% 
    dplyr::group_by(!!sym(annot)) %>% 
    dplyr::summarise(label_level = paste(count_level_cluster, collapse = ", "))
  
  prop_plot <- ggplot(plot_data, aes(x = prop_level_cluster, y = !!sym(annot))) + 
    geom_bar(aes(fill = !!sym(lvl_analysis)), stat = "identity") + 
    geom_text(data = plot_labels,
              aes(x  = 1.03, y = !!sym(annot), label = label_level), hjust = 0,
              size = 2.5) +
    coord_cartesian(clip = "off") + 
    scale_fill_manual(values = col_palette) +  
    theme(panel.background = element_blank(), 
          legend.position = "none") + 
    xlim(0, 1.5) +
    ylab("Cell type") +
    xlab(x_axis_title)  +
    guides(fill = guide_legend(reverse = T))
  
  return(prop_plot)
  
}

# Visualization 


# Helper and utility functions for plotting



#' Get default ggplot2/Seurat colours
#'
#' Get evenly spaced colours from around the colour wheel, which are the default
#' colours assigned to clusters by Seurat. The output of this function can be
#' passed to the \code{scale_colour_manual()} and \code{scale_fill_manual()} functions
#' from ggplot2, as the \code{values} argument. (\code{\link{ggColors}} points
#' to this function.)
#'
#' @param n Number of colours to return
#'
#' @return Named character vector, where names are the names of clusters, from
#' 0 to n-1, and values are the hex codes for the colours.
#' @export
#'
#' @examples
#'
#' n_clust <- 5
#' ggColours(n_clust)
#'
#' @references https://stackoverflow.com/a/8197703
#' @aliases ggColors
#' @importFrom grDevices hcl
ggColours <- function(n) {
  
  hues <- seq(15, 375, length = n + 1)
  colours <- hcl(h = hues, l = 65, c = 100)[1:n]
  names(colours) <- seq(0, n - 1) # Since the first cluster in Seurat is 0
  
  return(colours)
  
}

#' @export
ggColors <- ggColours




#' Get a vector of cluster colours, optionally highlighting select clusters
#'
#' Given a seurat object, get a named character vector of cluster colours, where
#' names are cluster names (coresponding to \code{levels(seurat@@ident)}), and
#' values are hex codes of the colours, either the default colours from Seurat,
#' or colours you specify. This is trivial to make this yourself,
#' so this is used as a utility function for retaining the colours of
#' clusters you want to highlight, and setting the colours of all other clusters
#' to grey, or another default, non-intrusive colour.
#'
#' @param seurat Seurat object
#' @param clusters Vector of one or more clusters to highlight, matching the levels at
#' \code{levels(seurat@@ident)}. If "none", returns \code{default_colour}
#' for every clusters. Default: all clusters, obtained from \code{levels(seurat@@ident)}.
#' @param original_colours  (Optional) Vector of colours to use. Either one colour
#' per cluster, in the order of \code{levels(seurat@@ident)}, or one colour per
#' cluster passed to \code{clusters}, in the other they were provided.
#' Default: default ggplot2 colours used by Seurat.
#' @param default_colour Colour to use for non-highlighted clusters
#' Default: gray80 (light grey).
#'
#' @return Named character vector
#' @export
#' @author Selin Jessa
#' @examples
#'
#' # Trivial: get named character vector with default colours
#' getClusterColours(pbmc)
#'
#' # Highlight clusters 2 and 3
#' getClusterColours(pbmc, clusters = c(2, 3))
#'
#' # Highlight clusters 2 and 3, set all other cluster colours to white
#' getClusterColours(pbmc, clusters = c(2, 3), default_colour = "white")
getClusterColours <- function(seurat, clusters = NULL,
                              original_colours = NULL, default_colour = "gray80") {
  
  # Handle clusters argument
  if (is.numeric(clusters)) clusters <- clusters + 1
  if (is.null(clusters)) clusters <- levels(Idents(seurat))
  
  n_clust <- length(levels(Idents(seurat)))
  
  highlight_colours <- rep(default_colour, n_clust)
  names(highlight_colours) <- levels(Idents(seurat))
  
  if (length(clusters) == 1) {
    if (clusters == "none") return(highlight_colours)
  }
  
  # Check if enough orig colours provided
  if (!is.null(original_colours) && length(original_colours) != length(levels(Idents(seurat)))) {
    
    if (length(original_colours) != length(clusters)) {
      
      stop("Not enough 'original_colours'! Please provide as many colours ",
           "as clusters in the dataset, or one per cluster specified in the ",
           "'clusters' argument.")
      
    } else highlight_colours[clusters] <- original_colours
    
  } else if (length(original_colours) == length(levels(Idents(seurat)))) {
    
    highlight_colours[clusters] <- original_colours[clusters]
    
  } else if (is.null(original_colours)) {
    
    original_colours <- ggColours(n_clust)
    names(original_colours) <- levels(Idents(seurat))
    highlight_colours[clusters] <- original_colours[clusters]
    
  }
  
  return(highlight_colours)
  
}



#' Rotate the x axis labels in a ggplot
#'
#' @param angle Integer, value in degrees to rotate labels. Default: 90.
#'
#' @return A theme element to rotate labels
#' @export
#'
#' @author Selin Jessa
#'
#' @examples
#' # gg <- mpg %>%
#' # filter(class %in% c("compact", "suv")) %>%
#' # ggplot(aes(x = displ, y = hwy)) +
#' # geom_point(aes(colour = factor(year))) +
#' # facet_wrap(~ class, ncol = 2)
#'
#' # gg
#' # gg + rotateX()
rotateX <- function(angle = 90) {
  
  theme(axis.text.x = element_text(angle = angle, hjust = 1))
  
}



#' Remove the legend in a ggplot
#'
#' @return A theme element to hide legend
#' @export
#'
#' @author Selin Jessa
#'
#' @examples
#' # gg <- mpg %>%
#' # filter(class %in% c("compact", "suv")) %>%
#' # ggplot(aes(x = displ, y = hwy)) +
#' # geom_point(aes(colour = factor(year))) +
#' # facet_wrap(~ class, ncol = 2)
#'
#' # gg
#' # gg + noLegend()
noLegend <- function() {
  
  theme(legend.position = "none")
  
}



#' Remove axis ticks and tick labels from a ggplot
#'
#' @return A theme element to remove ticks
#' @export
#'
#' @author Selin Jessa
noTicks <- function() {
  
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
  
}


#' clusterCenters
#'
#' Get centers of clusters given a Seurat object, to use for labelling
#' in tSNE space. The cluster center is defined as the median X and Y coordinate
#' across cells in each cluster.
#'
#' @param seurat Seurat object, where dimensionality reduction has been applied,
#' i.e. (after applying Seurat::RunTSNE() to the object).
#'
#' @return Data frame with three columns: Cluster, mean_tSNE_1, and mean_tSNE_2
#' @export
#'
#' @author Selin Jessa
#' @examples
#'
#' clusterCenters(pbmc, reduction = "pca", dim1 = 1, dim2 = 3)
clusterCenters <- function(seurat, reduction, dim1, dim2) {
  
  n_clusters <- length(unique(Idents(seurat)))
  
  # Attempts at tidyeval...
  # vars <- colnames(seurat@reductions[[reduction]]@cell.embeddings)[c(1, 2)]
  # col_names <- paste0("mean_", vars)
  
  # Get the embedding
  df <- as.data.frame(seurat@reductions[[reduction]]@cell.embeddings[, c(dim1, dim2)]) %>%
    mutate(Cell = names(Idents(seurat)),
           Cluster = Idents(seurat))
  
  # Generalize these
  colnames(df)[c(1, 2)] <- c("Dim_1", "Dim_2")
  
  # Compute cluster centers
  centers <- df %>%
    group_by(Cluster) %>%
    summarise(mean_x = median(Dim_1),
              mean_y = median(Dim_2))
  
  return(centers)
  
}


#' Add cluster labels to a tSNE ggplot2 plot
#'
#' @param centers Data frame with at least three columns: "mean_x", "mean_y",
#' and "Cluster", as returned by \code{\link{clusterCenters}}
#' @param label_repel Logical, whether to plot cluster
#' labels repelled from the center, on a slightly transparent white background and
#' with an arrow pointing to the cluster center. If FALSE, simply plot the
#' cluster label at the cluster center. Default: TRUE.
#' @param label_size Numeric, controls the size of text labels. Default: 4.
#' @param label_short (Optional/Experimental!!) Logical, if TRUE, assumes clusters
#' (at \code{seurat@@ident}) consist of a prefix and a suffix separated by a non-alpha
#' numeric character (\code{"[^[:alnum:]]+"}), and tries to separate these names
#' and only plot the prefix, for shorter labels and a cleaner plot. Default: FALSE.
#' @param clusters (Optional) Clusters for which labels should be plot (if only
#' a subset of clusters should be labelled). Default: NULL (Label all clusters).
#'
#'
#' @author Selin Jessa
#' @export
addLabels <- function(centers, label_repel = FALSE, label_size = 4, label_short = FALSE, clusters = NULL) {
  
  if (!is.null(clusters)) centers <- filter(centers, Cluster %in% clusters)
  
  if (label_short) centers <- suppressWarnings(
    tidyr::separate(centers, Cluster, into = c("Cluster", "Cluster_long"), extra = "drop"))
  
  if (label_repel) {
    
    ggrepel::geom_label_repel(data = centers,
                              aes(x = mean_x, y = mean_y),
                              label = centers$Cluster,
                              size = label_size,
                              segment.color = 'grey50',
                              fontface = 'bold',
                              alpha = 0.8,
                              segment.alpha = 0.8,
                              label.size = NA,
                              force = 2,
                              # Leaving these unspecified for now, since it really depends on
                              # the dimensionality reduction
                              # nudge_x = 5, nudge_y = 5,
                              segment.size = 0.5,
                              arrow = arrow(length = unit(0.01, 'npc')))
    
  } else {
    
    geom_text(data = centers,
              aes(x = mean_x, y = mean_y, label = Cluster),
              size = label_size)
    
  }
  
}



#' Get the limits of a the first two dimensions in a dimensionality reduction
#'
#' When plotting an embedding, we may want to plot specific cells, but
#' constrain the scale to match plots of the whole dataset. Given a dim.
#' reduction, this function extracts the x and y limits to use for plotting.
#'
#' @param seurat Seurat object for which a dimensionality reduction has been
#' computed (e.g. PCA or tSNE)
#' @param reduction String, corresponding to the dimensionality reduction to use.
#' Default: "tsne".
#'
#' @return A list with two elements: "xlim", which is a character vector of
#' the limits for the x-axis, and "ylim", correspondingly for the y-axis
#' @export
#' @author Selin Jessa
#'
#' @examples
#' drLims(pbmc)
drLims <- function(seurat, reduction = "tsne") {
  
  dim1 <- seurat@reductions[[reduction]]@cell.embeddings[,1]
  dim2 <- seurat@reductions[[reduction]]@cell.embeddings[,2]
  
  return(list(xlim = c(min(dim1), max(dim1)),
              ylim = c(min(dim2), max(dim2))))
  
}




#' Constrain the scale of the plot to the dimensionality reduction limits
#'
#' @inheritParams drLims
#'
#' @export
#' @author Selin Jessa
constrainScale <- function(seurat, reduction = "tsne")  {
  
  limits <- drLims(seurat = seurat, reduction = reduction)
  lims(x = limits$xlim, y = limits$ylim)
  
  
}



#' Apply a clean theme to a ggplot2 object
#'
#' @references https://github.com/sjessa/ggmin
#'
#' @importFrom ggplot2 theme_light theme
#' @author Selin Jessa
#' @export
theme_min <- function(base_size = 11, base_family = "",
                      border_colour = "grey90",
                      border_size = 1) {
  
  theme_light(base_size = 11, base_family = "") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = border_colour, linewidth = border_size),
      axis.ticks = element_line(colour = border_colour),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black", size = rel(1.2)),
      strip.text.y = element_text(colour = "black", size = rel(1.2)),
      title = element_text(size = rel(0.9)),
      axis.text = element_text(colour = "black", size = rel(0.8)),
      axis.title = element_text(colour = "black", size = rel(0.9)),
      legend.title = element_text(colour = "black", size = rel(0.9)),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA)
    )
}



# Functions for basic visualization of the data
# These plots borrow the idea of the Seurat "FeaturePlot" function,
# with additional options for customization of the visualizations.

#' Plot a reduced dimensionality embedding for a Seurat object, colouring
#' by cluster or arbitrary variables.
#'
#' @param seurat Seurat object, where Seurat::RunTSNE() has been applied
#' @param reduction String, specifying a lot of \code{seurat@@reductions}, which
#' indicates which embedding to plot. Default: "tsne". (Can also take "umap" or "pca").
#' @param colour_by (Optional) String, specifying the column in \code{seurat@@meta.data}
#' by which to colour cells. Default: NULL, colour cells by cluster (in \code{seurat@@ident}).
#' @param colours (Optional) Character vector of colours for points. If \code{colour_by}
#' is NULL, cells will be coloured by cluster; this should be a named character vector of colours for points. Names should
#' correspond to cluster names (e.g. \code{levels(seurat@@ident)}). If
#' specifying \code{colour_by}, and the variable is discrete, this should be a character vector,
#' with either names or order corresponding to categorical
#' values in the column of \code{seurat@@meta.data} specified. Otherwise, if the variable
#' is continuous, pass the gradient to use, or a few colours (from low to high) from which a gradient
#' should be created, and specify \code{colour_by_type = "continuous"}. The default is to
#' use ggplot2 colours.
#' @param colour_by_type (Optional) String, one of "discrete" or "continuous".
#' If specifying \code{colour_by} and providing colours to the \code{colours}
#' argument, specify whether the \code{colour_by} variable is discrete or continuous.
#' Default: discrete. The function is pretty good at setting the palette correctly on
#' its own, but try modifying this in case of errors.
#' @param label Logical, whether to plot cluster labels. Default: TRUE
#' @param point_size Numeric, size of points in scatter plot. Default: 0.6 for datasets
#' with < 300 cells, and 1.3 for datasets with otherwise.
#' @param alpha Numeric, fixed alpha value for points: Default: 0.8
#' @param legend Logical, whether to plot legend. Default: FALSE if \code{colour_by}
#' is NULL and \code{label} is TRUE, true otherwise.
#' @param label_repel Logical, if \code{label} is TRUE, whether to plot cluster
#' labels repelled from the center, on a slightly transparent white background and
#' with an arrow pointing to the cluster center. If FALSE, simply plot the
#' cluster label at the cluster center. Default: TRUE.
#' @param label_size Numeric, controls the size of text labels. Default: 4.
#' @param title (Optional) String specifying title.
#' @param hide_ticks Logical, whether to hide axis ticks, i.e. both the text and the
#' small lines indicating the breaks along the x- and y-axes. Default: FALSE
#' @param hide_axes Logical, whether to hide axis labels. Default: TRUE
#' @param label_short (Optional/Experimental!!) Logical, if TRUE, assumes cluster
#' names (at seurat@@ident) consist of a prefix and a suffix separated by a non-alpha
#' numeric character (\code{"[^[:alnum:]]+"}), and tries to separate these names
#' and only plot the prefix, for shorter labels and a cleaner plot. Default: FALSE.
#' @param highlight_cells Character vector of cell names if only a subset of cells should be
#' coloured in the plot (these should correspond to Cells(seurat)). Default: Plot all cells.
#' See the argument \code{clusters_to_label} for only labelling certain clusters.
#' See the \code{constrain_scale} argument for controlling the scales of the plot.
#' @param show_all_cells Logical. When passing cells to \code{highlight_cells},
#' if TRUE, plot the remaining cells in \code{na_colour}, if FALSE,
#' only plot the highlighted cells. Default: TRUE.
#' @param order_by String, corresponding to a column in \code{seurat@@meta.data}, specifying
#' a variable to control the order in which cells are plot. (Thus, you can manually
#' specify the order, add it as a new column in \code{seurat@@meta.data}, and pass that).
#' If numeric, cells with high values are plot on top. If not, the column must
#' be a factor, and cells will be ordered according to the levels, with cells
#' in the first level plot on top. Default: if a numeric column is specified
#' to \code{colour_by}, sort by that variable, otherwise, use the ordering of the cells
#' in the Seurat object.
#' @param clusters_to_label (Optional.) If \code{label} is TRUE,
#' clusters for which labels should be plot (if only a subset of clusters should be labelled).
#' Default: NULL (Label all clusters).
#' @param na_colour String, specifying the colour (built-in or hex code) to use to
#' plot points which have an NA value, for example
#' in the variable specified in \code{colour_by}. Default: light gray ("gray80),
#' change to "white" to purposely hide those cells. If you do not want to plot
#' certain cells at all, pass names of cells to plot to the \code{cells} argument.
#' @param limits Numeric vector of length two providing the lower and upper limits of
#' the colour scale, if colouring by a continuous variable. Default: min and max
#' of the values the variable takes on in the data.
#' @param constrain_scale Logical, if plotting a subset of cells, whether to
#' use the limits of the tSNE embedding computed on the whole dataset (useful
#' for constraining scales across plots while only plotting specific cells).
#' Default: TRUE
#' @param dim1 Numeric, dimension of embedding to use for x-axis. Default = 1.
#' @param dim2 Numeric, dimension of embedding to use for y-axis. Default = 2.
#'
#' @return A ggplot2 object
#' @export
#' @aliases plot_pca plot_tsne plot_umap
#'
#' @author Selin Jessa
#'
#' @examples
#' tsne(pbmc)
#'
#' # Demonstrate label_short:
#' # Set cluster IDs to be too long
#' pbmc2 <- pbmc
#' levels(pbmc2@ident) <- c("1-Cell type A", "2-Cell type B", "3-Cell type C", "4-Cell type D")
#' tsne(pbmc2)
#'
#' # Plot the prefixes only
#' tsne(pbmc2, label_short = TRUE)
#'
#' # Only colour cells in cluster 3
#' plot_dr(pbmc, highlight_cells = whichCells(pbmc, clusters = 3))
#'
#' # Plot the # of genes for only cells in cluster 3, but display all cells on the plot
#' tsne(pbmc, highlight_cells = whichCells(pbmc, clusters = 3), colour_by = "nGene",
#' show_all_cells = FALSE, colours = viridis::viridis(100), colour_by_type = "continuous")
plot_dr <- function(seurat,
                    reduction = "tsne",
                    colour_by = NULL,
                    colours = NULL,
                    colour_by_type = "discrete",
                    label = TRUE,
                    point_size = ifelse(length(Cells(seurat)) > 300, 0.6, 1.3),
                    alpha = 0.8,
                    legend = ifelse((is.null(colour_by)) && (label), FALSE, TRUE),
                    label_repel = TRUE,
                    label_size = 4,
                    highlight_cells = NULL,
                    show_all_cells = TRUE,
                    order_by = NULL,
                    clusters_to_label = NULL,
                    hide_ticks = TRUE,
                    title = NULL,
                    label_short = FALSE,
                    na_colour = "gray80",
                    limits = NULL,
                    constrain_scale = TRUE,
                    hide_axes = FALSE,
                    dim1 = 1,
                    dim2 = 2,
                    border_colour = NULL,
                    border_size = NULL) {
  
  if (!(reduction %in% names(seurat@reductions))) stop(reduction, " reduction has not been computed.")
  
  # Get the data
  embedding <- data.frame(Cell = Cells(seurat),
                          dim1 = seurat@reductions[[reduction]]@cell.embeddings[, dim1],
                          dim2 = seurat@reductions[[reduction]]@cell.embeddings[, dim2],
                          Cluster = Idents(seurat),
                          stringsAsFactors = FALSE)
  
  # Control cell ordering (in the z-axis)
  if (!is.null(order_by)) {
    
    # Check the variable is usable for sorting
    if (!is.numeric(seurat@meta.data[[order_by]]) && !is.factor(seurat@meta.data[[order_by]])) {
      
      stop("The variable specified in 'order_by' is neither numeric ",
           "nor a factor. If the column is of type character, consider ",
           "converting it to a factor. Otherwise, pass the name of a numeric column.")
      
    }
    
    embedding[[order_by]] <- seurat@meta.data[[order_by]]
    
  } else if ((!is.null(colour_by)) && is.numeric(seurat@meta.data[[colour_by]])) {
    
    # If order_by is not specified but colour_by is, and is numeric,
    # by default, order the cells by that variable
    order_by <- colour_by
    
  }
  
  if (!is.null(colour_by)) embedding[[colour_by]] <- seurat@meta.data[[colour_by]]
  
  # Deal with which cells to plot
  if (!is.null(highlight_cells)) {
    
    # Show all cells, but only colour the highlighted ones, by setting the colour_by
    # value to NA for all other cells
    if (show_all_cells) {
      
      if (is.null(colour_by)) embedding[!(embedding$Cell %in% highlight_cells), ]$Cluster <- NA
      else embedding[!(embedding$Cell %in% highlight_cells), ][[colour_by]] <- NA
      
    } else { # Otherwise, only display and colour the highlighted ones
      
      embedding <- embedding %>% filter(Cell %in% highlight_cells)
      
    }
  }
  
  # We want to sort point such that any NAs will be plot first/underneath
  colour_by2 <- ifelse(is.null(colour_by), "Cluster", colour_by)
  if (!is.null(order_by)) embedding <- embedding %>% arrange(paste0("!is.na(", colour_by2, ")"), order_by)
  else embedding <- embedding %>% arrange(paste0("!is.na(", colour_by2, ")"))
  
  # Start the plot!
  gg <- ggplot(embedding, aes(x = dim1, y = dim2))
  
  if (label && all(is.na(Idents(seurat)))) {
    label <- FALSE
    message("NOTE: identity of all cells is NA, setting 'label' to FALSE.")
  }
  
  # Deal with the palette
  if (is.null(colour_by)) {
    
    if (is.null(colours)) {
      
      # Check if there's a cluster palette stored with the object
      if (!is.null(seurat@misc$colours)) colours <- seurat@misc$colours
      else {
        
        # Assuming that the order of the levels is correct in the seurat object,
        # find the default colours for the clusters
        colours <- ggColors(length(levels(Idents(seurat))))
        names(colours) <- levels(Idents(seurat))
        
      }
      
    }
    
    gg <- gg +
      geom_point(aes(colour = Cluster), size = point_size, alpha = alpha) +
      scale_color_manual(values = colours, na.value = na_colour)
    
  } else {
    
    if (is.null(limits)) lims <- c(NA, NA)
    else lims <- limits
    
    gg <- gg +
      geom_point(aes_string(colour = colour_by), size = point_size, alpha = alpha)
    
    if (!is.null(colours)) {
      
      if (colour_by_type == "discrete") gg <- gg + scale_color_manual(values = colours, na.value = na_colour)
      else if (colour_by_type == "continuous") {
        
        gg <- gg + scale_color_gradientn(colours = colours,
                                         na.value = na_colour,
                                         limits = lims)
      }
      
    } else {
      
      if (colour_by_type == "continuous") { # Otherwise for discrete, default ggplot2 colours are used
        
        gg <- gg + scale_color_gradientn(colours = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "OrRd"))(n = 100),
                                         na.value = na_colour,
                                         limits = lims)
        
      }
    }
  }
  
  # Label clusters
  if (label) {
    
    centers <- clusterCenters(seurat, reduction = reduction, dim1 = dim1, dim2 = dim2)
    gg <- gg + addLabels(centers     = centers,
                         label_repel = label_repel,
                         label_size  = label_size,
                         label_short = label_short,
                         clusters    = clusters_to_label)
    
  }
  
  gg <- gg + theme_min(border_colour = border_colour, border_size = border_size)
  
  # Set the right axes titles
  if (hide_axes) gg <- gg + xlab(NULL) + ylab(NULL)
  else if (reduction == "tsne") gg <- gg + xlab(glue("tSNE {dim1}")) + ylab(glue("tSNE {dim2}"))
  else if (reduction == "umap") gg <- gg + xlab(glue("UMAP {dim1}")) + ylab(glue("UMAP {dim2}"))
  else if (reduction == "phate") gg <- gg + xlab(glue("PHATE {dim1}")) + ylab(glue("PHATE {dim2}"))
  else if (reduction == "pca") {
    
    var_exp <- getVarianceExplained(seurat)
    gg <- gg +
      xlab(glue("PC{dim1} ({round(var_exp$percent.var.explained[dim1], 1)}%)")) +
      ylab(glue("PC{dim2} ({round(var_exp$percent.var.explained[dim2], 1)}%)"))
    
  } else {
    gg <- gg + xlab(glue("{reduction} {dim1}")) + ylab(glue("{reduction} {dim2}"))
  }
  
  # More aesthetics
  if (!legend) gg <- gg + noLegend()
  else if (!is.null(colour_by)) {
    if (colour_by == "orig.ident") gg <- gg + labs(colour = "Sample")
  }
  
  if (!is.null(title)) gg <- gg + ggtitle(title)
  if (hide_ticks) gg <- gg + noTicks()
  if (constrain_scale) gg <- gg + constrainScale(seurat, reduction = reduction)
  
  return(gg)
  
}

#' @describeIn plot_dr Plot a tSNE embedding
#' @export
tsne <- function(seurat, ...) {
  
  plot_dr(seurat, reduction = "tsne", ...)
  
}

#' @describeIn plot_dr Plot a PCA embedding
#' @export
pca <- function(seurat, ...) {
  
  plot_dr(seurat, reduction = "pca", ...)
  
}


#' @describeIn plot_dr Plot a UMAP embedding
#' @export
umap <- function(seurat, ...) {
  
  plot_dr(seurat, reduction = "umap", ...)
  
}

#' @describeIn plot_dr Plot a PHATE embedding
#' @export
phate <- function(seurat, ...) {
  
  plot_dr(seurat, reduction = "phate", ...)
  
}




# Cell Cucle Assignment with custom marker genes 

compute_cell_cycle_genes <- function(seurat, 
                                     facets=TRUE, 
                                     legend=FALSE, 
                                     return_scores = FALSE,
                                     cell.cycle.genes) {
  
  
  g1.s.genes <- dplyr::filter(cell.cycle.genes, phase=="G1/S") %>% 
    .$hsapiens.gene.symbol
  g2.m.genes <- dplyr::filter(cell.cycle.genes, phase=="G2/M") %>% 
    .$hsapiens.gene.symbol
  
  expression.data <- as.data.frame(as.matrix(GetAssayData(seurat, 
                                                          assay = "RNA",
                                                          layer = "data")))
  
  expression.data.g1.s.genes <- dplyr::filter(expression.data, rownames(expression.data) %in% g1.s.genes)
  expression.data.g2.m.genes <- dplyr::filter(expression.data, rownames(expression.data) %in% g2.m.genes)
  
  expression.data.g1.s.scores <- colMeans(expression.data.g1.s.genes)
  expression.data.g2.m.scores <- colMeans(expression.data.g2.m.genes)
  
  cell.cycle.scores <- as.data.frame(rbind(expression.data.g1.s.scores, expression.data.g2.m.scores))
  
  rownames(cell.cycle.scores) <- gsub("expression.data.", "", rownames(cell.cycle.scores))
  
  cell.cycle.scores.tidy <- as.data.frame(t(cell.cycle.scores))
  cell.cycle.scores.tidy <- tibble::rownames_to_column(cell.cycle.scores.tidy, "cell")
  cell.cycle.scores.tidy <- tibble::add_column(cell.cycle.scores.tidy, cluster=Idents(seurat), .after="cell")
  
  if (return_scores) return(cell.cycle.scores.tidy)
  
  # Plots
  p <- ggplot(cell.cycle.scores.tidy, aes(x=g1.s.scores, y=g2.m.scores)) +
    geom_point(aes(color=cluster)) + xlab("G1/S score") + ylab("G2/M score")
  
  if(facets) {
    p <- p + facet_grid(~cluster) }
  
  if(legend==FALSE) {
    p <- p + theme(legend.position="none")
  }
  
  return(p)
}



whichCells <- function(seurat, clusters) {
  
  names(Idents(seurat))[Idents(seurat) %in% clusters]
  
}


rdbu <- rev(c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7", 
          "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC"))

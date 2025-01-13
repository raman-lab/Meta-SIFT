library(cluster)
library(ggplot2)
library(factoextra)
library(RColorBrewer)
library(gplots)
library(pheatmap)
library(dendextend)
library(grid)
library(ComplexHeatmap)
library(dplyr)
library(magrittr)
library(ggrepel)
library(Rtsne)
library(circlize)
library(RColorBrewer)
library(magrittr)
library(xlsx)

set.seed(123)

#importing the dataset
clustering <- read.csv("clusteringmotif.csv", row.names = 1)
df <- clustering

#Agglomerative Hierarchical Clustering using Complex Heatmaps - Blue White Red
#Spearman - Ward.D2 is used for clustering
#ensure the colors used for the legend in colorRamp2 and the scale used for the heatmap_legend_param is accurate for each dataset
p = Heatmap(as.matrix(df),
            column_title = "Hosts", #literally the column title
            clustering_distance_rows = "spearman", #distance method
            clustering_distance_column = "spearman", #distance method
            clustering_method_rows = "ward.D2", #clustering method
            clustering_method_columns = "ward.D2", #clustering method
            show_row_names = FALSE, #No row names are shown as there are too many rows to be legible
            col = colorRamp2(c(-10, 0, 10), c("#2166AC", "#F7F7F7", "#B2182B")), #Choose the color scale size and color here
            row_dend_width = unit(4, "cm"), #defines the width of the dendrogram
            width = unit(4, "cm"), #defines the width of the heatmap for each host
            row_split = 50, #critically, splits the dendrogram into the specified number of clusters
            row_gap = unit(3, "mm"), #defines the gap between each cluster to make it more interpretable
            heatmap_legend_param = list(title = "log2 Fn", at = c(-7, -3, 0, 3, 7)), #defines the title and numeric scale in the legend
            row_title = "Cluster %s", #specifies the title for each cluster
            row_title_rot = 0, #rotates the cluster title for legibility
            row_title_side = "right", #puts the cluster titles on the right for legibility
            row_title_gp = gpar(fontsize = 8)) #defines the font size of the cluster titles
finalheatmap = draw(p, background = "light gray", row_title = "Variant Dendrogram") #draws the heatmap, adds a background color and a title to the dendrogram


# Individual variables so these can be pulled out later if needed, although this does look hideous

#Writing different distance values
spear.cor <- get_dist(df, method = "spearman") #this is a correlation-based distance measurement that should handle outliers better and is our default choice for the phage datasets.
pear.cor <- get_dist(df, method = "pearson") #this is a correlation-based distance but handles outliers more poorly than spearman.
euc.cor <- get_dist(df, method = "euclidean") #this is not a correlation-based distance and should not be ideal for this data set, but is included here for comparison.
man.cor <- get_dist(df, method = "manhattan") #this is not a correlation-based distance and should not be ideal for this data set, but is included here for comparison.

#Writing different clustering approaches. Repeated for each method.
spearward.hc <- hclust(d = spear.cor, method = "ward.D2") #Ward is used to cluster as it is considered best practice to minimize distance between clusters, this is the default choice for phage datasets.
spearmax.hc <- hclust(d = spear.cor, method = "complete") #Maximum or complete is similar to ward in producing compact clusters and is an alternative approach.
spearsingle.hc <- hclust(d = spear.cor, method = "single") #Minimum or single linkage makes loose clusters and is provided for comparison.
spearavg.hc <- hclust(d = spear.cor, method = "average") #Average produces, well, average clustering and is provided for comparison.
spearcen.hc <- hclust(d = spear.cor, method = "centroid") #Centroid is another method for producing 'average' clusters and is provided for comparison.

#Calculating cophenetic score to verify clustering
spearward.coph <- cophenetic(spearward.hc)
spearmax.coph <- cophenetic(spearmax.hc)
spearsingle.coph <- cophenetic(spearsingle.hc)
spearavg.coph <- cophenetic(spearavg.hc)
spearcen.coph <- cophenetic(spearcen.hc)

sink('cophenetic_analysis.csv')
#correlation between cophenetic distance and original distance to verify clustering, values over 0.75 are considered good.
writeLines('spearman - ward')
cor(spear.cor, spearward.coph)
writeLines('spearman - max')
cor(spear.cor, spearmax.coph)
writeLines('spearman - single')
cor(spear.cor, spearsingle.coph)
writeLines('spearman - average')
cor(spear.cor, spearavg.coph)
writeLines('spearman - centroid')
cor(spear.cor, spearcen.coph)

# Plotting cluster calculations
p1 <- fviz_nbclust(df, FUN = hcut, method = "wss", 
                   k.max = 25) +
  ggtitle("(A) Elbow method")
p2 <- fviz_nbclust(df, FUN = hcut, method = "silhouette", 
                   k.max = 25) +
  ggtitle("(B) Silhouette method")
p3 <- fviz_nbclust(df, FUN = hcut, method = "gap_stat", k.max = 50, nboot = 500) +
  ggtitle("(C) Gap statistic")

p3 + theme(axis.text.x = element_text(size=10, angle=90))

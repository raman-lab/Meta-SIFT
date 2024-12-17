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
            col = colorRamp2(c(-7, 0, 7), c("#2166AC", "#F7F7F7", "#B2182B")), #Choose the color scale size and color here
            row_dend_width = unit(4, "cm"), #defines the width of the dendrogram
            width = unit(4, "cm"), #defines the width of the heatmap for each host
            row_split = 21, #critically, splits the dendrogram into the specified number of clusters
            row_gap = unit(3, "mm"), #defines the gap between each cluster to make it more interpretable
            heatmap_legend_param = list(title = "log2 Fn", at = c(-7, -3, 0, 3, 7)), #defines the title and numeric scale in the legend
            row_title = "Cluster %s", #specifies the title for each cluster
            row_title_rot = 0, #rotates the cluster title for legibility
            row_title_side = "right", #puts the cluster titles on the right for legibility
            row_title_gp = gpar(fontsize = 8)) #defines the font size of the cluster titles
finalheatmap = draw(p, background = "light gray", row_title = "Variant Dendrogram") #draws the heatmap, adds a background color and a title to the dendrogram


#output
rows <- row_order(finalheatmap)
df2 <- df
for (i in 1:length(rows[]))
{
  r = unlist(rows[i])
  df2[r,6] = i
}

colnames(df2)[6]="clusterid"
write.csv(df2,"clustersummary21.csv")


#out output
r.dend <- row_dend(finalheatmap)
rcl.list <- row_order(finalheatmap)  #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x))  #check/confirm size clusters


clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(GeneID = rownames(mat[rcl.list[[i]],]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
  do.call(rbind, .)

clu.df

# loop to extract genes for each cluster.


# Determining the best clustering approach using cophenetic correlation for Agglomerative Nesting (Hierarchical Clustering)
# I wrote each as individual variables so these can be pulled out later if needed, although this does look hideous

#Writing different distance values
spear.cor <- get_dist(df, method = "spearman") #this is a correlation-based distance measurement that should handle outliers better and is our default choice for the phage datasets.
pear.cor <- get_dist(df, method = "pearson") #this is a correlation-based distance but handles outliers more poorly than spearman.
euc.cor <- get_dist(df, method = "euclidean") #this is not a correlation-based distance and should not be ideal for this data set, but is included here for comparison.
man.cor <- get_dist(df, method = "manhattan") #this is not a correlation-based distance and should not be ideal for this data set, but is included here for comparison.

#Writing different clustering approaches
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

#Writing different clustering approaches
pearward.hc <- hclust(d = pear.cor, method = "ward.D2") #Ward is used to cluster as it is considered best practice to minimize distance between clusters, this is the default choice for phage datasets.
pearmax.hc <- hclust(d = pear.cor, method = "complete") #Maximum or complete is similar to ward in producing compact clusters and is an alternative approach.
pearsingle.hc <- hclust(d = pear.cor, method = "single") #Minimum or single linkage makes loose clusters and is provided for comparison.
pearavg.hc <- hclust(d = pear.cor, method = "average") #Average produces, well, average clustering and is provided for comparison.
pearcen.hc <- hclust(d = pear.cor, method = "centroid") #Centroid is another method for producing 'average' clusters and is provided for comparison.

#Calculating cophenetic score to verify clustering
pearward.coph <- cophenetic(pearward.hc)
pearmax.coph <- cophenetic(pearmax.hc)
pearsingle.coph <- cophenetic(pearsingle.hc)
pearavg.coph <- cophenetic(pearavg.hc)
pearcen.coph <- cophenetic(pearcen.hc)

#correlation between cophenetic distance and original distance to verify clustering, values over 0.75 are considered good.
writeLines('pearson - ward')
cor(pear.cor, pearward.coph)
writeLines('pearson - max')
cor(pear.cor, pearmax.coph)
writeLines('pearson - single')
cor(pear.cor, pearsingle.coph)
writeLines('pearson - average')
cor(pear.cor, pearavg.coph)
writeLines('pearson - centroid')
cor(pear.cor, pearcen.coph)

#Writing different clustering approaches
eucward.hc <- hclust(d = euc.cor, method = "ward.D2") #Ward is used to cluster as it is considered best practice to minimize distance between clusters, this is the default choice for phage datasets.
eucmax.hc <- hclust(d = euc.cor, method = "complete") #Maximum or complete is similar to ward in producing compact clusters and is an alternative approach.
eucsingle.hc <- hclust(d = euc.cor, method = "single") #Minimum or single linkage makes loose clusters and is provided for comparison.
eucavg.hc <- hclust(d = euc.cor, method = "average") #Average produces, well, average clustering and is provided for comparison.
euccen.hc <- hclust(d = euc.cor, method = "centroid") #Centroid is another method for producing 'average' clusters and is provided for comparison.

#Calculating cophenetic score to verify clustering
eucward.coph <- cophenetic(eucward.hc)
eucmax.coph <- cophenetic(eucmax.hc)
eucsingle.coph <- cophenetic(eucsingle.hc)
eucavg.coph <- cophenetic(eucavg.hc)
euccen.coph <- cophenetic(euccen.hc)

#correlation between cophenetic distance and original distance to verify clustering, values over 0.75 are considered good.
writeLines('euclidean - ward')
cor(euc.cor, eucward.coph)
writeLines('euclidean - max')
cor(euc.cor, eucmax.coph)
writeLines('euclidean - single')
cor(euc.cor, eucsingle.coph)
writeLines('euclidean - average')
cor(euc.cor, eucavg.coph)
writeLines('euclidean - centroid')
cor(euc.cor, euccen.coph)

#Writing different clustering approaches
manward.hc <- hclust(d = man.cor, method = "ward.D2") #Ward is used to cluster as it is considered best practice to minimize distance between clusters, this is the default choice for phage datasets.
manmax.hc <- hclust(d = man.cor, method = "complete") #Maximum or complete is similar to ward in producing compact clusters and is an alternative approach.
mansingle.hc <- hclust(d = man.cor, method = "single") #Minimum or single linkage makes loose clusters and is provided for comparison.
manavg.hc <- hclust(d = man.cor, method = "average") #Average produces, well, average clustering and is provided for comparison.
mancen.hc <- hclust(d = man.cor, method = "centroid") #Centroid is another method for producing 'average' clusters and is provided for comparison.

#Calculating cophenetic score to verify clustering
manward.coph <- cophenetic(manward.hc)
manmax.coph <- cophenetic(manmax.hc)
mansingle.coph <- cophenetic(mansingle.hc)
manavg.coph <- cophenetic(manavg.hc)
mancen.coph <- cophenetic(mancen.hc)

#correlation between cophenetic distance and original distance to verify clustering, values over 0.75 are considered good.
writeLines('manhattan - ward')
cor(man.cor, manward.coph)
writeLines('manhattan - max')
cor(man.cor, manmax.coph)
writeLines('manhattan - single')
cor(man.cor, mansingle.coph)
writeLines('manhattan - average')
cor(man.cor, manavg.coph)
writeLines('manhattan - centroid')
cor(man.cor, mancen.coph)
sink()


#Assessing different clustering methods using the agglomerative coefficient (AC), closer to 1 indicates a more balanced clustered structure
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
# function to compute coefficient
ac <- function(x) {
  agnes(df, method = x)$ac
}
# get agglomerative coefficient for each linkage method
purrr::map_dbl(m, ac)


# Plot cluster results
p1 <- fviz_nbclust(df, FUN = hcut, method = "wss", 
                   k.max = 25) +
  ggtitle("(A) Elbow method")
p2 <- fviz_nbclust(df, FUN = hcut, method = "silhouette", 
                   k.max = 25) +
  ggtitle("(B) Silhouette method")
p3 <- fviz_nbclust(df, FUN = hcut, method = "gap_stat", k.max = 50, nboot = 500) +
  ggtitle("(C) Gap statistic")


p3 + theme(axis.text.x = element_text(size=10, angle=90))

# Display plots side by side
gridExtra::grid.arrange(p1, p2, nrow = 1)

gridExtra::grid.arrange(p3, nrow = 1, cex   =0.5)


dend_plot <- fviz_dend(spearward.hc)
dend_data <- attr(dend_plot, "dendrogram")
dend_cuts <- cut(dend_data, h = 21)
fviz_dend(dend_cuts$lower[[2]])
sub_grp <- cutree(spearward.hc, k = 21)
table(sub_grp)
# Plot full dendogram
fviz_dend(
  spearward.hc,
  k = 21,
  horiz = TRUE,
  rect = TRUE,
  rect_fill = TRUE,
  rect_border = "jco",
  k_colors = "jco",
  cex = 0.1
)

#plot sub-dendrogram
dendsub1 <- fviz_dend(dend_cuts$lower[[1]])
dendsub2 <- fviz_dend(dend_cuts$lower[[1]], type = 'circular')
# Side by side plots
gridExtra::grid.arrange(dendsub1, dendsub2, nrow = 1)
fviz_cluster(list(data = df, cluster = sub_grp))


clusGap(df, FUNcluster, K.max, B = 100, d.power = 1,
        spaceH0 = c("scaledPCA", "original"),
        verbose = interactive())


fviz_cluster(list(data = df, cluster = grp),
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"), 
             ellipse.type = "convex", # Concentration ellipse
             repel = TRUE, # Avoid label overplotting (slow)
             show.clust.cent = FALSE, ggtheme = theme_minimal())

rowfinal = row_order(final)

sorted <- as.data.frame(rowfinal)


sort(cutree(ht1, k=21))
sorted <- as.data.frame(sorted)

nbclust


round(as.matrix(spear.cor)[1:3, 1:3], 1)

write.csv(as.matrix(spear.cor,"spear_cor.csv"))
write.csv(pear.cor,"pear_cor.csv")
write.csv(euc.cor,"euc.cor.csv")
write.csv(man.cor,"man.cor.csv")






                   
                   
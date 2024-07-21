rm(list=ls())
library(Biobase)
library(GEOquery)
library(limma) 
library(mogene10sttranscriptcluster.db) 
library(RColorBrewer) 
library(gProfileR)
library(genefilter) 
library(openxlsx)
library(gplots) 

# Load the data from GEO for GSE76999
GEO_List <- getGEO(GEO="GSE62628") 
# Check the type of this variable 
class(GEO_List)

# Check the length of this list
length(GEO_List)

# An element of a list always comes together with a name, right?
names(GEO_List)

# Since there is only one ExpressionSet, let's assign it into another variable
eset <- GEO_List[[1]]

varLabels(eset)

eset$title
eset$organism_ch1

# Select the samples we need
sel <- 1:10
eset <- eset[ ,sel]
# We also need to modify its expression matrix 
exprs(eset) <- exprs(eset)[ ,sel]

# Create labels for later steps
labels <- c("Exercise","Sedentary") 
sel_groups <- factor(rep(labels, each=5))

palette(c("red","blue", "green"))


boxplot(exprs(eset), col=sel_groups, las=2, ces.axis=0.7, boxwex=0.6, outline=FALSE)

# normalizeBetweenArrays is a function in LIMMA
exprs_matrix.norm <- normalizeBetweenArrays(exprs(eset))
boxplot(exprs_matrix.norm, col=sel_groups, las=2, cex.axis = 0.7, boxwex=0.6, outline=FALSE)

# transform into log2 scale
ex <- exprs(eset) 
exprs(eset) <- log2(ex)

# set up the data and proceed with analysis
eset$description <- sel_groups # Replace the long description by short labels

# Construct design matrices by model.matrix
# model.matrix loads the columns of eset and matches to its description.
design <- model.matrix(~ description + 0, eset) 
# 0 defines the way we compare the samples 
colnames(design) <- levels(sel_groups)

eset$description

design

# Fit linear model for each gene given a series of arrays
fit <- lmFit(eset, design)
# Build comparison and compute the satistics
cont.matrix <- makeContrasts(Exercise-Sedentary, levels=design) # Apply the contrast matrix to the linear model
fit2 <- contrasts.fit(fit, cont.matrix)
# Generate the statistics
fit2 <- eBayes(fit2, 0.05)

all_de_genes <- topTable(fit2, coef=1, adjust="fdr", p.value=0.05, number=Inf)
all_genes_exercise <- topTable(fit2, coef=1, adjust="fdr", p.value=1, number=Inf)

# Output in an excel file
wb <- createWorkbook()

addWorksheet(wb, "ExVSedDEGs")
writeData(wb, "ExVSedDEGs", exercise_v_sedentary_degs, rowNames = FALSE)
saveWorkbook(wb, "DE_genes_GSE62628.xlsx", overwrite = TRUE)

columns(mogene10sttranscriptcluster.db)

# Get all the probe ID in our ExpressionSet
ID <- featureNames(eset)

# Select the symbols according to the IDs we just got
tmp <- select(mogene10sttranscriptcluster.db, ID, c("SYMBOL", "ENTREZID"))

# Then we can insert a column for symbols which match the ID for each row

all_genes_exercise$ENTREZID <- tmp$ENTREZID[match(all_genes_exercise$ID, tmp$PROBEID)]
all_de_genes$symbol <- tmp$SYMBOL[match(all_de_genes$ID, tmp$PROBEID)]

# By Limma build-in function
ma <- fit2[,"Exercise - Sedentary"] 
limma::plotMA(ma)

volcanoplot(fit2, coeff=1, xlim=c(-2,2), col="gray") 
title("Exercise - Sedentary")

# Performs a principal components analysis on the given data matrix
pca <- prcomp(t(ex))
# Take PC1 and PC2 for the plot
plot(pca$x[,1:2],col=sel_groups, pch=19)
# include a legend for points
legend("topright", inset=.2, labels, pch=19, col=1:4, horiz=TRUE)

# ex is the matrix for all the expression values
# all_de_genes is all the DE genes
# Let's filter ex by the list in all_de_genes
sel <- match(all_de_genes$ID, rownames(ex)) # Get the index of DE genes in ex 
sel_ex <- ex[sel,] # Select those rows only

clust <- function(x) hclust(x, method="complete") # perfrom complete linkage clustering
dist <- function(x) as.dist((1-cor(t(x)))/2) # use the inverse of correlation as distanc
hmcol <- rev(brewer.pal(11, "RdBu"))

hp <- heatmap.2(sel_ex, # input matrix
                col=hmcol,
                scale='row',
                labRow=F,
                density.info='none',
                trace='none',
                margins = c(8,3),
                labCol = rep(labels,each=5),
                hclust = clust,
                distfun = dist)

# Extract the hierarchical cluster from heatmap to class "hclust"
hc <- as.hclust(hp$rowDendrogram)
# Cut the tree by the dedired number of groups
clusters <- cutree(hc, k=5)

clustered_genes <- list(all_de_genes[names(clusters[clusters==1]),]$symbol, 
                        all_de_genes[names(clusters[clusters==2]),]$symbol, 
                        all_de_genes[names(clusters[clusters==3]),]$symbol, 
                        all_de_genes[names(clusters[clusters==4]),]$symbol, 
                        all_de_genes[names(clusters[clusters==5]),]$symbol)

wb <- createWorkbook()
addWorksheet(wb, "cluster 1")
addWorksheet(wb, "cluster 2")
addWorksheet(wb, "cluster 3")
addWorksheet(wb, "cluster 4")
addWorksheet(wb, "cluster 5")
writeData(wb, "cluster 1", clustered_genes[[1]], rowNames = FALSE) 
writeData(wb, "cluster 2", clustered_genes[[2]], rowNames = FALSE) 
writeData(wb, "cluster 3", clustered_genes[[3]], rowNames = FALSE) 
writeData(wb, "cluster 4", clustered_genes[[4]], rowNames = FALSE) 
writeData(wb, "cluster 5", clustered_genes[[5]], rowNames = FALSE) 
saveWorkbook(wb, "gene_clusters.xlsx", overwrite = TRUE)


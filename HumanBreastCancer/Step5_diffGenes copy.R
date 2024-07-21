# Introduction to this script -----------
# the goal of this script is to identify differentially expressed genes (DEGs) and differential transcript usage (DTU)
# you should already know which pairwise comparisons are most important to you
# whether you look for differential expression at the gene or transcript level depends on how you read the Kallisto output into R using TxImport back in Step 1
# if you have no biological replicates, you will NOT be able to leverage statistical tools for differential expression analysis
# instead, you will ONLY rely on fold changes, and can use the dplyr 'verbs' we discussed in Step 3 and 4 to identify genes based on log fold-change

# Load packages -----
library(tidyverse) # you know it well by now!
library(limma) # venerable package for differential gene expression using linear modeling
library(edgeR)
library(gt) 
library(DT) 
library(plotly) 

# Set up Irwin study data like Steps 1-3 of DIYTransciptomics ----

targets <- studyDesign
sampleLabels <- targets$Sample_title
log2.cpm <- archs4.filtered.norm.log2.cpm

#make into dataframe
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")

# use the tidy package to 'pivot' your dataframe (from wide to long)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = GSM3714577:GSM3714640, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)


# note it is easy to plot this pivoted data
p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="Filtered, TMM-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()


#Hierarchical Clustering

distance <- dist(t(log2.cpm), method = "maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters <- hclust(distance, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(clusters, labels=sampleLabels)

# Load packages ------
library(tidyverse) # you're familiar with this fromt the past two lectures
library(DT) # for making interactive tables
library(plotly) # for making interactive plots
library(gt) # A layered 'grammar of tables' - think ggplot, but for tables

# Use dplyr 'verbs' to modify our dataframe ----
# use dplyr 'mutate' function to add new columns based on existing data

PostExvMind.df <- log2.cpm.df %>% 
  mutate(PostExercise.AVG = (GSM3714580	+ GSM3714582 + GSM3714584	+ GSM3714590 +	GSM3714592	+ GSM3714594	+ GSM3714598	+ GSM3714602	+ GSM3714604	+ GSM3714608	+ GSM3714614	+ GSM3714620	+ GSM3714624	+ GSM3714632	+ GSM3714636	+ GSM3714640)/16,
         PostMind.AVG = (GSM3714578	+ GSM3714586	+ GSM3714588	+ GSM3714596	+ GSM3714600	+ GSM3714606	+ GSM3714610	+ GSM3714612	+ GSM3714616	+ GSM3714618	+ GSM3714622	+ GSM3714626	+ GSM3714628	+ GSM3714630	+ GSM3714634	+ GSM3714638)/16,
         #now make columns comparing each of the averages above that you're interested in
         LogFC = (PostExercise.AVG - PostMind.AVG)) %>% 
  mutate_if(is.numeric, round, 2)

PostExvPreEx.df <- log2.cpm.df %>% 
  mutate(PostExercise.AVG = (GSM3714580	+ GSM3714582 + GSM3714584	+ GSM3714590 +	GSM3714592	+ GSM3714594	+ GSM3714598	+ GSM3714602	+ GSM3714604	+ GSM3714608	+ GSM3714614	+ GSM3714620	+ GSM3714624	+ GSM3714632	+ GSM3714636	+ GSM3714640)/16,
         PreExercise.AVG = (GSM3714579	+ GSM3714581	+ GSM3714583	+ GSM3714589	+ GSM3714591	+ GSM3714593	+ GSM3714597	+ GSM3714601	+ GSM3714603	+ GSM3714607	+ GSM3714613	+ GSM3714619	+ GSM3714623	+ GSM3714631	+ GSM3714635	+ GSM3714639)/16,
         #now make columns comparing each of the averages above that you're interested in
         LogFC = (PostExercise.AVG - PreExercise.AVG)) %>% 
  mutate_if(is.numeric, round, 2)


#now look at this modified data table
PostExvMind.df

# Use dplyr 'arrange' and 'select' to sort your dataframe based on any variable
# first, we'll use dplyr "arrange" function to sort rows based on the values in a column of interest
# then we'll display 'select' only the columns we're interested in seeing
PostExvMind.df.sort <- PostExvMind.df %>%
  dplyr::arrange(desc(LogFC)) %>% 
  dplyr::select(geneID, LogFC)

PostExvPreEx.df.sort <- PostExvPreEx.df %>%
  dplyr::arrange(desc(LogFC)) %>% 
  dplyr::select(geneID, LogFC)

# Use dplyr "filter" and "select" functions to pick out genes of interest 
# ways to tweak the 'select' function:
# use ':' between two column names to select all columns between
# use 'contains', 'starts_with' or 'ends_with' to modify how you select
# can refer to columns using exact name or numerical indicator
# use boolean operators such as '&' (and), '|' (or), '==' (equal to), '!' (not)
#mydata.filter <- mydata.df %>%
#  dplyr::filter(geneID=="SLC2A1" | geneID=="SLC2A3" | geneID=="IL1B" | geneID=="GNLY" | geneID=="IFNG"
#                | geneID=="CCL4" | geneID=="KIR2DL4" | geneID=="PRF1" | geneID=="APOBEC3A" | geneID=="UNC13A" ) %>%
#  dplyr::select(geneID, healthy.AVG, disease.AVG, LogFC) %>%
#  dplyr::arrange(desc(LogFC))

# Make an interactive table using the DT package ----
datatable(PostExvMind.df[,c(1,12:14)], 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100")))

# Make an interactive scatter plot with plotly -----
# begin by storing your ggplot object
PostExvMind.plot <- ggplot(PostExvMind.df) + 
  aes(x=PostExercise.AVG, y=PostMind.AVG) +
  geom_point(shape=16, size=1) +
  ggtitle("Exercise vs. Mind-Body Control") +
  theme_bw()

#now use the ggplotly function from the plotly package to convert this ggplot object into an interactive plot
ggplotly(PostExvMind.plot)

#let's customize this graphic by adding a more informative mouseover tooltip
PostExvMind.plot <- ggplot(PostExvMind.df) +
  aes(x=PostExercise.AVG, y=PostMind.AVG, 
      text = paste("Symbol:", geneID)) +
  geom_point(shape=16, size=1) +
  ggtitle("Exercise vs. Mind-Body Control") +
  theme_bw()

ggplotly(PostExvMind.plot)


# Start of New Block: ----

#Make new "group" column in targets with timepoint+condition

targets$group <- paste(targets$timepoint,targets$condition)

#Set up your design matrix 

group <- factor(targets$group)
design <- model.matrix(~0 + group)
colnames(design) <- c("postex", "postmind", "preex", "premind")

# NOTE: if you need a paired analysis (a.k.a.'blocking' design) or have a batch effect, the following design is useful
# design <- model.matrix(~block + treatment) 
# this is just an example. 'block' and 'treatment' would need to be objects in your environment

# Model mean-variance trend and fit linear model to data ----
# Use VOOM function from Limma package to model the mean-variance relationship
v.archs4.dgelist.filtered.norm <- voom(archs4.dgelist.filtered.norm, design, plot = TRUE)
# fit a linear model to your data
fit <- lmFit(v.archs4.dgelist.filtered.norm, design)

# Contrast matrix ----
contrast.matrix <- makeContrasts(exercise = postex-preex, 
                                 exvs.mind = postex-postmind, 
                                 mindbody = postmind-premind,
                                 levels=design)

# extract the linear model fit -----
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)
#write.fit(ebFit, file="lmfit_results.txt")

# TopTable to view DEGs -----
#coef is the contrast (from the contrast matrix) that you want to compare
exerciseTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=10, sort.by="logFC")
exvs.mindTopHits <- topTable(ebFit, adjust ="BH", coef=2, number=10, sort.by="logFC")
mindbodyTopHits <- topTable(ebFit, adjust ="BH", coef=3, number=10, sort.by="logFC")

# convert to a tibble
exvs.mindTopHits.df <- exvs.mindTopHits %>%
  as_tibble(rownames = "geneID")

gt(exvs.mindTopHits.df)
# TopTable (from Limma) outputs a few different stats:
# logFC, AveExpr, and P.Value should be self-explanatory
# adj.P.Val is your adjusted P value, also known as an FDR (if BH method was used for multiple testing correction)
# B statistic is the log-odds that that gene is differentially expressed. If B = 1.5, then log odds is e^1.5, where e is euler's constant (approx. 2.718).  So, the odds of differential expression os about 4.8 to 1 
# t statistic is ratio of the logFC to the standard error (where the error has been moderated across all genes...because of Bayesian approach)

# Volcano Plots ----
# in topTable function above, set 'number=40000' to capture all genes

# now plot
vplot <- ggplot(exvs.mindTopHits.df) +
  aes(y=-log10(P.Value), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=1) +
  #geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  #geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  #geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "Exercise vs. Mind Post-Intervention",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# Now make the volcano plot above interactive with plotly
ggplotly(vplot)

# decideTests to pull out the DEGs and make Venn Diagram ----
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1)

# take a look at what the results of decideTests looks like
head(results)
summary(results)
vennDiagram(results, include="up")

# retrieve expression data for your DEGs ----
head(v.DEGList.filtered.norm$E)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels

diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
head(diffGenes)
dim(diffGenes)
#convert your DEGs to a dataframe using as_tibble
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

# create interactive tables to display your DEGs ----
datatable(diffGenes.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: DEGs in cutaneous leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:11), digits=2)

#write your DEGs to a file
write_tsv(diffGenes.df,"DiffGenes.txt") #NOTE: this .txt file can be directly used for input into other clustering or network analysis tools (e.g., String, Clust (https://github.com/BaselAbujamous/clust, etc.)

# OPTIONAL: differential transcript usage (DTU) analysis ----
library(IsoformSwitchAnalyzeR)

# The IsoformSwitchAnalyzeR package looks for certain column headers in our study design
# So, the first step is to make sure our study design contains the following:
# unique sample IDs must be contained in column called 'sampleID'
# covariate(s) of interest must be in column labeled 'condition'
# remove extraneous columns
targets.mod <- targets %>%
  dplyr::rename(sampleID = sample, condition = group) %>%
  dplyr::select(sampleID, condition)

# import transcript Kallisto quant data
# using the same path variable we set way back in the step 1 script
Txi_trans <- importIsoformExpression(sampleVector = path)

# fix column headers of abundance and counts data to match sampleID in target.mod
colnames(Txi_trans$abundance) <- c("isoform_id", sampleLabels) 
colnames(Txi_trans$counts) <- c("isoform_id", sampleLabels) 

# import data
mySwitchList <- importRdata(
  isoformCountMatrix   = Txi_trans$counts,
  isoformRepExpression = Txi_trans$abundance,
  designMatrix         = targets.mod,
  removeNonConvensionalChr = TRUE,
  addAnnotatedORFs=TRUE,
  isoformExonAnnotation = "Homo_sapiens.GRCh38.100.chr_patch_hapl_scaff.gtf.gz", #need to download this file
  isoformNtFasta       = "Homo_sapiens.GRCh38.cdna.all.fa.gz",
  showProgress = TRUE)

# We'll do the isoform analysis in one step, but there's a lot to unpack here, so you should really read the package documentation at:
# https://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html
# Note that without additional manual work here (beyond the scope of this class), we'll only capture isoform annotations for 1) intron retention; 2) ORF sequence similarity; and 3) nonsense mediate decay (NMD)
mySwitchList <- isoformSwitchAnalysisCombined(
  switchAnalyzeRlist   = mySwitchList,
  pathToOutput = 'isoform_output', # directory must already exist
  outputSequences = TRUE)  
# now look at the directory that you just created above
# in case you missed the summary output from the function above
extractSwitchSummary(mySwitchList)

# extract the top n isoform switching events 
extractTopSwitches(
  mySwitchList, 
  filterForConsequences = TRUE, # these 'consequences' related to the annotations I reference above.
  n = 50, 
  sortByQvals = FALSE) #change to TRUE if you want this list sorted by FDR-adusted Pval (a.k.a., q value)

# visualize by making a 'switch plot'
switchPlot(
  mySwitchList,
  gene='FCGR1B',
  condition1 = 'disease',
  condition2 = 'healthy',
  localTheme = theme_bw())

# the essentials ----
library(tidyverse)
library(limma)
library(edgeR) 
library(gt)
library(DT) 
library(plotly) 

group <- factor(targets$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.filtered.norm, design)
contrast.matrix <- makeContrasts(infection = disease - healthy,
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

vplot <- ggplot(myTopHits) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "Cutaneous leishmaniasis",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggplotly(vplot)

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")
datatable(diffGenes.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: DEGs in cutaneous leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:11), digits=2)



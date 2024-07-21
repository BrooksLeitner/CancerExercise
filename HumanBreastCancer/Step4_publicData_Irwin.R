# Introduction to this script -----------
# the goal of this script is to give you access to massive amounts of gene expression data without the need to download and analyze individual .fastq files
# this script allows you to access RNAseq data from GEO and SRA using the "All RNA-seq and CHIP-Seq Sample and Signature Search" (ARCHS4) project (https://amp.pharm.mssm.edu/archs4)

# Load packages ------
# nothing new here...you should already have all these packages in your R package library
library(tidyverse)
library(rhdf5)
library(edgeR)

# load ARCHS4 database -----
# you should have already downloaded the most recent versions of mouse and human RNAseq data from ARCHS4 in hdf5 format
# begin by creating file paths that point to the hdf5 archs4 files
archs4.human <- "/Users/brooksleitner/Desktop/Python/PerryLabData/RNA_Sequencing/DIYtranscriptomics/human_matrix_v10.h5" # if you placed the hdf5 file in your working directory, just use "human_matrix_v8.h5" as the path
#archs4.mouse <- "mouse_matrix_v10.h5" # if you placed the hdf5 file in your working directory, just use "human_matrix_v8.h5" as the path
# use the h5 list (h5ls) function from the rhdf5 package to look at the contents of these databases
h5ls(archs4.human)
#h5ls(archs4.mouse)

# 348,184 samples from human
all.samples.human <- h5read(archs4.human, name="meta/samples/geo_accession")

# 284,907 samples from mouse
#all.samples.mouse <- h5read(archs4.mouse, name="meta/samples/geo_accession")

# query ARCHS4 database ----
# choose your samples based on GEO or SRA ID

mySamples <- c("GSM3714577",      #`Study.patient._2_mind_pre-intervention`,
  "GSM3714578",      #"Study patient _2_mind_post-intervention",
  "GSM3714579",   #"Study patient _3_exercise_pre-intervention",
  "GSM3714580",  #"Study patient _3_exercise_post-intervention",
  "GSM3714581",   #"Study patient _4_exercise_pre-intervention",
  "GSM3714582",  #"Study patient _4_exercise_post-intervention",
  "GSM3714583",   #"Study patient _6_exercise_pre-intervention",
  "GSM3714584",  #"Study patient _6_exercise_post-intervention",
  "GSM3714585",       #"Study patient _7_mind_pre-intervention",
  "GSM3714586",      #"Study patient _7_mind_post-intervention",
  "GSM3714587",       #"Study patient _8_mind_pre-intervention",
  "GSM3714588",      #"Study patient _8_mind_post-intervention",
  "GSM3714589",   #"Study patient _9_exercise_pre-intervention",
  "GSM3714590",  #"Study patient _9_exercise_post-intervention",
  "GSM3714591",  #"Study patient _10_exercise_pre-intervention",
  "GSM3714592", #"Study patient _10_exercise_post-intervention",
  "GSM3714593",  #"Study patient _11_exercise_pre-intervention",
  "GSM3714594", #"Study patient _11_exercise_post-intervention",
  "GSM3714595",      #"Study patient _12_mind_pre-intervention",
  "GSM3714596",     #"Study patient _12_mind_post-intervention",
  "GSM3714597",  #"Study patient _15_exercise_pre-intervention",
  "GSM3714598", #"Study patient _15_exercise_post-intervention",
  "GSM3714599",      #"Study patient _16_mind_pre-intervention",
  "GSM3714600",     #"Study patient _16_mind_post-intervention",
  "GSM3714601",  #"Study patient _17_exercise_pre-intervention",
  "GSM3714602", #"Study patient _17_exercise_post-intervention",
  "GSM3714603",  #"Study patient _18_exercise_pre-intervention",
  "GSM3714604", #"Study patient _18_exercise_post-intervention",
  "GSM3714605",      #"Study patient _19_mind_pre-intervention",
  "GSM3714606",     #"Study patient _19_mind_post-intervention",
  "GSM3714607",  #"Study patient _13_exercise_pre-intervention",
  "GSM3714608", #"Study patient _13_exercise_post-intervention",
  "GSM3714609",      #"Study patient _24_mind_pre-intervention",
  "GSM3714610",     #"Study patient _24_mind_post-intervention",
  "GSM3714611",      #"Study patient _21_mind_pre-intervention",
  "GSM3714612",     #"Study patient _21_mind_post-intervention",
  "GSM3714613",  #"Study patient _27_exercise_pre-intervention",
  "GSM3714614", #"Study patient _27_exercise_post-intervention",
  "GSM3714615",      #"Study patient _28_mind_pre-intervention",
  "GSM3714616",     #"Study patient _28_mind_post-intervention",
  "GSM3714617",      #"Study patient _30_mind_pre-intervention",
  "GSM3714618",     #"Study patient _30_mind_post-intervention",
  "GSM3714619",  #"Study patient _31_exercise_pre-intervention",
  "GSM3714620", #"Study patient _31_exercise_post-intervention",
  "GSM3714621",      #"Study patient _36_mind_pre-intervention",
  "GSM3714622",     #"Study patient _36_mind_post-intervention",
  "GSM3714623",  #"Study patient _37_exercise_pre-intervention",
  "GSM3714624", #"Study patient _37_exercise_post-intervention",
  "GSM3714625",      #"Study patient _26_mind_pre-intervention",
  "GSM3714626",     #"Study patient _26_mind_post-intervention",
  "GSM3714627",      #"Study patient _34_mind_pre-intervention",
  "GSM3714628",     #"Study patient _34_mind_post-intervention",
  "GSM3714629",      #"Study patient _35_mind_pre-intervention",
  "GSM3714630",     #"Study patient _35_mind_post-intervention",
  "GSM3714631",  #"Study patient _39_exercise_pre-intervention",
  "GSM3714632", #"Study patient _39_exercise_post-intervention",
  "GSM3714633",      #"Study patient _41_mind_pre-intervention",
  "GSM3714634",     #"Study patient _41_mind_post-intervention",
  "GSM3714635",  #"Study patient _42_exercise_pre-intervention",
  "GSM3714636", #"Study patient _42_exercise_post-intervention",
  "GSM3714637",      #"Study patient _43_mind_pre-intervention",
  "GSM3714638",     #"Study patient _43_mind_post-intervention",
  "GSM3714639",  #"Study patient _44_exercise_pre-intervention",
  "GSM3714640","" #"Study patient _44_exercise_post-intervention"
)


# Identify columns to be extracted from ARCHS4 database
my.sample.locations <- which(all.samples.human %in% mySamples) # first time you've seen the %in% operator.
# extract gene symbols from the metadata
genes <- h5read(archs4.human, "meta/genes/genes")

# Extract expression data from ARCHS4 ----
expression <- h5read(archs4.human, "data/expression", 
                     index=list(my.sample.locations, 1:length(genes)))
#transpose expression data matrix
expression = t(expression)

rownames(expression) <- genes
colnames(expression) <- all.samples.human[my.sample.locations]
colSums(expression) #this shows the sequencing depth for each of the samples you've extracted
archs4.dgelist <- DGEList(expression)
archs4.cpm <- cpm(archs4.dgelist)
colSums(archs4.cpm)

# Filter and normalize the extracted data ----
table(rowSums(archs4.dgelist$counts==0)==64)
keepers <- rowSums(archs4.cpm>1)>=16
archs4.dgelist.filtered <- archs4.dgelist[keepers,]
dim(archs4.dgelist.filtered)
archs4.dgelist.filtered.norm <- calcNormFactors(archs4.dgelist.filtered, method = "TMM")

archs4.filtered.norm.log2.cpm <- cpm(archs4.dgelist.filtered.norm, log=TRUE)

#write to CSV

write.csv(archs4.filtered.norm.log2.cpm,"log2normalizedexpression_Irwinstudy.csv", row.names = TRUE)

# Extract sample metadata from ARCHS4 to create a study design file ----
# extract the sample source
Sample_source_name_ch1 <- h5read(archs4.human, "meta/samples/source_name_ch1")
# extract sample title
Sample_title <- h5read(archs4.human, name="meta/samples/title")
# extract sample characteristics
Sample_characteristics<- h5read(archs4.human, name="meta/samples/characteristics_ch1")

# let's try putting this all together in a study design file
studyDesign <- tibble(Sample_title = Sample_title[my.sample.locations], 
                      Sample_source = Sample_source_name_ch1[my.sample.locations],
                      Sample_characteristics = Sample_characteristics[my.sample.locations])

#based on what we extracted from ARCHS4 above, lets customize and clean-up this study design file
studyDesign <- tibble(Sample_title = Sample_title[my.sample.locations],
                      timepoint = c("pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post","pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post","pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post", "pre", "post"),
                      condition = c("mind",	"mind",	"exercise",	"exercise",	"exercise",	"exercise",	"exercise",	"exercise",	"mind",	"mind",	"mind",	"mind",	"exercise",	"exercise",	"exercise",	"exercise",	"exercise",	"exercise",	"mind",	"mind",	"exercise",	"exercise",	"mind",	"mind",	"exercise",	"exercise",	"exercise",	"exercise",	"mind",	"mind",	"exercise",	"exercise",	"mind",	"mind",	"mind",	"mind",	"exercise",	"exercise",	"mind",	"mind",	"mind",	"mind",	"exercise",	"exercise",	"mind",	"mind",	"exercise",	"exercise",	"mind",	"mind",	"mind",	"mind",	"mind",	"mind",	"exercise",	"exercise",	"mind",	"mind",	"exercise",	"exercise",	"mind",	"mind",	"exercise","exercise"))

#capture experimental variables as factors from this study design
timepoint <- factor(studyDesign$timepoint)
condition <- factor(studyDesign$condition)
sampleName <- studyDesign$Sample_title

# Principal component analysis (PCA) -------------
pca.res <- prcomp(t(archs4.filtered.norm.log2.cpm), scale.=F, retx=T)
#look at pca.res in environment
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
pca.res$rotation #$rotation shows you how much each gene influenced each PC (called 'scores')
pca.res$x #$x shows you how much each sample influenced each PC (called 'loadings')
#note that these loadings have a magnitude and a direction (this is the basis for making a PCA plot)
pc.var<-pca.res$sdev^2 #sdev^2 gives you the eigenvalues
pc.per<-round(pc.var/sum(pc.var)*100, 1)
pc.per

# Visualize your PCA result ------------------
#lets first plot any two PCs against each other
#We know how much each sample contributes to each PC (loadings), so let's plot
pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df) +
  aes(x=PC2, y=PC3, color=condition, shape=timepoint) +
  geom_point(size=4) +
  # geom_label() +
  # stat_ellipse() +
  xlab(paste0("PC2 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC3 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

# now try painting other variables from your study design file onto this PCA.
# can you determine the relationship between PC2 and your metadata?
# can we map one variable to point color and another to point shape?

# now create a small multiple PCA plot
pca.res.df <- pca.res$x[,1:4] %>% # note that this is the first time you've seen the 'pipe' operator from the magrittr package
  as_tibble() %>%
  add_column(condition) %>%
  add_column(timepoint) %>%
  add_column(sampleName)


pca.pivot <- pivot_longer(pca.res.df, # dataframe to be pivoted
                          cols = PC1:PC4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)

ggplot(pca.pivot) +
  aes(x=sampleName, y=loadings, fill=condition) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes. Try doing this with the genotype variable you created above.
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()



# the essentials ----
library(tidyverse)
library(rhdf5)
library(edgeR)

archs4.mouse <- "../../ARCHS4/mouse_matrix_v8.h5" # if you placed the hdf5 file in your working directory, just use "human_matrix_v8.h5" as the path
all.samples.mouse <- h5read(archs4.mouse, name="meta/Sample_geo_accession")
mySamples <- c("GSM2310941", # WT_unstim_rep1
               "GSM2310942", # WT_unstim_rep2
               "GSM2310943", # Ripk3_unstim_rep1
               "GSM2310944", # Ripk3_unstim_rep2
               "GSM2310945", # Ripk3Casp8_unstim_rep1
               "GSM2310946", # Ripk3Casp8_unstim_rep2
               "GSM2310947", # WT_LPS.6hr_rep1
               "GSM2310948", # WT_LPS.6hr_rep2
               "GSM2310949", # Ripk3_LPS.6hr_rep1
               "GSM2310950", # Ripk3_LPS.6hr_rep2
               "GSM2310951", # Ripk3Casp8_LPS.6hr_rep1
               "GSM2310952") # Ripk3Casp8_LPS.6hr_rep2

my.sample.locations <- which(all.samples.mouse %in% mySamples)
genes <- h5read(archs4.mouse, "meta/genes")
expression <- h5read(archs4.mouse, "data/expression", 
                     index=list(1:length(genes), my.sample.locations))

rownames(expression) <- genes
colnames(expression) <- all.samples.mouse[my.sample.locations]
archs4.dgelist <- DGEList(expression)
archs4.cpm <- cpm(archs4.dgelist)

keepers <- rowSums(archs4.cpm>1)>=2
archs4.dgelist.filtered <- archs4.dgelist[keepers,]
archs4.dgelist.filtered.norm <- calcNormFactors(archs4.dgelist.filtered, method = "TMM")
archs4.filtered.norm.log2.cpm <- cpm(archs4.dgelist.filtered.norm, log=TRUE)

Sample_source_name_ch1 <- h5read(archs4.mouse, "meta/Sample_source_name_ch1")
Sample_title <- h5read(archs4.mouse, name="meta/Sample_title")
Sample_characteristics<- h5read(archs4.mouse, name="meta/Sample_characteristics_ch1")

studyDesign <- tibble(Sample_title = Sample_title[my.sample.locations], 
                      genotype = c("WT", "WT", "Ripk3", "Ripk3", "Ripk3Casp8", "Ripk3Casp8", "WT", "WT", "Ripk3", "Ripk3", "Ripk3Casp8", "Ripk3Casp8"),
                      treatment = c("unstim", "unstim", "unstim", "unstim", "unstim", "unstim", "LPS", "LPS", "LPS", "LPS", "LPS", "LPS"))

pca.res <- prcomp(t(archs4.filtered.norm.log2.cpm), scale.=F, retx=T)
pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, color=treatment, shape=genotype) +
  geom_point(size=4) +
  # geom_label() +
  # stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()


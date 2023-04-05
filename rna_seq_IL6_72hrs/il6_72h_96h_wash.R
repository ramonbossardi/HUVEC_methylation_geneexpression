setwd("E:/R/Ramon/IL6_72h_96h_wash_081522")



save.image(file = 'hg38_huvec_72h_96h.Rdata')
load(file = 'hg38_huvec_72h_96h.Rdata')

#Identify many cores have to use 
#library(parallel)
#detectCores() 
#Gabi computer has 28 cores to use

library(Rsubread)
library(gplots)
library(openxlsx)
library(edgeR)
library(limma)
library(DESeq2)
library(genefilter)
library(RColorBrewer)
library(org.Hs.eg.db)
library(cluster)
library(factoextra)
library(clusterProfiler)
library(pathview)
library(sva)
library(systemPipeR)
library(rtracklayer)
library(stringr)
library(GenomicFeatures)

# Read the sample information into R
sampleTable <- read.csv('il6_metadata.csv', header = TRUE, sep = ",")
sampleTable

#how to build
#https://bioinformatics-core-shared-training.github.io/RNAseq-R/align-and-count.nb.html
#http://combine-australia.github.io/RNAseq-R/07-rnaseq-day2.html
buildindex(basename="hg38",reference="hg38.fa")

#Path for the reference
hg38.genome_ref.path   <- "E:/R/Human_rnaseq/hg38_index/hg38.fa"

#Load FASTQ Files to align
fastq_dir <- "E:/R/Ramon/IL6_72h_96h_wash_081522/samples"

#Upload to environment FASTQ files 
reads1 <- list.files(path = "E:/R/Ramon/IL6_72h_96h_wash_081522/samples", pattern = "_1.fq.gz$", full.names = TRUE)
reads2 <- list.files(path = "E:/R/Ramon/IL6_72h_96h_wash_081522/samples", pattern = "_2.fq.gz$", full.names = TRUE)

align( index          = "E:/R/Human_rnaseq/hg38_index/hg38", 
       readfile1      = reads1,
       readfile2      = reads2,
       type           = 'rna', 
       input_format   = "gzFASTQ",
       output_format  = "BAM",
       PE_orientation = 'fr', 
       nthreads       = 28 ) #use the parallel to identify how many cores have to use

#Have a variable with all the bam files
bam.files <- list.files(path = "./samples", pattern = ".BAM$", full.names = TRUE)
bam.files

#The function propmapped returns the proportion of mapped reads in the output SAM file: 
#total number of input reads, number of mapped reads and proportion of mapped reads.
# promapped proportion ranging from 0.05 to 0.150
# https://qfab-bioinformatics.github.io/workshops-RNAseq-analysis-with-R/mapping.html - If need to look a specific sample with prolem (examining the mapped reads)
props <- propmapped(files=bam.files)
props

#Rsubread provides a read summarization function featureCounts, which takes two inputs:
# 1. the aligned reads (BAM or SAM) and assigns them to
# 2. genomic features (GTF annotation file)
seqdata <- featureCounts( bam.files, useMetaFeatures = TRUE, annot.inbuilt = "hg38", isPairedEnd = TRUE,
                          nthreads = 28)

# Examine the attributes in the returned mini.counts object:
summary(seqdata)

#Find the dimension of the counts table:
dim(seqdata$counts)

seqdata$counts[1:6,]

#Look at the annotations, which corresponds to the rows of the counts table:
head(seqdata$annotation)

#featureCounts also returns a very hand summary of the number of reads that were assigned or unassiged:
seqdata$stat

#matrix with the counts
gene.counts <- seqdata$counts
#Save the Ids for the lines
gene.ids    <- seqdata$annotation$GeneID

#Convert counts to DGEList object
y <- DGEList(gene.counts)
# have a look at y
y
# See what slots are stored in y
names(y)
# Library size information is stored in the samples slot
y$samples

#Add groups for the samples
#We can also store the groups for the samples in the DGEList object.

group <- paste(sampleTable$CellType)
# Take a look
group

# Convert to factor
group <- factor(group)

# Add the group information into the DGEList
y$samples$group <- group
y$samples

#Adding annotation - Take care (mouse and Human are different)
columns(org.Hs.eg.db)

#build up our annotation information in a separate data frame using the select function
ann <- AnnotationDbi::select(org.Hs.eg.db,keys=rownames(y$counts),columns=c("ENTREZID","SYMBOL","GENENAME"))

#double check that the ENTREZID column matches exactly to our y$counts rownames.
table(ann$ENTREZID==rownames(y$counts))

#merge data witn the annotation (ann)
y$genes <- ann


################################
#Filtering lowly expressed genes#
###############################

#Genes with very low counts across all libraries provide little evidence for differential expression and they interfere with some of the 
#statistical approximations that are used later in the pipeline. They also add to the multiple testing burden when estimating false discovery rates, 
#reducing power to detect differentially expressed genes. These genes should be filtered out prior to further analysis.

#counts-per-million (CPM) above 0.5 in at least two samples
# Obtain CPMs
myCPM <- cpm(gene.counts)
# Have a look at the output
head(myCPM)

# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

# Summary of how many TRUEs there are in each row
# There are 11433 genes that have TRUEs in all 12 samples.
table(rowSums(thresh))

# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
summary(keep)

# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- gene.counts[keep,]
dim(gene.counts)

dim(counts.keep)


#A CPM of 0.5 is used as it corresponds to a count of 10-15 for the library sizes in this data set. If the count is any smaller, it is considered to 
#be very low, indicating that the associated gene is not expressed in that sample. A requirement for expression in two or more libraries is used as 
#each group contains two replicates. This ensures that a gene will be retained if it is only expressed in one group. Smaller CPM thresholds are 
#usually appropriate for larger libraries. As a general rule, a good threshold can be chosen by identifying the CPM that corresponds to a count of 10, 
#which in this case is about 0.5. You should filter with CPMs rather than filtering on the counts directly, as the latter does not account for
#differences in library sizes between samples.


# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
plot(myCPM[,1],gene.counts[,1])

# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],gene.counts[,1],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)

#filter the DGEList object
y <- y[keep, keep.lib.sizes=FALSE]

#Save csv file
write.csv(y, file="./results/counts_per_sample.csv")

##################
#Quality control#
#################
#Library size and distribution plots

y$samples$lib.size

## The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(y$samples$lib.size,names=colnames(y),las=2)
# Add a title to the plot
title("Barplot of library sizes")

# we can also adjust the labelling if we want
barplot(y$samples$lib.size/1e06, names=colnames(y), las=2, ann=FALSE, cex.names=0.75)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Library size (millions)", line = 3)
title("Barplot of library sizes")


#Count data is not normally distributed, so if we want to examine the distributions of the raw counts we need to log the counts. 
#Next we’ll use box plots to check the distribution of the read counts on the log2 scale. We can use the cpm function to get log2 counts per million, 
#which are corrected for the different library sizes. The cpm function also adds a small offset to avoid taking log of zero.

# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

#From the boxplots we see that overall the density distributions of raw log-intensities are not identical but still not very different. 
#If a sample is really far above or below the blue horizontal line we may need to investigate that sample further.

###############################
#Multidimensional scaling plots
###############################

plotMDS(y)

# We specify the option to let us plot two plots side-by-sde
par(mfrow=c(1,1))
# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
sampleTable$CellType <- as.factor(sampleTable$CellType)
levels(sampleTable$CellType)

## Let's choose purple for basal and orange for luminal
col.cell <- c("purple","orange", "blue", "red")[sampleTable$CellType]
data.frame(sampleTable$CellType,col.cell)

# Redo the MDS with cell type colouring
plotMDS(y,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend("topleft",fill=c("purple","orange", "blue", "red"),legend=levels(sampleTable$CellType))
# Add a title
title("HUVEC IL6 treatment")

#Another alternative is to generate an interactive MDS plot using the Glimma package. This allows the user to interactively explore the different 
#dimensions.
library("Glimma")
labels <- paste(sampleTable$SampleName, sampleTable$CellType)
glMDSPlot(y, labels=labels, groups=group, folder="mds")


#Hierarchical clustering with heatmaps
#An alternative to plotMDS for examining relationships between samples is using hierarchical clustering. Heatmaps are a nice visualisation to examine 
#hierarchical clustering of your samples. We can do this using the heatmap.2 function from the gplots package. In this example heatmap.2 calculates a 
#matrix of euclidean distances from the logCPM (logcounts object) for the 500 most variable genes. (Note this has more complicated code than plotting 
#principle components using plotMDS.)

# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)

# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)

# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange", "blue", "red")[sampleTable$CellType]

# Plot the heatmap
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row", srtCol=45)

# Save the heatmap
png(file="E:/R/Human_rnaseq/Results/High_var_genes.heatmap.png")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row", srtCol=45)
dev.off()

#Normalisation for composition bias
# Apply normalisation to DGEList object
y <- calcNormFactors(y)

y$samples

#########################################
#Differential expression with limma-voom#
#########################################

# Specify a design matrix without an intercept term
design <- model.matrix(~ 0 + group)
design
## Make the column names of the design matrix a bit nicer
colnames(design) <- levels(group)
design

#Each column of the design matrix tells us which samples correspond to each group. The samples which come from basal cells from a lactating mouse 
#correspond to columns 5 and 6 in the counts matrix, i.e. the samples which have 1s.

#Voom transform the data
par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)

# What is contained in this object?
names(v)

#repeat the box plots for the normalised data to compare to before normalisation. 
#The expression values in v$E are already log2 values so we don’t need to log-transform.

par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")

#####################################
#Testing for differential expression
#####################################
#Now that we have the voom transformed data we can use limma to test for differential expression. First we fit a linear model for each gene using the 
#lmFit function in limma. lmFit needs the voom object and the design matrix that we have already specified, which is stored within the voom object.

# Fit the linear model
fit <- lmFit(v)
names(fit)

#lmFit estimates group means according to the design matrix, as well as gene-wise variances. There are a number of items stored in the fit object, 
#most of which are specific to the statistical testing, and we won’t be discussing these in detail today.

cont.matrix <- makeContrasts(pbs72Vsil672=IL6_72h - PBS_72h, pbs96Vsil696=IL6_96h - PBS_96h , IL672VsIL696=IL6_96h - IL6_72h  ,levels=design)
#Take a look at the contrast matrix. The contrast matrix tells limma which columns of the design matrix we are interested in testing our comparison. 
#Note that here we have specified only one comparison to test, but we can specify as many as we want in one go.
cont.matrix

#apply the contrasts matrix to the fit object to get the statistics and estimated parameters of our comparison that we are interested in. 
#Here we call the contrasts.fit function in limma.

fit.cont <- contrasts.fit(fit, cont.matrix)

# final step is to call the eBayes function, which performs empirical Bayes shrinkage on the variances, and estimates moderated t-statistics and the 
#associated p-values.

fit.cont <- eBayes(fit.cont)
dim(fit.cont)

#We can use the limma decideTests function to generate a quick summary of DE genes for the contrasts
summa.fit <- decideTests(fit.cont)
summary(summa.fit)

#If two comparations - Sort.by = P
top.table <- topTable(fit.cont, sort.by = "F", n = Inf)
head(top.table, 20)

length(which(top.table$adj.P.Val < 0.05))
write.csv(top.table,file = "./results/DEG_Il6_PBS_72h_96h.csv")

###########################
#Plots after testing for DE
###########################

# We want to highlight the significant genes. We can get this from decideTests.
par(mfrow=c(1,2))
plotMD(fit.cont,coef=1,status=summa.fit[,"pbs96Vsil696"], values = c(-1, 1), hl.col=c("blue","red"))

# For the volcano plot we have to specify how many of the top genes to highlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
volcanoplot(fit.cont,coef=1,highlight=100,names=fit.cont$genes$SYMBOL, main="pbs96Vsil696")

#interactive version of the volcano plot
group2 <- group
levels(group2) <- c("IL6_72h  ","IL6_96h ", "PBS_72h", "PBS_96h")
glXYPlot(x=fit.cont$coefficients[,1], y=fit.cont$lods[,1],
         xlab="logFC", ylab="B", main="pbs72Vsil672",
         counts=v$E, groups=group2, status=summa.fit[,1],
         anno=fit.cont$genes, side.main="ENTREZID", folder="volcano")

#########################################
#Testing relative to a threshold (TREAT)#
#########################################


#################################
#Gene ontology testing with goana
#################################
go <- goana(fit.cont, coef="pbs72Vsil672",species = "Hs")
topGO(go, n=10)

write.csv(go,file = "./results/GO_pbs_vs_IL6_72hrs.csv")

#The row names of the output are the universal identifiers of the GO terms, with one term per row. The Term column gives the names of the GO terms. 
#These terms cover three domains - biological process (BP), cellular component (CC) and molecular function (MF), as shown in the Ont column. The N 
#column represents the total number of genes that are annotated with each GO term. The Up and Down columns represent the number of differentially 
#expressed genes that overlap with the genes in the GO term. The P.Up and P.Down columns contain the p-values for over-representation of the GO term 
#across the set of up- and down-regulated genes, respectively. The output table is sorted by the minimum of P.Up and P.Down by default.

#Gene Set Enrichment Analysis
library("msigdbr")
Hs.GOBP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
Hs.Reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
Hs.Hallmark <- msigdbr(species = "Homo sapiens", category = "H")

Hs.GOBP.Entrez <- Hs.GOBP %>% split(x = .$entrez_gene, f = .$gs_name)
Hs.Hallmark.Entrez <- Hs.Hallmark %>% split(x = .$entrez_gene, f = .$gs_name)
Hs.GOBP.Symbol <- Hs.GOBP %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.Hallmark.Symbol <- Hs.Hallmark %>% split(x = .$gene_symbol, f = .$gs_name)

#Camera competitive enrichment test

#create an index with the corresponding gene set database we generated.
idx <- ids2indices(Hs.GOBP.Entrez,id=v$genes$ENTREZID)

#perform camera enrichment - specifying the condition2 vs condition1 contrast, which is contrast 1 in our matrix
cam.HsGO <- camera(v,idx,design,contrast=cont.matrix[,1])
#Top 40 results
head(cam.HsGO,10)

write.csv(cam.HsGO,file = "./0_Figures fromRNAseq_manipulation_081722/r_results/GO_results_Camera competitive enrichment test.csv")

#create a barcodeplot of the enrichment of the geneset, plotting the t-value of each gene in the condition2vscondition1 contrast within this geneset
barcodeplot(fit.cont$t[,1], 
            index=idx$GOBP_CELL_CELL_JUNCTION_ORGANIZATION, 
            main="condition2vscondition1 - GOBP_CELL_CELL_JUNCTION_ORGANIZATION")


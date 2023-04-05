source("http://rnbeads.org/data/install.R")

library(RnBeads)
library(RnBeads.mm10)
library(RnBeads.mm10)
library(dplyr)
library(writexl)
library(foreach)
library(annotables)
library(tidyverse)
library(grid)

setwd("~/Analises R/r_array_mouse_013122/knockout_vs_wild_lps")
save.image(file = "knockout_vs_control.RData")
load(file = "knockout_vs_het_lps.RData")

save.rnb.set(rnb.set, path="~/Analises R/r_array_mouse_013122/Knockout_vs_het_lps/rnb_set", archive=F)

#Code for cellular deconvolution using methylation data-------
rnb.path <- "Analysis/reports_preprocessing/rnbSet_preprocessed/"#Specify the path to your preprocessed RnBSet
rnb.set <- load.rnb.set(rnb.path)

# Directory where your data is located
data.dir <- "~/Analises R/r_array_mouse_013122/knockout_vs_wild_lps"
idat.dir <- file.path(data.dir, "idat_files")
sample.annotation <- file.path(data.dir, "sample_annotation.csv")
# Directory where the output should be written to
analysis.dir <- "~/Analises R/r_array_mouse_013122/knockout_vs_wild_lps"
# Directory where the report files should be written to
report.dir <- file.path(analysis.dir, "reports")

show(sample.annotation)

#
rnb.options(assembly = "mm10", filtering.sex.chromosomes.removal=TRUE, filtering.snp = "no", qc.cnv=TRUE, 
            enforce.memory.management = TRUE, identifiers.column="Sample_ID",  export.to.csv = TRUE)

#
report.dir <- file.path(analysis.dir, "reports_details")
rnb.initialize.reports(report.dir)
logger.start(fname=NA)

#
data.source <- c(idat.dir, sample.annotation)
data.source

#Data Import
result <- rnb.run.import(data.source= data.source ,
                         data.type="idat.dir", dir.reports=report.dir)
rnb.set <- result$rnb.set

rnb.set
summary(pheno(rnb.set))
dim(M(rnb.set))
summary(M(rnb.set))
dim(meth(rnb.set))
summary(meth(rnb.set))
dim(meth(rnb.set, type="promoters"))
summary(meth(rnb.set, type="promoters"))
summary(dpval(rnb.set))
str(qc(rnb.set))
samples(rnb.set)
## # Visualize the bimodal distribution of beta values in sample 5
hist(meth(rnb.set)[, 5])

#QC quality control 
rnb.run.qc(rnb.set, report.dir)

#Pre processing sample standard 
rnb.set.unfiltered <- rnb.set
result <- rnb.run.preprocessing(rnb.set.unfiltered, dir.reports=report.dir)
save.rnb.set(rnb.set, path="~/Analises R/r_array_mouse_013122/knockout_vs_wild_lps/data_analysis/rnb.set.norm", archive=F)
rnb.set <- result$rnb.set
nrow(meth(rnb.set)) # the number of sites in the unfiltered object

#Normalization 
rnb.set.norm <- rnb.execute.normalization(rnb.set, method="bmiq",
                                          bgcorr.method="methylumi.noob")
save.rnb.set(rnb.set, path="~/Analises R/r_array_mouse_013122/knockout_vs_wild_lps/data_analysis/rnb.set.norm_2", archive=F)

nrow(meth(rnb.set.norm))

#running the exploratory analysis module
rnb.run.exploratory(rnb.set.norm, report.dir)
save.rnb.set(rnb.set, path="~/Analises R/r_array_mouse_013122/knockout_vs_wild_lps/data_analysis/rnb.set.expl", archive=F)
#Differential methylation analysis in \cs{RnBeads} can be conducted on the site level as well as on the genomic region level.
rnb.run.differential(rnb.set.norm, report.dir)

#Differential methylation
comp_cols <- colnames("Sample_Group")
# Specify the region types
reg_types <- c("genes", "promoters")
# Conduct the analysis
diffmeth <- rnb.execute.computeDiffMeth(rnb.set.norm, pheno.cols = comp_cols, region.types = reg_types)

####################################
#Save beta

annotations <- pheno(rnb.set.norm)
bvalues <- meth(rnb.set.norm, row.names = TRUE)
colnames(bvalues) <- annotations$Sample_ID
bvalues <- data.frame(CpG = rownames(bvalues), bvalues)
x <- list(bvalues, annotations)
names(x) <- c("bvalues", "annotations")
unlink(fileNameOut2)
write_xlsx(x, path = "~/Analises R/r_array_mouse_013122/knockout_vs_wild_lps/bvalue_norm.xlsx")

#Anottation package
# Prepate gene annotations, remove non-canonical chromosome names
#https://github.com/stephenturner/annotables
gene_annotations <- grcm38[ !(grepl("_", grcm38$chr) | grepl("GL", grcm38$chr)), c("ensgene", "symbol", "biotype", "description")]
gene_annotations <- gene_annotations[ !duplicated(gene_annotations) & !is.na(gene_annotations$symbol) & gene_annotations$description != "", ]

##############
#Annotation tables

comparison <- get.comparisons(diffmeth)[2]

diff.meth.proms <- get.table(diffmeth, comparison, "promoters", return.data.frame = TRUE)
diff.meth.genes <- get.table(diffmeth, comparison, "genes", return.data.frame = TRUE)
diff.meth.sites <- get.table(diffmeth, comparison, "sites", return.data.frame = TRUE)

sites.annot     <- annotation(rnb.set.norm, type="sites")
promoter.annot  <- annotation(rnb.set.norm, type="promoters")
genes.annot     <- annotation(rnb.set.norm, type ="genes")

annotated.diff.meth.sites <- data.frame(sites.annot, diff.meth.sites)
annotated.diff.meth.proms <- data.frame(promoter.annot, diff.meth.proms)
annotated.diff.meth.genes <- data.frame(genes.annot, diff.meth.genes) 

#Look how many rows have p < 0.05
p_adj_cutoff <- 0.05
nrow(annotated.diff.meth.sites[annotated.diff.meth.sites$diffmeth.p.val < p_adj_cutoff, ])
nrow(annotated.diff.meth.proms[annotated.diff.meth.proms$comb.p.val < p_adj_cutoff, ])
nrow(annotated.diff.meth.genes[annotated.diff.meth.genes$comb.p.val < p_adj_cutoff, ])

# All significant sites
de.sites <- annotated.diff.meth.sites[annotated.diff.meth.sites$diffmeth.p.val < p_adj_cutoff, ]
de.sites <- data.frame(CpG = rownames(de.sites), de.sites)
de.sites <- de.sites[order(de.sites$mean.diff), ]

#Illumina annotation csv file MouseMethylation-12v1-0_A1_Annotation_Mus_musculus.csv
illumina.annotation <- read.csv(file="./MouseMethylation_Annotation_Mus_musculus.csv")

de.sites.annot <- left_join(de.sites, illumina.annotation, by = c("Name"))
de.sites.annot <- de.sites.annot[order(de.sites.annot$mean.diff), ]
de.sites.annot_nodupli <- distinct(de.sites.annot,CpG, .keep_all= TRUE) #Remove duplicate values based on CpG

write_xlsx(de.sites.annot, "Differential_Methylation_cpg.xlsx" )

# All significant promoters
de.promoters <- annotated.diff.meth.proms[annotated.diff.meth.proms$comb.p.val < p_adj_cutoff, ]
# Append annotations
de.promoters <- data.frame(ensgene = rownames(de.promoters), de.promoters)
de.promoters <- left_join(de.promoters, gene_annotations, by = c("ensgene"))
# Order by most significant
de.promoters <- de.promoters[order(de.promoters$comb.p.val), ]
# Select protein-coding only
#de.promoters.coding <- de.promoters[de.promoters$biotype == "protein_coding" & !is.na(de.promoters$biotype), ]
# Select protein-coding upregulated, sorted by largest difference
de.promoters.up <- de.promoters[de.promoters$mean.mean.diff > 0, ]
de.promoters.up <- de.promoters.up[order(abs(de.promoters.up$mean.mean.diff), decreasing = TRUE), ]
# Select protein-coding upregulated, sorted by largest difference
de.promoters.dn <- de.promoters[de.promoters$mean.mean.diff < 0, ]
de.promoters.dn <- de.promoters.dn[order(abs(de.promoters.dn$mean.mean.diff), decreasing = TRUE), ]

# All significant DEGs 
de.genes <- annotated.diff.meth.genes[annotated.diff.meth.genes$comb.p.val < p_adj_cutoff, ]
# Append annotations
de.genes <- data.frame(ensgene = rownames(de.genes), de.genes)
de.genes <- left_join(de.genes, gene_annotations, by = c("ensgene"))
# Order by most significant
de.genes <- de.genes[order(de.genes$comb.p.val), ]
# Select protein-coding only
#de.genes.coding <- de.genes[de.genes$biotype == "protein_coding" & !is.na(de.genes$biotype), ]
# Select protein-coding upregulated, sorted by largest difference
de.genes.up <- de.genes[de.genes$mean.mean.diff > 0, ]
de.genes.up <- de.genes.up[order(abs(de.genes.up$mean.mean.diff), decreasing = TRUE), ]
# Select protein-coding upregulated, sorted by largest difference
de.genes.dn <- de.genes[de.genes$mean.mean.diff < 0, ]
de.genes.dn <- de.genes.dn[order(abs(de.genes.dn$mean.mean.diff), decreasing = TRUE), ]

# Save all results
# Note the full results make ~200Mb xlsx
x <- list(de.genes, de.genes.up, de.genes.dn,
          de.promoters, de.promoters.up, de.promoters.dn,
          de.sites)
names(x) <- c("DEGs", "DEGs.up", "DEGs.dn",
              "Promoters", "Promoters.up", "Promoters.dn",
              "CpGs")
fileNameOut1 = file.path(data.dir, "Differential_Methylation.xlsx")


unlink(fileNameOut1)
write_xlsx(x, path = fileNameOut1)

# Save methylation beta values
annotations <- pheno(rnb.set.norm)
bvalues <- meth(rnb.set.norm, row.names = TRUE)
colnames(bvalues) <- annotations$Sample_ID
bvalues <- data.frame(CpG = rownames(bvalues), bvalues)
x <- list(bvalues, annotations)
names(x) <- c("bvalues", "annotations")
fileNameOut2 = file.path(data.dir, "Methylation_bvalues.xlsx")

unlink(fileNameOut2)
write_xlsx(x, path = fileNameOut2)

#Load excel file if necessary

# Ramon's path
data.dir           = "/Users/mdozmorov/Documents/Data/VCU_work/Lathika"
fileNameIn1        = file.path(data.dir, "Differential_Methylation.xlsx")
fileSheet          = "Promoters.coding.dn" # "DEGs.coding"
fileNameOut1       = file.path(data.dir, paste0("Enrichment_", fileSheet, ".xlsx") )

mtx <- read_xlsx(fileNameIn1, sheet = fileSheet)
#Selecionar colunas específicas
res <- data.frame(symbol = mtx$symbol.y, logFC = mtx$mean.mean.diff, p.adj = mtx$comb.p.adj.fdr)




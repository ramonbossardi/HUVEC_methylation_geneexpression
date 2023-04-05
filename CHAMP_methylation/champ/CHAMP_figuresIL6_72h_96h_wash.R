# EB2022 Figures 
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(ggrepel)
library(Vennerable)
#Volcano Plot f
setwd("~/Analises R/r_array_huvec_120121/all_comparison_72hrs_96hrs_champ/EB2022/03_12_22")

raw_72hrs <- read.csv(file="./raw_result_72hrs.csv")
champ_72hrs <- read.csv(file="./DMP_IL6_72hrs_031222.csv")

champ_72hrs_volc <- champ_72hrs


save.image(file = "huvec_mice_ebfigures_031322.RData")
load(file = "huvec_mice_ebfigures_031322.RData")

#List with position p.adj <0.05

#### Select row with adj P value < 0.05 ###########
champ_72hrs_list <- filter(champ_72hrs, adj.P.Val <0.05 )
#dmp1 <- subset(dmp1, select = -diffexpressed)
champ_72hrs_list$diffexpressed <- "Non-difference CpGs"
# if log2Foldchange > 0.15 and pvalue < 0.001, set as "UP" 
champ_72hrs_list$diffexpressed[champ_72hrs_list$deltaBeta > 0 & champ_72hrs_list$adj.P.Val < 0.05] <- "Hypermethylated CpGs"
# if log2Foldchange < 0.15 and pvalue < 0.001, set as "DOWN"
champ_72hrs_list$diffexpressed[champ_72hrs_list$deltaBeta < -0 & champ_72hrs_list$adj.P.Val < 0.05] <- "Hypomethylated CpGs"



####### Start Volcano 72hrs Il6######
# add a column of NAs
#dmp1 <- subset(dmp1, select = -diffexpressed)
champ_72hrs_volc$diffexpressed <- "Non-difference CpGs"
# if log2Foldchange > 0.15 and pvalue < 0.001, set as "UP" 
champ_72hrs_volc$diffexpressed[champ_72hrs_volc$deltaBeta > 0 & champ_72hrs_volc$adj.P.Val < 0.05] <- "Hypermethylated CpGs"
# if log2Foldchange < 0.15 and pvalue < 0.001, set as "DOWN"
champ_72hrs_volc$diffexpressed[champ_72hrs_volc$deltaBeta < -0 & champ_72hrs_volc$adj.P.Val < 0.05] <- "Hypomethylated CpGs"

## Change point color 
mycolors <- c("blue", "black", "red")
names(mycolors) <- c("Hypomethylated CpGs", "Hypermethylated CpGs", "Non-difference CpGs")

#Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
champ_72hrs_volc$delabel <- NA
champ_72hrs_volc$delabel[champ_72hrs_volc$diffexpressed != "Non-difference CpGs"] <- champ_72hrs_volc$gene[champ_72hrs_volc$diffexpressed != "Non-difference CpGs"]

#Remove columns if necessary
#set1_72hrs_il <- subset(set1_72hrs_il, select=-c(diffexpressed, delabel))

champ_72hrs_volc_1<- filter(champ_72hrs_volc, diffexpressed %in% c("Hypermethylated CpGs", "Hypomethylated CpGs"))

# Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
# load library
#https://ggrepel.slowkow.com/articles/examples.html
library(ggrepel)
# plot adding up all layers we have seen so far
plot1 <- ggplot(data=champ_72hrs_volc_1, aes(x=deltaBeta, y=-log10(adj.P.Val), col=diffexpressed, label=delabel)) +
  geom_point() + 
  #geom_label_repel(min.segment.length = 0, seed = 32, box.padding = 0.5, max.overlaps = 30, segment.curvature = -0.1,
  #                 segment.ncp = 3,
  #                 segment.angle = 20) + 
  theme_minimal() +
  xlab("Delta Beta") + ylab("log10 adjusted p value") +
  #geom_text_repel() +
  scale_color_manual(values=c("blue", "red", "black")) +
  #geom_vline(xintercept=c(-0.05, 0.05), col="red") +
  #geom_hline(yintercept=-log10(0.05), col="red") +
  labs(title = "Control vs. IL6 72hrs HUVEC")
plot1

####### Heat map 72hrs ############
# Just to change the SampleID name, don't need to repeat
write.csv(champ_72hrs_list,file="./champ_72hrs_list.csv")
rm(champ_72hrs_heat, champ_72hrs_heatm)
champ_72hrs_list <- read.csv(file="./champ_72hrs_list.csv")

#Merge the raw data with the champ list with different positions
champ_72hrs_heatm <- merge(x= raw_72hrs, y= champ_72hrs_list, by=c("SampleID"), all=T)

#Select the samples that will be in the heatmap
champ_72hrs_heatm<- filter(champ_72hrs_heatm, diffexpressed %in% c("Hypermethylated CpGs", "Hypomethylated CpGs"))
write.csv(champ_72hrs_heatm,file="./champ_72hrs_heatm.csv")
champ_72hrs_heatm <- read.csv(file="./champ_72hrs_heatm.csv")

# Select collumns that are necessary to the heatmap
sample_heat <- champ_72hrs_heatm [, c( 24, 29, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)]

# Pivot make the matrix to heatmap
sample_heat_pivot <- sample_heat %>% 
  mutate_if(is.numeric,as.character, is.factor, as.character)
sapply(sample_heat_pivot, class)

sample_heat_pivot <- pivot_longer(data = sample_heat_pivot,
                                  cols = -c(1,2), 
                                  names_to = "Sample", 
                                  values_to = "Beta")

# Class of variables, Beta need to be numeric
sapply(sample_heat_pivot, class)
sample_heat_pivot$Beta <- as.numeric(sample_heat_pivot$Beta)
sample_heat_pivot$deltaBeta <- as.numeric(sample_heat_pivot$deltaBeta)

#Order the values on the X
head(sample_heat)
level_order <- c('PBS_exp2_1', 'PBS_exp3_1', 'PBS_exp2_2', 'PBS_exp3_2', 'PBS_exp1_1', 'PBS_exp1_2', 'IL6_exp2_1', 
                 'IL6_exp3_1', 'IL6_exp2_2', 'IL6_exp3_2', 'IL6_exp1_1', 'IL6_exp1_2', 'IL6_exp1_3') #this vector might be useful for other plots/analyses

p1 <- ggplot(data = sample_heat_pivot, mapping = aes(x = factor(Sample, level= level_order), y = reorder(gene, deltaBeta), fill = Beta)) +
  geom_tile() +
  xlab(label = "Sample") + 
  #remove x and y axis labels
  labs(x="",y="")+
  theme_grey(base_size=8)+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+ 
  ggtitle("CpG control vs. 72hrs IL6") +
  scale_alpha( trans = "log" ) +
  scale_fill_gradient2(low = "#353436",
                       mid = "#f6f805",
                       high = "red",
                       guide = "colorbar") +
  coord_equal() + #square 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

p1

p1 + coord_flip()

###### Heatmap 72hrs IL6 + 96hrs wash
raw_72hrs_96hrs_wash <- read.csv(file="./raw_result_72hrs_il6_96hrs_wash.csv")

champ_72hrs_96wash_heatm <- merge(x= champ_72hrs_heatm, y= raw_72hrs_96hrs_wash, by=c("SampleID"))
write.csv(champ_72hrs_96wash_heatm,file="./champ_72hrs_96wash.csv")

# Select collumns that are necessary to the heatmap
sample_heat_96wash <- champ_72hrs_96wash_heatm [, c( 29, 24, 78:94)]

# Pivot make the matrix to heatmap
sample_heat_96wash_pivot <- sample_heat_96wash %>% 
  mutate_if(is.numeric,as.character, is.factor, as.character)
sapply(sample_heat_96wash_pivot, class)

sample_heat_96wash_pivot <- pivot_longer(data = sample_heat_96wash_pivot,
                                  cols = -c(1,2), 
                                  names_to = "Sample", 
                                  values_to = "Beta")

# Class of variables, Beta need to be numeric
sapply(sample_heat_96wash_pivot, class)
sample_heat_96wash_pivot$Beta <- as.numeric(sample_heat_96wash_pivot$Beta)
sample_heat_96wash_pivot$deltaBeta <- as.numeric(sample_heat_96wash_pivot$deltaBeta)

head(sample_heat_96wash)
#Order the values on the X
head(sample_heat)
level_order_96wash <- c('hPBS_exp2_1_72', 'hPBS_exp3_1_72', 'hPBS_exp2_2_72', 'hPBS_exp3_2_72', 'hPBS_exp1_1_72', 'hPBS_exp1_2_72',
                        'hIL6_exp2_1_72', 'hIL6_exp3_1_72', 'hIL6_exp2_2_72', 'hIL6_exp3_2_72', 'hIL6_exp1_1_72', 'hIL6_exp1_2_72', 'hIL6_exp1_3_72',
                        'hIL6_exp2_1_96', 'hIL6_exp3_1_96', 'hIL6_exp2_2_96', 'hIL6_exp3_2_96')
                        
p2 <- ggplot(data = sample_heat_96wash_pivot, mapping = aes(x = factor(Sample, level= level_order_96wash), y = reorder(gene, deltaBeta), fill = Beta)) +
  geom_tile() +
  xlab(label = "Sample") + 
  #remove x and y axis labels
  labs(x="",y="")+
  theme_grey(base_size=8)+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+ 
  ggtitle("CpG control vs. 72hrs IL6 vs. 96hrs wash") +
  scale_alpha( trans = "log" ) +
  scale_fill_gradient2(low = "#353436",
                       mid = "#f6f805",
                       high = "red",
                       guide = "colorbar") +
  coord_equal() + #square 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

p2

p2 + coord_flip()

save.image(file = "huvec_ebfigures_031322.RData")
load("huvec_ebfigures_031322.RData")
#Filter difference on delta 72hrs and 96hrs -5 + 5

diff_72hrs_96hrs <- read.csv(file="~/Analises R/r_array_huvec_120121/all_comparison_72hrs_96hrs_champ/EB2022/03_12_22/champ_72hrs_96wash.csv")

diff_filter_72hrs_96hrs <- diff_72hrs_96hrs %>% filter(between (average_76_96, -0.02, 0.02))
write.csv(diff_filter_72hrs_96hrs,file="./cpg_filter_02_72hrs_96hrs.csv")




#Filter Hypermethylation genes at 72hrs ttr

cpg_methy<- filter(set2_72hrs_il, diffexpressed %in% c("Hypermethylated CpGs", "Hypomethylated CpGs"))

## Heatmap IL6 72hrs + 96 wash just values diff minor then 0.2

diff_filter_72hrs_96hrs <- diff_72hrs_96hrs %>% filter(between (average_76_96, -0.02, 0.02))

diff_heat_72_96 <- diff_filter_72hrs_96hrs [, c( 2, 98)]


###### Heatmap 72hrs IL6 + 96hrs wash

diff_heat_72_96 <- merge(x= diff_heat_72_96, y= champ_72hrs_96wash_heatm, by=c("SampleID"))

# Select collumns that are necessary to the heatmap
diff_heat_72_96 <- diff_heat_72_96 [, c( 25, 30, 79:95)]

# Pivot make the matrix to heatmap
diff_heat_96wash_pivot <- diff_heat_72_96 %>% 
  mutate_if(is.numeric,as.character, is.factor, as.character)
sapply(diff_heat_96wash_pivot, class)

diff_heat_96wash_pivot <- pivot_longer(data = diff_heat_96wash_pivot,
                                         cols = -c(1,2), 
                                         names_to = "Sample", 
                                         values_to = "Beta")

# Class of variables, Beta need to be numeric
sapply(diff_heat_96wash_pivot, class)
diff_heat_96wash_pivot$Beta <- as.numeric(diff_heat_96wash_pivot$Beta)
diff_heat_96wash_pivot$deltaBeta <- as.numeric(diff_heat_96wash_pivot$deltaBeta)

head(diff_heat_72_96)
#Order the values on the X
head(sample_heat)
level_order_96wash <- c('hPBS_exp2_1_72', 'hPBS_exp3_1_72', 'hPBS_exp2_2_72', 'hPBS_exp3_2_72', 'hPBS_exp1_1_72', 'hPBS_exp1_2_72',
                        'hIL6_exp2_1_72', 'hIL6_exp3_1_72', 'hIL6_exp2_2_72', 'hIL6_exp3_2_72', 'hIL6_exp1_1_72', 'hIL6_exp1_2_72', 'hIL6_exp1_3_72',
                        'hIL6_exp2_1_96', 'hIL6_exp3_1_96', 'hIL6_exp2_2_96', 'hIL6_exp3_2_96')

p8 <- ggplot(data = diff_heat_96wash_pivot, mapping = aes(x = factor(Sample, level= level_order_96wash), y = reorder(gene, deltaBeta), fill = Beta)) +
  geom_tile() +
  xlab(label = "Sample") + 
  #remove x and y axis labels
  labs(x="",y="")+
  theme_grey(base_size=8)+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+ 
  ggtitle("CpG control vs. 72hrs IL6 vs. 96hrs wash") +
  scale_alpha( trans = "log" ) +
  scale_fill_gradient2(low = "blue",
                       mid = "white ",
                       high = "red",
                       guide = "colorbar") +
  coord_equal() + #square 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

p8

p8 + coord_flip()

###### Heatmap 72hrs IL6 + 96hrs wash with 96hrs control 

diff_72hrs_96hrs_control96 <- read.csv(file="~/Analises R/r_array_huvec_120121/all_comparison_72hrs_96hrs_champ/EB2022/03_12_22/champ_72hrs_96wash.csv")

diff_filter_72hrs_96hrs_c <- diff_72hrs_96hrs %>% filter(between (average_76_96, -0.02, 0.02))

diff_heat_72_96_c <- diff_filter_72hrs_96hrs_c [, c( 2, 98)]

###### Heatmap 72hrs IL6 + 96hrs wash

diff_heat_72_96_control96 <- merge(x= diff_heat_72_96_c, y= diff_72hrs_96hrs_control96, by=c("SampleID"))

# Select collumns that are necessary to the heatmap
diff_heat_72_96_control96 <- diff_heat_72_96_control96 [, c( 26, 31, 80:100)]

# Pivot make the matrix to heatmap
diff_heat_96wash_pivot_control96 <- diff_heat_72_96_control96 %>% 
  mutate_if(is.numeric,as.character, is.factor, as.character)
sapply(diff_heat_96wash_pivot_control96, class)

diff_heat_96wash_pivot_control96 <- pivot_longer(data = diff_heat_96wash_pivot_control96,
                                       cols = -c(1,2), 
                                       names_to = "Sample", 
                                       values_to = "Beta")

# Class of variables, Beta need to be numeric
sapply(diff_heat_96wash_pivot_control96, class)
diff_heat_96wash_pivot_control96$Beta <- as.numeric(diff_heat_96wash_pivot_control96$Beta)
diff_heat_96wash_pivot_control96$deltaBeta <- as.numeric(diff_heat_96wash_pivot_control96$deltaBeta)

head(diff_heat_72_96)
#Order the values on the X
head(sample_heat)
level_order_96wash_control96 <- c('hPBS_exp2_1_72', 'hPBS_exp3_1_72', 'hPBS_exp2_2_72', 'hPBS_exp3_2_72', 'hPBS_exp1_1_72', 'hPBS_exp1_2_72',
                                  'hPBS_exp2_1_96', 'hPBS_exp3_1_96', 'hPBS_exp2_2_96', 'hPBS_exp3_2_96',
                        'hIL6_exp2_1_72', 'hIL6_exp3_1_72', 'hIL6_exp2_2_72', 'hIL6_exp3_2_72', 'hIL6_exp1_1_72', 'hIL6_exp1_2_72', 'hIL6_exp1_3_72',
                        'hIL6_exp2_1_96', 'hIL6_exp3_1_96', 'hIL6_exp2_2_96', 'hIL6_exp3_2_96')

p9 <- ggplot(data = diff_heat_96wash_pivot_control96, mapping = aes(x = factor(Sample, level= level_order_96wash_control96), y = reorder(gene, deltaBeta), fill = Beta)) +
  geom_tile() +
  xlab(label = "Sample") + 
  #remove x and y axis labels
  labs(x="",y="")+
  theme_grey(base_size=8)+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+ 
  ggtitle("CpG control vs. 72hrs IL6 vs. 96hrs wash") +
  scale_alpha( trans = "log" ) +
  scale_fill_gradient2(low = "blue",
                       mid = "white ",
                       high = "red",
                       guide = "colorbar") +
  coord_equal() + #square 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

p9

p9 + coord_flip()


##############Animal results #####################################

#Salive vs. LPS wild
fileNameIn = file.path("~/Analises R/r_array_huvec_120121/all_comparison_72hrs_96hrs_champ/EB2022/03_12_22/mice", "DMP_saline_lpswild.xlsx")
mice_saline_wildlps <- read_xlsx(fileNameIn)
#Remove duplicate values based on CpG
mice_saline_wildlps <- distinct(mice_saline_wildlps,CpG, .keep_all= TRUE) 

#load B norm values
fileNameIn1 = file.path("~/Analises R/r_array_huvec_120121/all_comparison_72hrs_96hrs_champ/EB2022/03_12_22/mice", "bvalue__saline_lpswild.xlsx")
mice_b_saline_wildlps <- read_xlsx(fileNameIn1)
#Select p-value <0.05
mice_saline_wildlps_heat <- filter(mice_b_saline_wildlps, diffmeth.p.val <0.05 )
write.csv(mice_saline_wildlps_heat,file="~/Analises R/r_array_huvec_120121/all_comparison_72hrs_96hrs_champ/EB2022/03_12_22/mice/saline_lpswild_heat.csv")

#Merge the raw data with the champ list with different positions
mice_b_saline_wildlps <- merge(x= mice_b_saline_wildlps, y= mice_saline_wildlps, by=c("CpG"), all=T)

#mice_b_saline_wildlps <- subset(mice_b_saline_wildlps, select = -diffexpressed)
mice_b_saline_wildlps$diffexpressed <- "Non-difference CpGs"
# if log2Foldchange > 0.15 and pvalue < 0.001, set as "UP" 
mice_b_saline_wildlps$diffexpressed[mice_b_saline_wildlps$mean.diff > 0 & mice_b_saline_wildlps$diffmeth.p.val < 0.05] <- "Hypermethylated CpGs"
# if log2Foldchange < 0.15 and pvalue < 0.001, set as "DOWN"
mice_b_saline_wildlps$diffexpressed[mice_b_saline_wildlps$mean.diff < -0 & mice_b_saline_wildlps$diffmeth.p.val < 0.05] <- "Hypomethylated CpGs"

#Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
mice_b_saline_wildlps$delabel <- NA
mice_b_saline_wildlps$delabel[mice_b_saline_wildlps$diffexpressed != "Non-difference CpGs"] <- mice_b_saline_wildlps$Gene[mice_b_saline_wildlps$diffexpressed != "Non-difference CpGs"]

#Remove columns if necessary
#set1_72hrs_il <- subset(set1_72hrs_il, select=-c(diffexpressed, delabel))

# Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
# load library
#https://ggrepel.slowkow.com/articles/examples.html
library(ggrepel)
# plot adding up all layers we have seen so far
mice_plot1 <- ggplot(data=mice_b_saline_wildlps, aes(x=mean.diff, y=-log10(diffmeth.p.val), col=diffexpressed, label=delabel)) +
  geom_point() + 
  #geom_label_repel(min.segment.length = 0, seed = 32, box.padding = 0.5, max.overlaps = 30, segment.curvature = -0.1,
  #                 segment.ncp = 3,
  #                 segment.angle = 20) + 
  theme_minimal() +
  xlab("Delta Beta") + ylab("log10 adjusted p value") +
  #geom_text_repel() +
  scale_color_manual(values=c("blue", "red", "black")) +
  #geom_vline(xintercept=c(-0.05, 0.05), col="red") +
  #geom_hline(yintercept=-log10(0.05), col="red") +
  labs(title = "Saline vs. LPS Wild Type")
mice_plot1

# LPS wild  vs. LPS FL/FL
fileNameIn2 = file.path("~/Analises R/r_array_huvec_120121/all_comparison_72hrs_96hrs_champ/EB2022/03_12_22/mice", "DMP_lpswild_knockout.xlsx")
mice_wildlps_knockout <- read_xlsx(fileNameIn2)
#Remove duplicate values based on CpG
mice_wildlps_knockout <- distinct(mice_wildlps_knockout,CpG, .keep_all= TRUE) 

#load B norm values
fileNameIn3 = file.path("~/Analises R/r_array_huvec_120121/all_comparison_72hrs_96hrs_champ/EB2022/03_12_22/mice", "bvalue_norm_lps_wild_knockout.xlsx")
mice_b_wildlps_knockout <- read_xlsx(fileNameIn3)

#Merge the raw data with the champ list with different positions
mice_b_wildlps_knockout <- merge(x= mice_b_wildlps_knockout, y= mice_wildlps_knockout, by=c("CpG"), all=T)

#Select p-value <0.05
mice_b_wildlps_knockout_heat <- filter(mice_b_wildlps_knockout, diffmeth.p.val <0.05 )
write.csv(mice_b_wildlps_knockout_heat,file="~/Analises R/r_array_huvec_120121/all_comparison_72hrs_96hrs_champ/EB2022/03_12_22/mice/lpswild_knockout_heat.csv")

#mice_b_wildlps_knockout <- subset(mice_b_wildlps_knockout, select = -diffexpressed)
mice_b_wildlps_knockout$diffexpressed <- "Non-difference CpGs"
# if log2Foldchange > 0.15 and pvalue < 0.001, set as "UP" 
mice_b_wildlps_knockout$diffexpressed[mice_b_wildlps_knockout$mean.diff > 0 & mice_b_wildlps_knockout$diffmeth.p.val < 0.05] <- "Hypermethylated CpGs"
# if log2Foldchange < 0.15 and pvalue < 0.001, set as "DOWN"
mice_b_wildlps_knockout$diffexpressed[mice_b_wildlps_knockout$mean.diff < -0 & mice_b_wildlps_knockout$diffmeth.p.val < 0.05] <- "Hypomethylated CpGs"

#Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
mice_b_wildlps_knockout$delabel <- NA
mice_b_wildlps_knockout$delabel[mice_b_wildlps_knockout$diffexpressed != "Non-difference CpGs"] <- mice_b_wildlps_knockout$Gene[mice_b_wildlps_knockout$diffexpressed != "Non-difference CpGs"]

#Remove columns if necessary
#set1_72hrs_il <- subset(set1_72hrs_il, select=-c(diffexpressed, delabel))

# Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
# load library
#https://ggrepel.slowkow.com/articles/examples.html
library(ggrepel)
# plot adding up all layers we have seen so far
mice_plot2 <- ggplot(data=mice_b_wildlps_knockout, aes(x=mean.diff, y=-log10(diffmeth.p.val), col=diffexpressed, label=delabel)) +
  geom_point() + 
  #geom_label_repel(min.segment.length = 0, seed = 32, box.padding = 0.5, max.overlaps = 30, segment.curvature = -0.1,
  #                 segment.ncp = 3,
  #                 segment.angle = 20) + 
  theme_minimal() +
  xlab("Delta Beta") + ylab("log10 adjusted p value") +
  #geom_text_repel() +
  scale_color_manual(values=c("blue", "red", "black")) +
  #geom_vline(xintercept=c(-0.05, 0.05), col="red") +
  #geom_hline(yintercept=-log10(0.05), col="red") +
  labs(title = "LPS Wild vs. Knockout")
mice_plot2

#Select p-value <0.05
mice_wildlps_knockout_heat <- filter(mice_b_wildlps_knockout, diffmeth.p.val <0.05 )
write.csv(mice_wildlps_knockout_heat,file="~/Analises R/r_array_huvec_120121/all_comparison_72hrs_96hrs_champ/EB2022/03_12_22/mice/lpswild_knockout_heat.csv")

mice_b_saline_wildlps_knockout <- merge(x= mice_b_saline_wildlps, y= mice_saline_wildlps, by=c("CpG"), all=T)

###### Venn Diagram saline/wildtype vs. wildtype/knockout ##############
library(Vennerable)

set1 <- mice_saline_wildlps_heat
set1 <- set1$CpG %>% unique()
set2 <- mice_wildlps_knockout_heat
set2 <- set2$CpG %>% unique()

Venn(list(set1, set2), SetNames = c("saline/wildtype", "LPS wildtype/knockout")) %>% plot

######Heatmap mice##############

mice_b_heat_saline_wildlps <- read.csv(file="~/Analises R/r_array_huvec_120121/all_comparison_72hrs_96hrs_champ/EB2022/03_12_22/mice/heatmap_cpg_031622.csv")

fileNameIn4 = file.path("~/Analises R/r_array_huvec_120121/all_comparison_72hrs_96hrs_champ/EB2022/03_12_22/mice", "bvalue_norm_saline_lpswild_lpsikeo.xlsx")
mice_rawdata_saline_wild_knockout <- read_xlsx(fileNameIn4)

mice_heat_1 <- merge(x= mice_b_heat_saline_wildlps, y= mice_rawdata_saline_wild_knockout, by=c("CpG"))

write.csv(mice_heat_1,file="~/Analises R/r_array_huvec_120121/all_comparison_72hrs_96hrs_champ/EB2022/03_12_22/mice/lpswild_knockout_heat.csv")


# Select collumns that are necessary to the heatmap
mice_heat_2 <- mice_heat_1 [, c( 4, 3, 21:35)]

# Pivot make the matrix to heatmap
mice_heat_pivot <- mice_heat_2 %>% 
  mutate_if(is.numeric,as.character, is.factor, as.character)
sapply(mice_heat_pivot, class)

mice_heat_pivot <- pivot_longer(data = mice_heat_pivot,
                                  cols = -c(1,2), 
                                  names_to = "Sample", 
                                  values_to = "Beta")

# Class of variables, Beta need to be numeric
sapply(mice_heat_pivot, class)
mice_heat_pivot$Beta <- as.numeric(mice_heat_pivot$Beta)
mice_heat_pivot$mean.diff <- as.numeric(mice_heat_pivot$mean.diff)

#Order the values on the X
head(mice_heat_2)
level_order <- c( 'hsaline_1',   'hsaline_2',   'hsaline_3',   'hsaline_4',   'hsaline_5',   'hsaline_6',  'hlpswild_1',  'hlpswild_2',  'hlpswild_3',  'hlpswild_4', 'hlpsiko_1...28', 'hlpsiko_1...29', 'hlpsiko_1...30', 'hlpsiko_1...31', 'hlpsiko_1...32') 

p3 <- ggplot(data = mice_heat_pivot, mapping = aes(x = factor(Sample, level= level_order), y = reorder(Gene, mean.diff), fill = Beta)) +
  geom_tile() +
  xlab(label = "Sample") + 
  #remove x and y axis labels
  labs(x="",y="")+
  theme_grey(base_size=8)+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+ 
  ggtitle("Control vs. LPS wildtype vs. LPSiKEO") +
  scale_alpha( trans = "log" ) +
  scale_fill_gradient2(low = "#353436",
                       mid = "#f6f805",
                       high = "blue",
                       guide = "colorbar") +
  coord_equal() + #square 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

p3

p3 + coord_flip()



#Trying to organize the heatmap by delta between groups
write.csv(mice_heat_1,file="./raw__72hrs_il6_96hrs_heatmap.csv")

mice_heat_3 <- read.csv(file="./raw__72hrs_il6_96hrs_heatmap.csv")


mice_heat_4 <- mice_heat_3 [, c(5, 22:36, 40)]

# Pivot make the matrix to heatmap
mice_heat_pivot_1 <- mice_heat_4 %>% 
  mutate_if(is.numeric,as.character, is.factor, as.character)
sapply(mice_heat_pivot_1, class)

mice_heat_pivot_1 <- pivot_longer(data = mice_heat_pivot_1,
                                cols = -c(1, 17), 
                                names_to = "Sample", 
                                values_to = "Beta")

# Class of variables, Beta need to be numeric
sapply(mice_heat_pivot_1, class)
mice_heat_pivot_1$Beta <- as.numeric(mice_heat_pivot_1$Beta)
mice_heat_pivot_1$delta_saline_wt <- as.numeric(mice_heat_pivot_1$delta_saline_wt)

#Order the values on the X
head(mice_heat_2)
level_order <- c( 'hsaline_1',   'hsaline_2',   'hsaline_3',   'hsaline_4',   'hsaline_5',   'hsaline_6',  'hlpswild_1',  'hlpswild_2',  'hlpswild_3',  'hlpswild_4', 'hlpsiko_1...28', 'hlpsiko_1...29', 'hlpsiko_1...30', 'hlpsiko_1...31', 'hlpsiko_1...32') 

p4 <- ggplot(data = mice_heat_pivot_1, mapping = aes(x = factor(Sample, level= level_order), y = reorder(Gene, delta_saline_wt), fill = Beta)) +
  geom_tile() +
  xlab(label = "Sample") + 
  #remove x and y axis labels
  labs(x="",y="")+
  theme_grey(base_size=8)+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+ 
  ggtitle("Control vs. LPS wildtype vs. LPSiKEO") +
  scale_alpha( trans = "log" ) +
  scale_fill_gradient2(low = "#353436",
                       mid = "#f6f805",
                       high = "blue",
                       guide = "colorbar") +
  coord_equal() + #square 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

p4

p4 + coord_flip()

##### Heatmap with CpGs that changed


set3 <- merge(x= set1, y= set2, by=c("CpG"))
write.csv(set3,file="./genes_same_salinevswtlps_lpssocs3ieko.csv")


set4 <- subset(set3, select = -c(2:126))

set4 <- merge(x= set4, y= mice_heat_3, by=c("CpG"))
set4 <- distinct(set4,CpG, .keep_all= TRUE) 

set5 <- set4 [, c(6, 23:37, 41)]

# Pivot make the matrix to heatmap
mice_heat_pivot_2 <- set5 %>% 
  mutate_if(is.numeric,as.character, is.factor, as.character)
sapply(mice_heat_pivot_2, class)

mice_heat_pivot_2 <- pivot_longer(data = mice_heat_pivot_2,
                                  cols = -c(1, 17), 
                                  names_to = "Sample", 
                                  values_to = "Beta")

# Class of variables, Beta need to be numeric
sapply(mice_heat_pivot_2, class)
mice_heat_pivot_2$Beta <- as.numeric(mice_heat_pivot_2$Beta)
mice_heat_pivot_2$delta_saline_wt <- as.numeric(mice_heat_pivot_2$delta_saline_wt)

#Order the values on the X
head(mice_heat_2)
level_order <- c( 'hsaline_1',   'hsaline_2',   'hsaline_3',   'hsaline_4',   'hsaline_5',   'hsaline_6',  'hlpswild_1',  'hlpswild_2',  'hlpswild_3',  'hlpswild_4', 'hlpsiko_1...28', 'hlpsiko_1...29', 'hlpsiko_1...30', 'hlpsiko_1...31', 'hlpsiko_1...32') 

p5 <- ggplot(data = mice_heat_pivot_2, mapping = aes(x = factor(Sample, level= level_order), y = reorder(Gene, delta_saline_wt), fill = Beta)) +
  geom_tile() +
  xlab(label = "Sample") + 
  #remove x and y axis labels
  labs(x="",y="")+
  theme_grey(base_size=8)+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+ 
  ggtitle("Control vs. LPS wildtype vs. LPSiKEO") +
  scale_alpha( trans = "log" ) +
  scale_fill_gradient2(low = "#353436",
                       mid = "#f6f805",
                       high = "blue",
                       guide = "colorbar") +
  coord_equal() + #square 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

p5

p5 + coord_flip()


#### Heat map just saline vs. wt LPS

set6<- filter(mice_b_saline_wildlps, diffexpressed %in% c("Hypermethylated CpGs", "Hypomethylated CpGs"))
write.csv(set6,file="./mice_b_saline_wildlp_sigcpg.csv")
set6 <- read.csv(file="./mice_b_saline_wildlp_sigcpg.csv")

set6 <- distinct(set6,CpG, .keep_all= TRUE) 

set6 <- set6 [, c(13, 16:25, 28)]

# Pivot make the matrix to heatmap
mice_heat_pivot_3 <- set6 %>% 
  mutate_if(is.numeric,as.character, is.factor, as.character)
sapply(mice_heat_pivot_3, class)

mice_heat_pivot_3 <- pivot_longer(data = mice_heat_pivot_3,
                                  cols = -c(1, 12), 
                                  names_to = "Sample", 
                                  values_to = "Beta")

# Class of variables, Beta need to be numeric
sapply(mice_heat_pivot_3, class)
mice_heat_pivot_3$Beta <- as.numeric(mice_heat_pivot_3$Beta)
mice_heat_pivot_3$diff_meth <- as.numeric(mice_heat_pivot_3$diff_meth)

#Order the values on the X
head(mice_heat_2)
level_order <- c( 'hsaline_1',   'hsaline_2',   'hsaline_3',   'hsaline_4',   'hsaline_5',   'hsaline_6',  'hlpswild_1',  'hlpswild_2',  'hlpswild_3',  'hlpswild_4', 'hlpsiko_1...28', 'hlpsiko_1...29', 'hlpsiko_1...30', 'hlpsiko_1...31', 'hlpsiko_1...32') 
level_order_2 <- c( 'hsaline_1',   'hsaline_2',   'hsaline_3',   'hsaline_4',   'hsaline_5',   'hsaline_6',  'hlpswild_1',  'hlpswild_2',  'hlpswild_3',  'hlpswild_4')
p6 <- ggplot(data = mice_heat_pivot_3, mapping = aes(x = factor(Sample, level= level_order_2), y = reorder(Gene, diff_meth), fill = Beta)) +
  geom_tile() +
  xlab(label = "Sample") + 
  #remove x and y axis labels
  labs(x="",y="")+
  theme_grey(base_size=8)+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+ 
  ggtitle("Control vs. LPS WT") +
  scale_alpha( trans = "log" ) +
  scale_fill_gradient2(low = "#353436",
                       mid = "#f6f805",
                       high = "red",
                       guide = "colorbar") +
  coord_equal() + #square 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

p6

p6 + coord_flip()

###########LPS WT vs. LPS SOCS3iKEO

######Heatmap mice##############

set7 <- read.csv(file="~/Analises R/r_array_huvec_120121/all_comparison_72hrs_96hrs_champ/EB2022/03_12_22/mice/lpswild_knockout_heat.csv")


# Select collumns that are necessary to the heatmap
mice_heat_5 <- set7 [, c( 30, 19, 37:45)]

# Pivot make the matrix to heatmap
mice_heat_pivot_5 <- mice_heat_5 %>% 
  mutate_if(is.numeric,as.character, is.factor, as.character)
sapply(mice_heat_pivot_5, class)

mice_heat_pivot_5 <- pivot_longer(data = mice_heat_pivot_5,
                                cols = -c(1,2), 
                                names_to = "Sample", 
                                values_to = "Beta")

# Class of variables, Beta need to be numeric
sapply(mice_heat_pivot_5, class)
mice_heat_pivot_5$Beta <- as.numeric(mice_heat_pivot_5$Beta)
mice_heat_pivot_5$mean.diff <- as.numeric(mice_heat_pivot_5$mean.diff)

#Order the values on the X
head(mice_heat_2)
level_order_1 <- c( 'hlpswild_1',   'hlpswild_2',   'hlpswild_3',   'hlpswild_4',   'hikeo_1',   'hikeo_2',  'hikeo_3',  'hikeo_4',  'hikeo_5') 
p7 <- ggplot(data = mice_heat_pivot_5, mapping = aes(x = factor(Sample, level= level_order_1), y = reorder(Gene, mean.diff), fill = Beta)) +
  geom_tile() +
  xlab(label = "Sample") + 
  #remove x and y axis labels
  labs(x="",y="")+
  theme_grey(base_size=8)+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+ 
  ggtitle("LPS wildtype vs. LPSiKEO") +
  scale_alpha( trans = "log" ) +
  scale_fill_gradient2(low = "#353436",
                       mid = "#white",
                       high = "red",
                       guide = "colorbar") +
  coord_equal() + #square 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

p7

p7 + coord_flip()

## VennDiargram HUVEC vs Animals

###### Venn Diagram saline/wildtype vs. wildtype/knockout ##############
library(Vennerable)

set0 <- sample_heat #CpG 431 found in HUVEC
set0 <- set0$gene %>% unique()
#CpGs changed animals
set0.0 <- read.csv(file="~/Analises R/r_array_huvec_120121/all_comparison_72hrs_96hrs_champ/EB2022/03_12_22/mice/genes_compare_huvec.csv")
#Change column gene name to make identifal between both dataframes
names(set0.0)[names(set0.0) == "Gene"] <- "gene"
set0.0 <- set0.0$gene %>% unique()

Venn(list(set0, set0.0), SetNames = c("HUVEC", "Animal")) %>% plot

#Cpgs identical mouse and huvec
huvec_mice_CpGs <- merge(x= set0, y= set0.0, by=c("gene"))
write.csv(huvec_mice_CpGs,file="./genes_huvec_mice.csv")

### Venn diagram Hypomethylated 
set1_hypo<- filter(set1, diffexpressed %in% c("Hypomethylated CpGs"))
set2_hypo<- filter(set2, diffexpressed %in% c("Hypomethylated CpGs"))

set1_hypo <- set1_hypo$Gene %>% unique()
set2_hypo <- set2_hypo$Gene %>% unique()

Venn(list(set1_hypo, set2_hypo), SetNames = c("LPS Wt", "LPS SOCS3iEKO")) %>% plot

### Venn diagram Hypermethylated 
set1_hyper<- filter(set1, diffexpressed %in% c("Hypermethylated CpGs"))
set2_hyper<- filter(set2, diffexpressed %in% c("Hypermethylated CpGs"))

set1_hyper <- set1_hyper$Gene %>% unique()
set2_hyper <- set2_hyper$Gene %>% unique()

Venn(list(set1_hyper, set2_hyper, set1_hypo, set2_hypo), SetNames = c("set1_hyper", "set2_hyper", "set1_hypo", "set2_hypo")) %>% plot


#Heat map gene expression EB presentation 

#Merge the raw data with the champ list with different positions
eb_champ_72hrs_heatm <- merge(x= raw_72hrs, y= champ_72hrs_list, by=c("SampleID"), all=T)

#Select the samples that will be in the heatmap
eb_champ_72hrs_heatm<- filter(eb_champ_72hrs_heatm, SampleID %in% c("cg12857064", "cg13284285", "cg09648467"))

write.csv(eb_champ_72hrs_heatm,file="./eb_champ_72hrs_heatm.csv")

eb_champ_72hrs_heatm <- read.csv(file="./eb_champ_72hrs_heatm.csv")

# Select collumns that are necessary to the heatmap
eb_sample_heat <- eb_champ_72hrs_heatm [, c( 16, 38:51)]

# Pivot make the matrix to heatmap
eb_sample_heat_pivot <- eb_sample_heat %>% 
  mutate_if(is.numeric,as.character, is.factor, as.character)
sapply(eb_sample_heat_pivot, class)

eb_sample_heat_pivot <- pivot_longer(data = eb_sample_heat_pivot,
                                  cols = -c(1,2), 
                                  names_to = "Sample", 
                                  values_to = "Beta")

# Class of variables, Beta need to be numeric
sapply(eb_sample_heat_pivot, class)
eb_sample_heat_pivot$Beta <- as.numeric(eb_sample_heat_pivot$Beta)
eb_sample_heat_pivot$mean_pbs <- as.numeric(eb_sample_heat_pivot$mean_pbs)

#Order the values on the X
head(sample_heat)
level_order_eb <- c('PBS_exp2_1.1', 'PBS_exp3_1.1', 'PBS_exp2_2.1', 'PBS_exp3_2.1', 'PBS_exp1_1.1', 'PBS_exp1_2.1', 'IL6_exp2_1.1',
                 'IL6_exp3_1.1', 'IL6_exp2_2.1', 'IL6_exp3_2.1', 'IL6_exp1_1.1', 'IL6_exp1_2.1', 'IL6_exp1_3.1') #this vector might be useful for other plots/analyses

p10 <- ggplot(data = eb_sample_heat_pivot, mapping = aes(x = factor(Sample, level= level_order_eb), y = gene, fill = Beta)) +
  geom_tile() +
  xlab(label = "Sample") + 
  #remove x and y axis labels
  #labs(x="",y="")+
  theme_grey(base_size=8)+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+ 
  #ggtitle("CpG control vs. 72hrs IL6") +
  scale_alpha( trans = "log" ) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       guide = "colorbar") +
  coord_equal() + #square 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

p10

p10 + coord_flip()


## SETDB1 gene methylation 

setdb1_gene <- diff_72hrs_96hrs %>% filter(gene %in% c("SETDB1"))

setdb1_gene <- setdb1_gene [, c( 2, 98)]
setdb1_gene <- merge(x= setdb1_gene, y= diff_72hrs_96hrs_control96, by=c("SampleID"))

# Select collumns that are necessary to the heatmap
setdb1_gene <- setdb1_gene [, c( 31, 79, 80:100)]

# Pivot make the matrix to heatmap
setdb1_gene_pivot <- setdb1_gene %>% 
  mutate_if(is.numeric,as.character, is.factor, as.character)
sapply(setdb1_gene_pivot, class)

setdb1_gene_pivot <- pivot_longer(data = setdb1_gene_pivot,
                                         cols = -c(1,2), 
                                         names_to = "Sample", 
                                         values_to = "Beta")

# Class of variables, Beta need to be numeric
sapply(setdb1_gene_pivot, class)
setdb1_gene_pivot$Beta <- as.numeric(setdb1_gene_pivot$Beta)
setdb1_gene_pivot$deltaBeta_72hrs <- as.numeric(setdb1_gene_pivot$deltaBeta_72hrs)

head(sample_heat_96wash)
#Order the values on the X
head(sample_heat)
level_order_setdb1 <- c('hPBS_exp2_1_72', 'hPBS_exp3_1_72', 'hPBS_exp2_2_72', 'hPBS_exp3_2_72', 'hPBS_exp1_1_72', 'hPBS_exp1_2_72',
                                  'hIL6_exp2_1_72', 'hIL6_exp3_1_72', 'hIL6_exp2_2_72', 'hIL6_exp3_2_72', 'hIL6_exp1_1_72', 'hIL6_exp1_2_72', 'hIL6_exp1_3_72')

p10 <- ggplot(data = setdb1_gene_pivot, mapping = aes(x = factor(Sample, level= level_order_setdb1), y = reorder(gene, deltaBeta_72hrs), fill = Beta)) +
  geom_tile() +
  xlab(label = "Sample") + 
  #remove x and y axis labels
  labs(x="",y="")+
  theme_grey(base_size=8)+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+ 
  ggtitle("CpG control vs. 72hrs IL6 vs. 96hrs wash") +
  scale_alpha( trans = "log" ) +
  scale_fill_gradient2(low = "#353436",
                       mid = "#f6f805",
                       high = "red",
                       guide = "colorbar") +
  coord_equal() + #square 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

p10

p2 + coord_flip()

#Global B values
sapply(raw_72hrs, class)
colMeans(raw_72hrs[ , c(2:14)])

# Change name of 1st column of df to "a"
names(df)[1] <- "a"

# Change name of 2nd column of df to "b"
names(df)[2] <- "b"

## Venn genes

set1[, c("Feature")]

set1_gene<- filter(set1, Feature %in% c("tss_body"))
set2_gene<- filter(set2, Feature %in% c("tss_body"))

set1_gene <- set1_gene$Gene %>% unique()
set2_gene <- set2_gene$Gene %>% unique()

Venn(list(set1_gene, set2_gene), SetNames = c("LPS Wt", "LPS SOCS3iEKO")) %>% plot

# Hierar Cluster
setwd("~/Analises R/r_array_huvec_120121/all_comparison_72hrs_96hrs_champ/EB2022/03_12_22")

df <- read.csv(file="./champ_72hrs_cluster.csv")
df.t <- t(df)
# Euclidean distance
dist.t <- dist(df.t[ -1, ] , diag=TRUE) #If want to select specif columns [ , c(4:8)]

# Hierarchical Clustering with hclust
hc.t <- hclust(dist.t)

# Plot the result
plot(hc.t)

#Heatmap for genes that match with RNA-seq and methylation

#Select the samples that will be in the heatmap

raw_champ_72hrs_heatm <- read.csv(file="./champ_72hrs_heatm.csv")

nostrin_champ_72hrs_heatm<- filter(raw_champ_72hrs_heatm, SampleID %in% c("cg06597454", "cg25694915", "cg16780603"))

write.csv(eb_champ_72hrs_heatm,file="./eb_champ_72hrs_heatm.csv")

nostrin_champ_72hrs_heatm <- read.csv(file="./nostrin_serpina3_plce1_champ_72hrs_heatm.csv")

# Select collumns that are necessary to the heatmap
nostrin_sample_heat <- nostrin_champ_72hrs_heatm [, c( 1, 15:28)]

# Pivot make the matrix to heatmap
nostrin_sample_heat_pivot <- nostrin_sample_heat %>% 
  mutate_if(is.numeric,as.character, is.factor, as.character)
sapply(nostrin_sample_heat_pivot, class)

nostrin_sample_heat_pivot <- pivot_longer(data = nostrin_sample_heat_pivot,
                                     cols = -c(1,2), 
                                     names_to = "Sample", 
                                     values_to = "Beta")

# Class of variables, Beta need to be numeric
sapply(nostrin_sample_heat_pivot, class)
nostrin_sample_heat_pivot$Beta <- as.numeric(nostrin_sample_heat_pivot$Beta)
nostrin_sample_heat_pivot$average_pbs <- as.numeric(nostrin_sample_heat_pivot$average_pbs)

#Order the values on the X
level_order_nostrin <- c('hPBS_exp2_1', 'hPBS_exp3_1', 'hPBS_exp2_2', 'hPBS_exp3_2', 'hPBS_exp1_1', 'hPBS_exp1_2', 'hIL6_exp2_1',
                    'hIL6_exp3_1', 'hIL6_exp2_2', 'hIL6_exp3_2', 'hIL6_exp1_1', 'hIL6_exp1_2', 'hIL6_exp1_3') #this vector might be useful for other plots/analyses

p11 <- ggplot(data = nostrin_sample_heat_pivot, mapping = aes(x = factor(Sample, level= level_order_nostrin), y = SampleID, fill = Beta)) +
  geom_tile() +
  xlab(label = "Sample") + 
  #remove x and y axis labels
  #labs(x="",y="")+
  theme_grey(base_size=8)+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+ 
  #ggtitle("CpG control vs. 72hrs IL6") +
  scale_alpha( trans = "log" ) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       guide = "colorbar") +
  coord_equal() + #square 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

p11

p11 + coord_flip()

+ #square 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

p10

p10 + coord_flip()

#Genes Hypermethylated and changes expression

nav2_tnfsf4_champ_72hrs_heatm <- read.csv(file="./nav2_tnfsf4_champ_72hrs_heatm.csv")

# Select collumns that are necessary to the heatmap
nav2_tnfsf4_sample_heat <- nav2_tnfsf4_champ_72hrs_heatm [, c( 1, 15:28)]

# Pivot make the matrix to heatmap
nav2_tnfsf4_sample_heat_pivot <- nav2_tnfsf4_sample_heat %>% 
  mutate_if(is.numeric,as.character, is.factor, as.character)
sapply(nav2_tnfsf4_sample_heat_pivot, class)

nav2_tnfsf4_sample_heat_pivot <- pivot_longer(data = nav2_tnfsf4_sample_heat_pivot,
                                          cols = -c(1,2), 
                                          names_to = "Sample", 
                                          values_to = "Beta")

# Class of variables, Beta need to be numeric
sapply(nav2_tnfsf4_sample_heat_pivot, class)
nav2_tnfsf4_sample_heat_pivot$Beta <- as.numeric(nav2_tnfsf4_sample_heat_pivot$Beta)
nav2_tnfsf4_sample_heat_pivot$average_pbs <- as.numeric(nav2_tnfsf4_sample_heat_pivot$average_pbs)

#Order the values on the X
level_order_nav2_tnfsf4 <- c('hPBS_exp2_1', 'hPBS_exp3_1', 'hPBS_exp2_2', 'hPBS_exp3_2', 'hPBS_exp1_1', 'hPBS_exp1_2', 'hIL6_exp2_1',
                         'hIL6_exp3_1', 'hIL6_exp2_2', 'hIL6_exp3_2', 'hIL6_exp1_1', 'hIL6_exp1_2', 'hIL6_exp1_3') #this vector might be useful for other plots/analyses

p12 <- ggplot(data = nav2_tnfsf4_sample_heat_pivot, mapping = aes(x = factor(Sample, level= level_order_nav2_tnfsf4), y = SampleID, fill = Beta)) +
  geom_tile() +
  xlab(label = "Sample") + 
  #remove x and y axis labels
  #labs(x="",y="")+
  theme_grey(base_size=8)+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+ 
  #ggtitle("CpG control vs. 72hrs IL6") +
  scale_alpha( trans = "log" ) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       guide = "colorbar") +
  coord_equal()
p12

+ #square 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 


p11 + coord_flip()

+ #square 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

p10

p10 + coord_flip()




#####
# DADA2 pipeline for 16S data from a greenhouse experiment
# Project FAPESP 2021/10626-0
# Author: Anderson Freitas
####

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.17")


#Libraries
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(microbiome)
library(microeco)
library(file2meco)
library(magrittr)
library(igraph)
library(meconetcomp)


ze = c("#b556bd",
       "#60c451",
       "#7166d9",
       "#9dbb36",
       "#d84488",
       "#43c07c",
       "#cf3c40",
       "#40c0bc",
       "#df6f31",
       "#3e87c4",
       "#daad3b",
       "#7777c1",
       "#549434",
       "#d687c0",
       "#3e7f45",
       "#9b4669",
       "#73bc8a",
       "#d97577",
       "#308568",
       "#a1512e",
       "#60b3e4",
       "#b17f26",
       "#606d2b",
       "#da9767",
       "#898e2b",
       "#886b30",
       "#b1b36a")


####DADA2

#############################################################
#P1

library(dada2); packageVersion("dada2")

path <- "./all/P1" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

fnFs <- sort(list.files(path, pattern=".effective.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), ".e"), `[`, 1)
plotQualityProfile(fnFs[2:5])

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names
out <- filterAndTrim(fnFs, filtFs,
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=F) # On Windows set multithread=FALSE
head(out)

save.image(file = "./all/P1/Dada.RData")

errF <- learnErrors(filtFs, multithread=F)
save.image(file = "./all/P1/Dada2.RData")
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, nbases = 1e+08, multithread=F, randomize=TRUE)

save.image(file = "./all/P1/Dada3.RData")

derepFs <- derepFastq(filtFs, verbose=F)

save.image(file = "./all/P1/Dada4.RData")
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepFs)

seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=F, verbose=TRUE)
save.image(file = "./all/P1/Dada5.RData")
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "nonchim")
head(track)

taxa1 <- assignTaxonomy(seqtab.nochim,
                       "./silva_nr99_v138.1_train_set.fa.gz",
                       multithread=F)
taxa.print <- taxa1 # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(seqtab.nochim, "./all/P1/p1_seqtab_final.rds")
saveRDS(taxa1, "./all/P1/_p1_taxonomy_genera.rds")

#############################################################
gc()

#P2

library(dada2); packageVersion("dada2")

path <- "./all/P2" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

fnFs <- sort(list.files(path, pattern=".effective.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), ".e"), `[`, 1)
plotQualityProfile(fnFs[2:5])

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names
out <- filterAndTrim(fnFs, filtFs,
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=F) # On Windows set multithread=FALSE
head(out)

save.image(file = "./all/P2/Dada.RData")

errF <- learnErrors(filtFs, multithread=F)
save.image(file = "./all/P2/Dada2.RData")
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, nbases = 1e+08, multithread=F, randomize=TRUE)

save.image(file = "./all/P2/Dada3.RData")
gc()
derepFs <- derepFastq(filtFs, verbose=F)

save.image(file = "./all/P2/Dada4.RData")
# Name the derep-class objects by the sample names
gc()
names(derepFs) <- sample.names
names(derepFs)

seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
gc()
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=F, verbose=TRUE)
save.image(file = "./all/P2/Dada5.RData")
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "nonchim")
head(track)
gc()
taxa2 <- assignTaxonomy(seqtab.nochim,
                        "./silva_nr99_v138.1_train_set.fa.gz",
                        multithread=F)

taxa.print <- taxa2 # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(seqtab.nochim, "./all/P2/p2_seqtab_final.rds")
saveRDS(taxa2, "./all/P2/_p2_taxonomy_genera.rds")

#############################################################
gc()

#P3

library(dada2); packageVersion("dada2")

path <- "./all/P3" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

fnFs <- sort(list.files(path, pattern=".effective.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), ".e"), `[`, 1)
plotQualityProfile(fnFs[2:5])

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names
out <- filterAndTrim(fnFs, filtFs,
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=F) # On Windows set multithread=FALSE
head(out)

save.image(file = "./all/P3/Dada.RData")

errF <- learnErrors(filtFs, multithread=F)
save.image(file = "./all/P3/Dada2.RData")
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, nbases = 1e+08, multithread=F, randomize=TRUE)
gc()
save.image(file = "./all/P3/Dada3.RData")

derepFs <- derepFastq(filtFs, verbose=F)

save.image(file = "./all/P3/Dada4.RData")
gc()
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepFs)

seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
gc()
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=F, verbose=TRUE)
save.image(file = "./all/P3/Dada5.RData")
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "nonchim")
head(track)
gc()
taxa3 <- assignTaxonomy(seqtab.nochim,
                        "./silva_nr99_v138.1_train_set.fa.gz",
                        multithread=F)
taxa.print <- taxa3 # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(seqtab.nochim, "./all/P3/p3_seqtab_final.rds")
saveRDS(taxa3, "./all/P3/_p3_taxonomy_genera.rds")




####Phyloseq

#############################################################
#Phyloseq
library(phyloseq)
library(tidyverse)
library(dplyr)
library(purrr)
library(ggplot2)

#input1
setwd("C:/Users/andersonfreitas/OneDrive/Anderson-BackUp/FAPEAM/Paper4_Greenhouse")
seqtab1 = readRDS("./Data/all/P1/p1_seqtab_final.rds")
taxa1 = readRDS("./Data/all/P1/_p1_taxonomy_genera.rds")
map <- "./map07_23.txt"
ps1 <- phyloseq(otu_table(seqtab1, taxa_are_rows=F), 
               tax_table(taxa1))
sample_metadata = import_qiime_sample_data(map)
input1 =merge_phyloseq(ps1, sample_metadata)
input1

#input2
seqtab2 = readRDS("./Data/all/P2/p2_seqtab_final.rds")
taxa2 = readRDS("./Data/all/P2/_p2_taxonomy_genera.rds")
ps2 <- phyloseq(otu_table(seqtab2, taxa_are_rows=F), 
                tax_table(taxa2))
input2 =merge_phyloseq(ps2, sample_metadata)
input2

#input3
seqtab3 = readRDS("./Data/all/P3/p3_seqtab_final.rds")
taxa3 = readRDS("./Data/all/P3/_p3_taxonomy_genera.rds")
ps3 <- phyloseq(otu_table(seqtab3, taxa_are_rows=F), 
                tax_table(taxa3))
input3 =merge_phyloseq(ps3, sample_metadata)
input3

merged_physeq <- merge_phyloseq(input1, input2, input3)
merged_physeq
input.final   <- subset_samples(merged_physeq, Tree != "None")
input.final
input.final2  <- subset_samples(input.final, 
                                Tree != "None")
input.final3  <- subset_samples(input.final2,
                                SampleID != "TP1" & SampleID != "TP2" & SampleID != "TP3" & SampleID != "TP4")
input.final3


####Beta div

#############################################################
#Beta Diversity

#Control x ADE
ctrl.ade  <- subset_samples(input.final, Substrate == "ADE" | Substrate == "Control")
input.clr = microbiome::transform(ctrl.ade, "clr")
df        = as(sample_data(input.clr), "data.frame")
ds        = phyloseq::distance(input.clr, method = "euclidean")
permanova = vegan::adonis2(ds ~ Tree*Substrate, data = df, permutations = 999)
permanova
#And the plot
ze=c("#79b643", "#bb59bd", "#52ac7a", "#cf426b", "#4dadcf","#cf5234",
     "#7475cb", "#cea13d", "#be6c91", "#748037", "#be774e")

input_ord = ordinate(input.clr, "PCoA" , "euclidean") 
p4 = plot_ordination(input.clr, input_ord, color = "Substrate", shape = "Tree")
p1.ctrl.ade = p4 + geom_point(aes(shape = Tree), size = 6, alpha = 0.8) +
  theme(legend.position = "right") +
  annotate("text", x = 50, y = -48, hjust = 0.2 , 
           label = bquote('Substrate:'~R^2~'= 0.04  |  p = 0.001'), size = 3)+
  annotate("text", x = 50, y = -52, hjust = 0.2 , 
           label = bquote('Tree:'~R^2~'= 0.11  |  p = 0.001'), size = 3)+
  scale_fill_manual(values = ze) +
  theme(axis.text.x = element_text(size = 20))    +
  theme_bw()
p1.ctrl.ade

dev.print(tiff, "./Figures/Beta_euc_PCoA_Ctrl_ADE_with_bulk.tiff", compression = "lzw", res=600, height=5, width=8, units="in")

#Control x ADE + CS
ctrl.ade  <- subset_samples(input.final, Substrate == "ADE + CS" | Substrate == "Control")
input.clr = microbiome::transform(ctrl.ade, "clr")
df        = as(sample_data(input.clr), "data.frame")
ds        = phyloseq::distance(input.clr, method = "euclidean")
permanova = vegan::adonis2(ds ~ Tree*Substrate, data = df, permutations = 999)
permanova
#And the plot
ze=c("#79b643", "#bb59bd", "#52ac7a", "#cf426b", "#4dadcf","#cf5234",
     "#7475cb", "#cea13d", "#be6c91", "#748037", "#be774e")

input_ord = ordinate(input.clr, "PCoA" , "euclidean") 
p4 = plot_ordination(input.clr, input_ord, color = "Substrate", shape = "Tree")
p1.ctrl.adecs = p4 + geom_point(aes(shape = Tree), size = 6, alpha = 0.8) +
  theme(legend.position = "right") +
  annotate("text", x = -50, y = -48, hjust = 0.2 , 
           label = bquote('Substrate:'~R^2~'= 0.06  |  p = 0.001'), size = 3)+
  annotate("text", x = -50, y = -52, hjust = 0.2 , 
           label = bquote('Tree:'~R^2~'= 0.10  |  p = 0.001'), size = 3)+
  scale_fill_manual(values = ze) +
  theme_bw()
p1.ctrl.adecs

dev.print(tiff, "./Figures/Beta_euc_PCoA_Ctrl_ADECS_with_bulk.tiff", compression = "lzw", res=600, height=5, width=8, units="in")


#Control x CS
ctrl.ade  <- subset_samples(input.final, Substrate == "CS" | Substrate == "Control")
input.clr = microbiome::transform(ctrl.ade, "clr")
df        = as(sample_data(input.clr), "data.frame")
ds        = phyloseq::distance(input.clr, method = "euclidean")
permanova = vegan::adonis2(ds ~ Tree*Substrate, data = df, permutations = 999)
permanova
#And the plot
ze=c("#79b643", "#bb59bd", "#52ac7a", "#cf426b", "#4dadcf","#cf5234",
     "#7475cb", "#cea13d", "#be6c91", "#748037", "#be774e")

input_ord = ordinate(input.clr, "PCoA" , "euclidean") 
p4 = plot_ordination(input.clr, input_ord, color = "Substrate", shape = "Tree")
p1.ctrl.cs = p4 + geom_point(aes(shape = Tree), size = 6, alpha = 0.8) +
  theme(legend.position = "right") +
  annotate("text", x = 20, y = -48, hjust = 0.2 , 
           label = bquote('Substrate:'~R^2~'= 0.05  |  p = 0.001'), size = 3)+
  annotate("text", x = 20, y = -52, hjust = 0.2 , 
           label = bquote('Tree:'~R^2~'= 0.9  |  p = 0.001'), size = 3)+
  scale_fill_manual(values = ze) +
  theme_bw()
p1.ctrl.cs

dev.print(tiff, "./Figures/Beta_euc_PCoA_Ctrl_cs_with_bulk.tiff", compression = "lzw", res=600, height=5, width=8, units="in")

#microbiome::summarize_phyloseq(merged_physeq) #131139 minimum

#BETA for ALL -Initial
input.clr = microbiome::transform(input.final, "clr")
input.clr.all = microbiome::transform(merged_physeq, "clr")
df        = as(sample_data(input.clr), "data.frame")
ds        = phyloseq::distance(input.clr, method = "euclidean")
permanova = vegan::adonis2(ds ~ Tree*Substrate, data = df, permutations = 999)
permanova
#And the plot
ze=c("#79b643", "#bb59bd", "#52ac7a", "#cf426b", "#4dadcf","#cf5234",
     "#7475cb", "#cea13d", "#be6c91", "#748037", "#be774e")
input_ord = ordinate(input.clr, "PCoA" , "euclidean") 
p4 = plot_ordination(input.clr, input_ord, color = "Substrate", shape = "Tree")
p1.all.clr = p4 + geom_point(aes(shape = Tree), size = 6, alpha = 0.8) +
  theme(legend.position = "right") +
  annotate("text", x = -60, y = -48, hjust = 0.2 , 
           label = bquote('Substrate:'~R^2~'= 0.09  |  p = 0.001'), size = 3)+
  annotate("text", x = -60, y = -52, hjust = 0.2 , 
           label = bquote('Tree:'~R^2~'= 0.05  |  p = 0.001'), size = 3)+
  scale_fill_manual(values = ze) +
  theme_bw()
p1.all.clr

dev.print(tiff, "./Figures/Beta_euc_allsamples_with_bulk.tiff", compression = "lzw", res=600, height=7, width=12, units="in")

##### Subsetting
Schi  = subset_samples(input.final, Tree == "S. amazonicum")
Cec   = subset_samples(input.final, Tree == "Cecropia sp.")
Aca   = subset_samples(input.final, Tree == "A. mangium")
Ipe   = subset_samples(input.final, Tree == "H. avellaneadae")
Bri   = subset_samples(input.final, Tree == "U. brizantha")
Bulk  = subset_samples(input.final, Tree == "Bulk")
####

#Beta for Schizolobium
input.Schi.clr = microbiome::transform(Schi, "clr")
df        = as(sample_data(input.Schi.clr), "data.frame")
ds        = phyloseq::distance(input.Schi.clr, method = "euclidean")
permanova = vegan::adonis2(ds ~ Substrate, data = df, permutations = 999)
permanova
#And the plot
ze=c("#79b643", "#bb59bd", "#52ac7a", "#cf426b", "#4dadcf","#cf5234",
     "#7475cb", "#cea13d", "#be6c91", "#748037", "#be774e")
input_ord = ordinate(input.Schi.clr, "PCoA" , "euclidean") 
p4 = plot_ordination(input.Schi.clr, input_ord, color = "Substrate")
p1.schi = p4 + geom_point(size = 6, alpha = 0.8) +
  theme(legend.position = "right") +
  annotate("text", x = -50, y = -30, hjust = 0.2 , 
           label = bquote(''~R^2~'= 0.23  |  p = 0.001'), size = 3)+
  scale_fill_manual(values = ze) +
  labs(title = "S. amazonicum") +
  theme_bw()
p1.schi

dev.print(tiff, "./Figures/Beta_euc_parica.tiff", compression = "lzw", res=600, height=5, width=8, units="in")


#Beta diversity for Schizolobium with bray distance
df        = as(sample_data(Schi), "data.frame")
ds        = phyloseq::distance(Schi, method = "bray")
permanova = vegan::adonis2(ds ~ Substrate, data = df, permutations = 999)
permanova
#And the plot
ze=c("#79b643", "#bb59bd", "#52ac7a", "#cf426b", "#4dadcf","#cf5234",
     "#7475cb", "#cea13d", "#be6c91", "#748037", "#be774e")
input_ord = ordinate(Schi, "NMDS" , "bray") 
p4 = plot_ordination(Schi, input_ord, color = "Substrate")
p1.schi.bray = p4 + geom_point(size = 6, alpha = 0.8) +
  theme(legend.position = "right") +
  annotate("text", x = 0, y = 0, hjust = 0.2 , 
           label = bquote(''~R^2~'= 0.33  |  p = 0.001'), size = 3)+
  scale_fill_manual(values = ze) +
  labs(title = "S. amazonicum") +
  theme_bw()
p1.schi.bray
## was almost the same


#Beta diversity for Cecropia pachystachia
input.Cec.clr = microbiome::transform(Cec, "clr")
df        = as(sample_data(input.Cec.clr), "data.frame")
ds        = phyloseq::distance(input.Cec.clr, method = "euclidean")
permanova = vegan::adonis2(ds ~ Substrate, data = df, permutations = 999)
permanova
#And the plot
ze=c("#79b643", "#bb59bd", "#52ac7a", "#cf426b", "#4dadcf","#cf5234",
     "#7475cb", "#cea13d", "#be6c91", "#748037", "#be774e")
input_ord = ordinate(input.Cec.clr, "PCoA" , "euclidean") 
p4 = plot_ordination(input.Cec.clr, input_ord, color = "Substrate")
p1.cec = p4 + geom_point(size = 6, alpha = 0.8) +
  theme(legend.position = "right") +
  annotate("text", x = -55, y = 40, hjust = 0.2 , 
           label = bquote(''~R^2~'= 0.23  |  p = 0.006'), size = 3)+
  scale_fill_manual(values = ze) +
  labs(title = "C. pachystachya") +
  theme_bw()
p1.cec

dev.print(tiff, "./Figures/Beta_euc_embauba.tiff", compression = "lzw", res=600, height=5, width=8, units="in")


#Beta diversity for Urochloa brizantha
input.Bri.clr = microbiome::transform(Bri, "clr")
df        = as(sample_data(input.Bri.clr), "data.frame")
ds        = phyloseq::distance(input.Bri.clr, method = "euclidean")
permanova = vegan::adonis2(ds ~ Substrate, data = df, permutations = 999)
permanova
#And the plot
ze=c("#79b643", "#bb59bd", "#52ac7a", "#cf426b", "#4dadcf","#cf5234",
     "#7475cb", "#cea13d", "#be6c91", "#748037", "#be774e")
input_ord = ordinate(input.Bri.clr, "PCoA" , "euclidean") 
p4 = plot_ordination(input.Bri.clr, input_ord, color = "Substrate")
p1.bri = p4 + geom_point(size = 6, alpha = 0.8) +
  theme(legend.position = "right") +
  annotate("text", x = 50, y = -40, hjust = 0.2 , 
           label = bquote(''~R^2~'= 0.28  |  p = 0.001'), size = 3)+
  scale_fill_manual(values = ze) +
  labs(title = "U. brizantha") +
  theme_bw()
p1.bri

dev.print(tiff, "./Figures/Beta_euc_braquiaria_with_initial.tiff", compression = "lzw", res=600, height=5, width=8, units="in")


#Beta diversity for Bulk soil
input.Bulk.clr = microbiome::transform(Bulk, "clr")
df        = as(sample_data(input.Bulk.clr), "data.frame")
ds        = phyloseq::distance(input.Bulk.clr, method = "euclidean")
permanova = vegan::adonis2(ds ~ Substrate, data = df, permutations = 999)
permanova
#And the plot
ze=c("#79b643", "#bb59bd", "#52ac7a", "#cf426b", "#4dadcf","#cf5234",
     "#7475cb", "#cea13d", "#be6c91", "#748037", "#be774e")
input_ord = ordinate(input.Bulk.clr, "PCoA" , "euclidean") 
p4 = plot_ordination(input.Bulk.clr, input_ord, color = "Substrate")
p1.bulk = p4 + geom_point(size = 6, alpha = 0.8) +
  theme(legend.position = "right") +
  annotate("text", x = 80, y = -40, hjust = 0.2 , 
           label = bquote(''~R^2~'= 0.23  |  p = 0.001'), size = 3)+
  scale_fill_manual(values = ze) +
  labs(title = "Bulk Soil") +
  theme_bw()
p1.bulk

dev.print(tiff, "./Figures/Beta_euc_bulk_with_initial.tiff", compression = "lzw", res=600, height=5, width=8, units="in")


#Beta diversity for Handroanthus avellanedae
input.ipe.clr = microbiome::transform(Ipe, "clr")
df        = as(sample_data(input.ipe.clr), "data.frame")
ds        = phyloseq::distance(input.ipe.clr, method = "euclidean")
permanova = vegan::adonis2(ds ~ Substrate, data = df, permutations = 999)
permanova
#And the plot
ze=c("#79b643", "#bb59bd", "#52ac7a", "#cf426b", "#4dadcf","#cf5234",
     "#7475cb", "#cea13d", "#be6c91", "#748037", "#be774e")
input_ord = ordinate(input.ipe.clr, "PCoA" , "euclidean") 
p4 = plot_ordination(input.ipe.clr, input_ord, color = "Substrate")
p1.ipe = p4 + geom_point(size = 6, alpha = 0.8) +
  theme(legend.position = "right") +
  annotate("text", x = 35, y = -40, hjust = 0.2 , 
           label = bquote(''~R^2~'= 0.25  |  p = 0.001'), size = 3)+
  scale_fill_manual(values = ze) +
  labs(title = "H. avellanedae") +
  theme_bw()
p1.ipe
dev.print(tiff, "./Figures/Beta_euc_ipe.tiff", compression = "lzw", res=600, height=5, width=8, units="in")


#Beta diversity for Acacia mangium
input.aca.clr = microbiome::transform(Aca, "clr")
df        = as(sample_data(input.aca.clr), "data.frame")
ds        = phyloseq::distance(input.aca.clr, method = "euclidean")
permanova = vegan::adonis2(ds ~ Substrate, data = df, permutations = 999)
permanova
#And the plot
ze=c("#79b643", "#bb59bd", "#52ac7a", "#cf426b", "#4dadcf","#cf5234",
     "#7475cb", "#cea13d", "#be6c91", "#748037", "#be774e")
input_ord = ordinate(input.aca.clr, "PCoA" , "euclidean") 
p4 = plot_ordination(input.aca.clr, input_ord, color = "Substrate")
p1.aca = p4 + geom_point(size = 6, alpha = 0.8) +
  theme(legend.position = "right") +
  annotate("text", x = 35, y = -40, hjust = 0.2 , 
           label = bquote(''~R^2~'= 0.23  |  p = 0.001'), size = 3)+
  scale_fill_manual(values = ze) +
  labs(title = "A. mangium") +
  theme_bw()
p1.aca
dev.print(tiff, "./Figures/Beta_euc_aca.tiff", compression = "lzw", res=600, height=5, width=8, units="in")

#######

####RF

#############################################################
#Randon Forest
library(randomForest)

#RF for all
all.genus = microbiome::aggregate_rare(merged_physeq, level = "Genus",
                                 detection = 1/100, prevalence = 1/100)
all.genus.clr = microbiome::transform(all.genus, "clr")
predictors <- t(otu_table(all.genus.clr))
dim(predictors)
response <- as.factor(sample_data(all.genus.clr)$Treatment)
rf.data <- data.frame(response, predictors)
set.seed(2)
erie.classify <- randomForest(response~., data = rf.data, ntree = 1000)
print(erie.classify)
#not good at all


#RF for all -Initial
final.genus = microbiome::aggregate_rare(input.final, level = "Genus",
                                       detection = 1/100, prevalence = 1/100)
final.genus.clr = microbiome::transform(final.genus, "clr")
predictors <- t(otu_table(final.genus.clr))
response <- as.factor(sample_data(final.genus.clr)$Substrate)
rf.data <- data.frame(response, predictors)
set.seed(234)
erie.classify <- randomForest(response~., data = rf.data, ntree = 10000)
print(erie.classify)
#OOB = 16.79%. Biggest error rates are in Control.


#RF for Bulk Soil
Bulk.ag = microbiome::aggregate_rare(Bulk, level = "Genus",
                                       detection = 1/100, prevalence = 1/100)
Bulk.ag.clr = microbiome::transform(Bulk.ag, "clr")
predictors <- t(otu_table(Bulk.ag.clr))
dim(predictors)
response <- as.factor(sample_data(Bulk.ag.clr)$Substrate)
rf.data <- data.frame(response, predictors)
set.seed(2)
erie.classify <- randomForest(response~., data = rf.data, ntree = 10000)
print(erie.classify)
#OOB = 43.8%. Bad.


#RF for Schizolobium Amazonicum
Schi.ag = microbiome::aggregate_rare(Schi, level = "Genus",
                                 detection = 1/100, prevalence = 1/100)
Schi.ag.clr = microbiome::transform(Schi.ag, "clr")
predictors <- t(otu_table(Schi.ag.clr))
dim(predictors)
response <- as.factor(sample_data(Schi.ag.clr)$Substrate)
rf.data <- data.frame(response, predictors)
set.seed(2)
erie.classify <- randomForest(response~., data = rf.data, ntree = 10000)
print(erie.classify)
#OOB = 30%. All CS are correctly classified. But bad at all.


#RF for Cecropia pachystachya
Cec.ag = microbiome::aggregate_rare(Cec, level = "Genus",
                                     detection = 1/100, prevalence = 1/100)
Cec.ag.clr = microbiome::transform(Cec.ag, "clr")
predictors <- t(otu_table(Cec.ag.clr))
dim(predictors)
response <- as.factor(sample_data(Cec.ag.clr)$Substrate)
rf.data <- data.frame(response, predictors)
set.seed(2)
erie.classify <- randomForest(response~., data = rf.data, ntree = 10000)
print(erie.classify)
#OOB = 31.58%. All ADE+CS are correctly classified. But bad at all.


#RF for Urochloa brizantha
Bri.ag = microbiome::aggregate_rare(Bri, level = "Genus",
                                    detection = 1/100, prevalence = 1/100)
Bri.ag.clr = microbiome::transform(Bri.ag, "clr")
predictors <- t(otu_table(Bri.ag.clr))
dim(predictors)
response <- as.factor(sample_data(Bri.ag.clr)$Substrate)
rf.data <- data.frame(response, predictors)
set.seed(2)
erie.classify <- randomForest(response~., data = rf.data, ntree = 10000)
print(erie.classify)
#OOB = 31.43%. All ADE+CS are correctly classified. But bad at all.


#RF for Handroanthus avellanedae
Ipe.ag = microbiome::aggregate_rare(Ipe, level = "Genus",
                                    detection = 1/100, prevalence = 1/100)
Ipe.ag.clr = microbiome::transform(Ipe.ag, "clr")
predictors <- t(otu_table(Ipe.ag.clr))
dim(predictors)
response <- as.factor(sample_data(Ipe.ag.clr)$Substrate)
rf.data <- data.frame(response, predictors)
set.seed(2)
erie.classify <- randomForest(response~., data = rf.data, ntree = 10000)
print(erie.classify)
#OOB = 20%. It classified wrongly 1/5 in each treatment


#RF for Acacia mangium
Aca.ag = microbiome::aggregate_rare(Aca, level = "Genus",
                                    detection = 1/100, prevalence = 1/100)
Aca.ag.clr = microbiome::transform(Aca.ag, "clr")
predictors <- t(otu_table(Aca.ag.clr))
dim(predictors)
response <- as.factor(sample_data(Aca.ag.clr)$Substrate)
rf.data <- data.frame(response, predictors)
set.seed(2)
erie.classify <- randomForest(response~., data = rf.data, ntree = 10000)
print(erie.classify)
#OOB = 30%. Bad.


####FAPROTAX

#############################################################

###FAPROTAX
library(file2meco)
final.r6 <- phyloseq2meco(input.final)
final.r6
library(microeco)
t1 <- trans_func$new(final.r6)
t1$cal_spe_func(prok_database = "FAPROTAX")
t1$cal_spe_func_perc(abundance_weighted = TRUE)
genes = cbind(as.data.frame(t1$res_spe_func_perc),sample_data(input.final))
colnames(genes)
selgenes <- genes[c(1,5:8,10:12,17,18,20:22,24,34,46,50:53,57,58,59)]


#All of them heatmap
allgenes <- selgenes[-c(21,22)] %>% 
  group_by(Treatment) %>%
  summarise_all("mean")
allgenes$Treatment <- factor(allgenes$Treatment,
                             levels = c("Bulk Initial" , "U. brizantha Initial"  ,
                               "Bulk Control","U. brizantha Control","S. amazonicum Control","Cecropia sp. Control","H. avellaneadae Control", "A. mangium Control",
                                         "Bulk CS" ,"U. brizantha CS", "S. amazonicum CS","Cecropia sp. CS","H. avellaneadae CS", "A. mangium CS",
                                         "Bulk ADE", "U. brizantha ADE","S. amazonicum ADE","Cecropia sp. ADE","H. avellaneadae ADE", "A. mangium ADE",
                                          "Bulk ADE + CS","U. brizantha ADE + CS", "S. amazonicum ADE + CS","Cecropia sp. ADE + CS","H. avellaneadae ADE + CS", "A. mangium ADE + CS"
                                                            ))
data.all <- allgenes %>%
  tibble::column_to_rownames("Treatment") 
data.all = as.matrix(t(data.all))
heatmap(data.all,Colv = NA, Rowv = NA)
data_melt <- reshape2::melt(allgenes, id = c("Treatment"))
ggp <- ggplot(data_melt, aes(Treatment, variable)) +    
  geom_tile(aes(fill = value)) +
  labs(title = "ALL") + 
  theme_classic()
ggp + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradient(low = "lightblue", high = "red")     
#Saved the PDF

#S. amazonicum heatmap
Schi.genes <- filter(selgenes, Tree == "S. amazonicum")
Schi.gen.p <- Schi.genes[-22] %>% 
  group_by(Substrate) %>%
  summarise_all("mean")
Schi.gen.p$Substrate <- factor(Schi.gen.p$Substrate, levels = c("Control", "CS", "ADE", "ADE + CS"))
data.schi <- Schi.gen.p %>%
  tibble::column_to_rownames("Substrate") 
data.schi = as.matrix(t(data.schi))
heatmap(data.schi,Colv = NA, Rowv = NA)
data_melt <- reshape2::melt(Schi.gen.p)
ggp <- ggplot(data_melt, aes(Substrate, variable)) +  
  geom_tile(aes(fill = value)) +
  labs(title = "S. amazonicum") +
  theme_classic()
ggp + 
  scale_fill_gradient(low = "lightblue", high = "red")     
#Saved the PDF


#C. pachystachya heatmap
Cec.genes <- filter(selgenes, Tree == "Cecropia sp.")
Cec.gen.p <- Cec.genes[-22] %>% 
  group_by(Substrate) %>%
  summarise_all("mean")
Cec.gen.p$Substrate <- factor(Cec.gen.p$Substrate, levels = c("Control", "CS", "ADE", "ADE + CS"))
data.Cec <- Cec.gen.p %>%
  tibble::column_to_rownames("Substrate") 
data.Cec = as.matrix(t(data.Cec))
heatmap(data.Cec,Colv = NA, Rowv = NA)
data_melt <- reshape2::melt(Cec.gen.p)
ggp <- ggplot(data_melt, aes(Substrate, variable)) +    
  geom_tile(aes(fill = value)) +
  labs(title = "C. pachystachya") +
  theme_classic()
ggp + 
  scale_fill_gradient(low = "lightblue", high = "red")     
#Saved the PDF


#U. brizantha heatmap
Bri.genes <- filter(selgenes, Tree == "U. brizantha")
Bri.gen.p <- Bri.genes[-c(22:23)] %>% 
  group_by(Substrate) %>%
  summarise_all("mean")
Bri.gen.p$Substrate <- factor(Bri.gen.p$Substrate, levels = c("Initial","Control", "CS", "ADE", "ADE + CS"))
data.Bri <- Bri.gen.p %>%
  tibble::column_to_rownames("Substrate") 
data.Bri = as.matrix(t(data.Bri))
heatmap(data.Bri,Colv = NA, Rowv = NA)
data_melt <- reshape2::melt(Bri.gen.p)
ggp <- ggplot(data_melt, aes(Substrate, variable)) +    
  geom_tile(aes(fill = value)) +
  labs(title = "U. brizantha") +
  theme_classic()
ggp + 
  scale_fill_gradient(low = "lightblue", high = "red")     
#Saved the PDF


#H. avellanedae heatmap
Ipe.genes <- filter(selgenes, Tree == "H. avellaneadae")
Ipe.gen.p <- Ipe.genes[-c(22:23)] %>% 
  group_by(Substrate) %>%
  summarise_all("mean")
Ipe.gen.p$Substrate <- factor(Ipe.gen.p$Substrate, levels = c("Control", "CS", "ADE", "ADE + CS"))
data.Ipe <- Ipe.gen.p %>%
  tibble::column_to_rownames("Substrate") 
data.Ipe = as.matrix(t(data.Ipe))
heatmap(data.Ipe,Colv = NA, Rowv = NA)
data_melt <- reshape2::melt(Ipe.gen.p)
ggp <- ggplot(data_melt, aes(Substrate, variable)) +    
  geom_tile(aes(fill = value)) +
  labs(title = "H. avellaneadae") +
  theme_classic()
ggp + 
  scale_fill_gradient(low = "lightblue", high = "red")     
#Saved the PDF


#A. mangium heatmap
Aca.genes <- filter(selgenes, Tree == "A. mangium")
Aca.gen.p <- Aca.genes[-c(22:23)] %>% 
  group_by(Substrate) %>%
  summarise_all("mean")
Aca.gen.p$Substrate <- factor(Aca.gen.p$Substrate, levels = c("Control", "CS", "ADE", "ADE + CS"))
data.Aca <- Aca.gen.p %>%
  tibble::column_to_rownames("Substrate") 
data.Aca = as.matrix(t(data.Aca))
heatmap(data.Aca,Colv = NA, Rowv = NA)
data_melt <- reshape2::melt(Aca.gen.p)
ggp <- ggplot(data_melt, aes(Substrate, variable)) +    
  geom_tile(aes(fill = value)) +
  labs(title = "A. mangium") +
  theme_classic()
ggp + 
  scale_fill_gradient(low = "lightblue", high = "red")     
#Saved the PDF



####Diferential abundance

##########################
library(ALDEx2)

# Schizolobium Control x CS
Schi.ag
schi.ag.cs = subset_samples(Schi.ag, Substrate == "Control" | Substrate == "CS")
mi         = as.data.frame((otu_table(schi.ag.cs)))
var        = sample_data(schi.ag.cs)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
               include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
res.schi.cs= x[(x$we.ep<="0.01"),]
#View(res.schi.cs)
#colnames(res.schi.cs)
filt       = cbind(as.data.frame(row.names(res.schi.cs)),res.schi.cs[c(2,3,16,17,20)])
colnames(filt)
colnames(filt)[which(names(filt) == "row.names(res.schi.cs)")] <- "Genera"
write.csv(x = filt, file = "./Figures/Schi_CtrlxCS.csv")


# Schizolobium Control x CS + ADE
schi.ag.cs = subset_samples(Schi.ag, Substrate == "Control" | Substrate == "ADE + CS")
mi         = as.data.frame((otu_table(schi.ag.cs)))
var        = sample_data(schi.ag.cs)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
res.schi.cs= x[(x$we.ep<="0.01"),]
#View(res.schi.cs)
#colnames(res.schi.cs)
filt2       = cbind(as.data.frame(row.names(res.schi.cs)),res.schi.cs[c(2,3,16,17,20)])
colnames(filt2)
colnames(filt2)[which(names(filt2) == "row.names(res.schi.cs)")] <- "Genera"
write.csv(x = filt2, file = "./Figures/Schi_CtrlxCS+ADE.csv")


# Schizolobium Control x ADE
schi.ag.cs = subset_samples(Schi.ag, Substrate == "Control" | Substrate == "ADE")
mi         = as.data.frame((otu_table(schi.ag.cs)))
var        = sample_data(schi.ag.cs)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
res.schi.cs= x[(x$we.ep<="0.01"),]
#View(res.schi.cs)
#colnames(res.schi.cs)
filt3       = cbind(as.data.frame(row.names(res.schi.cs)),res.schi.cs[c(2,3,16,17,20)])
colnames(filt3)
colnames(filt3)[which(names(filt2) == "row.names(res.schi.cs)")] <- "Genera"
write.csv(x = filt3, file = "./Figures/Schi_CtrlxADE.csv")

# ------- ##### ---------

# Cecropia Control x CS
Cec.ag
cec.ag.cs = subset_samples(Cec.ag, Substrate == "Control" | Substrate == "CS")
mi         = as.data.frame((otu_table(cec.ag.cs)))
var        = sample_data(cec.ag.cs)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
res.cec.cs= x[(x$we.ep<="0.01"),]
#View(res.cec.cs)
#colnames(res.cec.cs)
filt4       = cbind(as.data.frame(row.names(res.cec.cs)),res.cec.cs[c(2,3,15,16,19)])
colnames(filt4)
colnames(filt4)[which(names(filt4) == "row.names(res.cec.cs)")] <- "Genera"
write.csv(x = filt4, file = "./Figures/Cec_CtrlxCS.csv")


# Cecropia Control x CS + ADE
cec.ag.cs = subset_samples(Cec.ag, Substrate == "Control" | Substrate == "ADE + CS")
mi         = as.data.frame((otu_table(cec.ag.cs)))
var        = sample_data(cec.ag.cs)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
res.cec.cs= x[(x$we.ep<="0.01"),]
#View(res.cec.cs)
#colnames(res.cec.cs)
filt5       = cbind(as.data.frame(row.names(res.cec.cs)),res.cec.cs[c(2,3,15,16,19)])
colnames(filt5)
colnames(filt5)[which(names(filt5) == "row.names(res.cec.cs)")] <- "Genera"
write.csv(x = filt5, file = "./Figures/Cec_CtrlxCS+ADE.csv")


# Cecropia Control x ADE
cec.ag.cs = subset_samples(Cec.ag, Substrate == "Control" | Substrate == "ADE")
mi         = as.data.frame((otu_table(cec.ag.cs)))
var        = sample_data(cec.ag.cs)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
res.cec.cs= x[(x$we.ep<="0.01"),]
#View(res.cec.cs)
#colnames(res.cec.cs)
filt6       = cbind(as.data.frame(row.names(res.cec.cs)),res.cec.cs[c(2,3,15,16,19)])
colnames(filt6)
colnames(filt6)[which(names(filt6) == "row.names(res.cec.cs)")] <- "Genera"
write.csv(x = filt6, file = "./Figures/Cec_CtrlxADE.csv")


# ------- ##### ---------

# Acacia Control x CS
Aca.ag
aca.ag.cs = subset_samples(Aca.ag, Substrate == "Control" | Substrate == "CS")
mi         = as.data.frame((otu_table(aca.ag.cs)))
var        = sample_data(aca.ag.cs)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
res.aca.cs= x[(x$we.ep<="0.01"),]
#View(res.aca.cs)
#colnames(res.aca.cs)
filt7       = cbind(as.data.frame(row.names(res.aca.cs)),res.aca.cs[c(2,3,16,17,20)])
colnames(filt7)
colnames(filt7)[which(names(filt7) == "row.names(res.aca.cs)")] <- "Genera"
write.csv(x = filt7, file = "./Figures/aca_CtrlxCS.csv")


# Acacia Control x ADE + CS
aca.ag.cs  = subset_samples(Aca.ag, Substrate == "Control" | Substrate == "ADE + CS")
mi         = as.data.frame((otu_table(aca.ag.cs)))
var        = sample_data(aca.ag.cs)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
res.aca.cs= x[(x$we.ep<="0.01"),]
#View(res.aca.cs)
#colnames(res.aca.cs)
filt8       = cbind(as.data.frame(row.names(res.aca.cs)),res.aca.cs[c(2,3,16,17,20)])
colnames(filt8)
colnames(filt8)[which(names(filt8) == "row.names(res.aca.cs)")] <- "Genera"
write.csv(x = filt8, file = "./Figures/aca_CtrlxCS+ADE.csv")


# Acacia Control x ADE
aca.ag.cs  = subset_samples(Aca.ag, Substrate == "Control" | Substrate == "ADE")
mi         = as.data.frame((otu_table(aca.ag.cs)))
var        = sample_data(aca.ag.cs)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
res.aca.cs= x[(x$we.ep<="0.01"),]
#View(res.aca.cs)
#colnames(res.aca.cs)
filt9       = cbind(as.data.frame(row.names(res.aca.cs)),res.aca.cs[c(2,3,16,17,20)])
colnames(filt9)
colnames(filt9)[which(names(filt9) == "row.names(res.aca.cs)")] <- "Genera"
write.csv(x = filt9, file = "./Figures/aca_CtrlxCS+ADE.csv")


# ------- ##### ---------

# Handroanthus Control x CS
Ipe.ag
ipe.ag.cs  = subset_samples(Ipe.ag, Substrate == "Control" | Substrate == "CS")
mi         = as.data.frame((otu_table(ipe.ag.cs)))
var        = sample_data(ipe.ag.cs)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
res.ipe.cs= x[(x$we.ep<="0.01"),]
#View(res.ipe.cs)
#colnames(res.ipe.cs)
filt10       = cbind(as.data.frame(row.names(res.ipe.cs)),res.ipe.cs[c(2,3,16,17,20)])
colnames(filt10)
colnames(filt10)[which(names(filt10) == "row.names(res.ipe.cs)")] <- "Genera"
write.csv(x = filt10, file = "./Figures/ipe_CtrlxCS.csv")


# Handroanthus Control x CS + ADE
ipe.ag.cs  = subset_samples(Ipe.ag, Substrate == "Control" | Substrate == "ADE + CS")
mi         = as.data.frame((otu_table(ipe.ag.cs)))
var        = sample_data(ipe.ag.cs)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
res.ipe.cs= x[(x$we.ep<="0.01"),]
#View(res.ipe.cs)
#colnames(res.ipe.cs)
filt11       = cbind(as.data.frame(row.names(res.ipe.cs)),res.ipe.cs[c(2,3,16,17,20)])
colnames(filt11)
colnames(filt11)[which(names(filt11) == "row.names(res.ipe.cs)")] <- "Genera"
write.csv(x = filt11, file = "./Figures/ipe_CtrlxCS+ADE.csv")


# Handroanthus Control x ADE
ipe.ag.cs  = subset_samples(Ipe.ag, Substrate == "Control" | Substrate == "ADE")
mi         = as.data.frame((otu_table(ipe.ag.cs)))
var        = sample_data(ipe.ag.cs)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
res.ipe.cs= x[(x$we.ep<="0.01"),]
#View(res.ipe.cs)
#colnames(res.ipe.cs)
filt12       = cbind(as.data.frame(row.names(res.ipe.cs)),res.ipe.cs[c(2,3,16,17,20)])
colnames(filt12)
colnames(filt12)[which(names(filt12) == "row.names(res.ipe.cs)")] <- "Genera"
write.csv(x = filt12, file = "./Figures/ipe_CtrlxADE.csv")


# ------- ##### ---------

# Urochloa Control x CS
Bri.ag
bri.ag.cs  = subset_samples(Bri.ag, Substrate == "Control" | Substrate == "CS")
mi         = as.data.frame((otu_table(bri.ag.cs)))
var        = sample_data(bri.ag.cs)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
res.bri.cs= x[(x$we.ep<="0.01"),]
#View(res.bri.cs)
#colnames(res.bri.cs)
filt13       = cbind(as.data.frame(row.names(res.bri.cs)),res.bri.cs[c(2,3,14,15,18)])
colnames(filt13)
colnames(filt13)[which(names(filt13) == "row.names(res.bri.cs)")] <- "Genera"
write.csv(x = filt13, file = "./Figures/bri_CtrlxCS.csv")


# Urochloa Control x ADE + CS
bri.ag.cs  = subset_samples(Bri.ag, Substrate == "Control" | Substrate == "ADE + CS")
mi         = as.data.frame((otu_table(bri.ag.cs)))
var        = sample_data(bri.ag.cs)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
res.bri.cs= x[(x$we.ep<="0.01"),]
#View(res.bri.cs)
#colnames(res.bri.cs)
filt14       = cbind(as.data.frame(row.names(res.bri.cs)),res.bri.cs[c(2,3,14,15,18)])
colnames(filt14)
colnames(filt14)[which(names(filt14) == "row.names(res.bri.cs)")] <- "Genera"
write.csv(x = filt14, file = "./Figures/bri_CtrlxCS+ADE.csv")


# Urochloa Control x ADE
bri.ag.cs  = subset_samples(Bri.ag, Substrate == "Control" | Substrate == "ADE")
mi         = as.data.frame((otu_table(bri.ag.cs)))
var        = sample_data(bri.ag.cs)
treat      = var$Substrate
x          <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
                    include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
#View(x)
res.bri.cs= x[(x$we.ep<="0.01"),]
#View(res.bri.cs)
#colnames(res.bri.cs)
filt15       = cbind(as.data.frame(row.names(res.bri.cs)),res.bri.cs[c(2,3,14,15,18)])
colnames(filt15)
colnames(filt15)[which(names(filt15) == "row.names(res.bri.cs)")] <- "Genera"
write.csv(x = filt15, file = "./Figures/bri_CtrlxADE.csv")

################################################################################


####Networks

###############################################################################




###
#Network for Cecropia Control (5)
Cec.ctrl <- subset_samples(input.final, Tree == "Cecropia sp." & Substrate == "Control")
cec.ctrl <- phyloseq2meco(Cec.ctrl)
cec.ctrl$tax_table %<>% tidy_taxonomy
cec.ctrl$cal_abund()
t1 <- trans_network$new(dataset = cec.ctrl, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
cec.ctrl.csv = t1$res_network_attr
write.csv(cec.ctrl.csv, "./Figures/Cec_Ctrl_atributes.csv")
t1$save_network(filepath = "./Figures/Cec_Control_network.gexf")


#Network for Cecropia CS (6)
Cec.cs <- subset_samples(input.final, Tree == "Cecropia sp." & Substrate == "CS")
cec.cs <- phyloseq2meco(Cec.cs)
cec.cs$tax_table %<>% tidy_taxonomy
cec.cs$cal_abund()
t1 <- trans_network$new(dataset = cec.cs, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
cec.cs.csv = t1$res_network_attr
write.csv(cec.cs.csv, "./Figures/Cec_CS_atributes.csv")
t1$save_network(filepath = "./Figures/Cec_CS_network.gexf")


#Network for Cecropia ADE (7)
Cec.ade <- subset_samples(input.final, Tree == "Cecropia sp." & Substrate == "ADE")
cec.ade <- phyloseq2meco(Cec.ade)
cec.ade$tax_table %<>% tidy_taxonomy
cec.ade$cal_abund()
t1 <- trans_network$new(dataset = cec.ade, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
cec.ade.csv = t1$res_network_attr
write.csv(cec.ade.csv, "./Figures/Cec_ADE_atributes.csv")
t1$save_network(filepath = "./Figures/Cec_ADE_network.gexf")


#Network for Cecropia CS+ADE (8)
Cec.csade <- subset_samples(input.final, Tree == "Cecropia sp." & Substrate == "CS + ADE")
cec.csade <- phyloseq2meco(Cec.csade)
cec.csade$tax_table %<>% tidy_taxonomy
cec.csade$cal_abund()
t1 <- trans_network$new(dataset = cec.csade, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
cec.csade.csv = t1$res_network_attr
write.csv(cec.csade.csv, "./Figures/Cec_CSADE_atributes.csv")
t1$save_network(filepath = "./Figures/Cec_CSADE_network.gexf")


###
#Network for Handroanthus Control (9)
Ipe.ctrl <- subset_samples(input.final, Tree == "H. avellaneadae" & Substrate == "Control")
ipe.ctrl <- phyloseq2meco(Ipe.ctrl)
ipe.ctrl$tax_table %<>% tidy_taxonomy
ipe.ctrl$cal_abund()
t1 <- trans_network$new(dataset = ipe.ctrl, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
ipe.ctrl.csv = t1$res_network_attr
write.csv(ipe.ctrl.csv, "./Figures/Ipe_Ctrl_atributes.csv")
t1$save_network(filepath = "./Figures/Ipe_Control_network.gexf")


#Network for Handroanthus CS (10)
Ipe.cs <- subset_samples(input.final, Tree == "H. avellaneadae" & Substrate == "CS")
ipe.cs <- phyloseq2meco(Ipe.cs)
ipe.cs$tax_table %<>% tidy_taxonomy
ipe.cs$cal_abund()
t1 <- trans_network$new(dataset = ipe.cs, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
ipe.cs.csv = t1$res_network_attr
write.csv(ipe.cs.csv, "./Figures/Ipe_cs_atributes.csv")
t1$save_network(filepath = "./Figures/Ipe_Control_network.gexf")


#Network for Handroanthus CS + ADE (11)
Ipe.csade <- subset_samples(input.final, Tree == "H. avellaneadae" & Substrate == "CS + ADE")
ipe.csade <- phyloseq2meco(Ipe.csade)
ipe.csade$tax_table %<>% tidy_taxonomy
ipe.csade$cal_abund()
t1 <- trans_network$new(dataset = ipe.csade, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
ipe.csade.csv = t1$res_network_attr
write.csv(ipe.csade.csv, "./Figures/Ipe_CSade_atributes.csv")
t1$save_network(filepath = "./Figures/Ipe_CSade_network.gexf")


#Network for Handroanthus ADE (12)
Ipe.ade <- subset_samples(input.final, Tree == "H. avellaneadae" & Substrate == "ADE")
ipe.ade <- phyloseq2meco(Ipe.ade)
ipe.ade$tax_table %<>% tidy_taxonomy
ipe.ade$cal_abund()
t1 <- trans_network$new(dataset = ipe.ade, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
ipe.ade.csv = t1$res_network_attr
write.csv(ipe.ade.csv, "./Figures/Ipe_ADE_atributes.csv")
t1$save_network(filepath = "./Figures/Ipe_ADE_network.gexf")



###
#Network for Acacia Control (13)
Aca.ctrl <- subset_samples(input.final, Tree == "A. mangium" & Substrate == "Control")
aca.ctrl <- phyloseq2meco(Aca.ctrl)
aca.ctrl$tax_table %<>% tidy_taxonomy
aca.ctrl$cal_abund()
t1 <- trans_network$new(dataset = aca.ctrl, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
aca.ctrl.csv = t1$res_network_attr
write.csv(aca.ctrl.csv, "./Figures/Aca_Ctrl_atributes.csv")
t1$save_network(filepath = "./Figures/Aca_Control_network.gexf")


#Network for Acacia CS (14)
Aca.cs <- subset_samples(input.final, Tree == "A. mangium" & Substrate == "CS")
aca.cs <- phyloseq2meco(Aca.cs)
aca.cs$tax_table %<>% tidy_taxonomy
aca.cs$cal_abund()
t1 <- trans_network$new(dataset = aca.cs, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
aca.cs.csv = t1$res_network_attr
write.csv(aca.cs.csv, "./Figures/Aca_cs_atributes.csv")
t1$save_network(filepath = "./Figures/Aca_Control_network.gexf")


#Network for Acacia CS + ADE (15)
Aca.csade <- subset_samples(input.final, Tree == "A.mangium" & Substrate == "CS + ADE")
aca.csade <- phyloseq2meco(Aca.csade)
aca.csade$tax_table %<>% tidy_taxonomy
aca.csade$cal_abund()
t1 <- trans_network$new(dataset = aca.csade, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
aca.csade.csv = t1$res_network_attr
write.csv(aca.csade.csv, "./Figures/Aca_CSade_atributes.csv")
t1$save_network(filepath = "./Figures/Aca_CSade_network.gexf")


#Network for Handroanthus ADE (16)
Aca.ade <- subset_samples(input.final, Tree == "A. mangium" & Substrate == "ADE")
aca.ade <- phyloseq2meco(Aca.ade)
aca.ade$tax_table %<>% tidy_taxonomy
aca.ade$cal_abund()
t1 <- trans_network$new(dataset = aca.ade, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
aca.ade.csv = t1$res_network_attr
write.csv(aca.ade.csv, "./Figures/Aca_ADE_atributes.csv")
t1$save_network(filepath = "./Figures/Aca_ADE_network.gexf")



###
#Network for Urochloa Control (17)
Braq.ctrl <- subset_samples(input.final, Tree == "U. brizantha" & Substrate == "Control")
Braq.ctrl.gen <- aggregate_rare(Braq.ctrl, level = "Phylum",
                                         detection = 1/100, prevalence = 1/100)
braq.ctrl <- phyloseq2meco(Braq.ctrl.gen)
braq.ctrl$tax_table %<>% base::subset(Phylum != "p__")
braq.ctrl$tidy_dataset()
braq.ctrl$tax_table %<>% tidy_taxonomy
braq.ctrl$cal_abund()

t1 <- trans_network$new(dataset = braq.ctrl, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1 <- trans_network$new(dataset = braq.ctrl, 
                        cor_method = "spearman", 
                        use_WGCNA_pearson_spearman = TRUE, filter_thres = 0.0001)
pdf("./AAAAAAAAAAAAAAAAAAAAAAAA.pdf")
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
t1$res_network_attr
write.csv(braq.ctrl.csv, "./Figures/WGCNA_Braq_Ctrl_atributes.csv")
t1$save_network(filepath = "./Figures/WGCNA_Braq_Control_network.gexf")
t1$cal_sum_links(taxa_level = "Phylum")
t1$plot_sum_links(plot_pos = TRUE, plot_num = 50, color_values = ze)
svg('filename.svg')
# make plot
dev.off()



network.ctrl <- t1$res_network
my_color <- coul[as.numeric(as.factor(V(network.ctrl)$label))]
plot(network.ctrl, edge.color=my_color)



# create data:
links <- data.frame(
  source=c("A","A", "A", "A", "A","J", "B", "B", "C", "C", "D","I"),
  target=c("B","B", "C", "D", "J","A","E", "F", "G", "H", "I","I"),
  importance=(sample(1:4, 12, replace=T))
)
nodes <- data.frame(
  name=LETTERS[1:10],
  carac=c( rep("young",3),rep("adult",2), rep("old",5))
)
# Turn it into igraph object
network <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
# Make a palette of 3 colors
library(RColorBrewer)
coul  <- brewer.pal(3, "Set1") 
# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(network)$carac))]
# Make the plot
plot(network, vertex.color=my_color)
# Add a legend
legend("bottomleft", legend=levels(as.factor(V(network)$carac))  , col = coul , bty = "n", pch=20 , pt.cex = 3, cex = 1.5, text.col=coul , horiz = FALSE, inset = c(0.1, 0.1))

















#Network for Urochloa CS (18)
Braq.cs <- subset_samples(input.final, Tree == "U. brizantha" & Substrate == "CS")
braq.cs <- phyloseq2meco(Braq.cs)
braq.cs$tax_table %<>% tidy_taxonomy
braq.cs$cal_abund()
t1 <- trans_network$new(dataset = braq.cs, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
braq.cs.csv = t1$res_network_attr
write.csv(braq.cs.csv, "./Figures/Braq_cs_atributes.csv")
t1$save_network(filepath = "./Figures/Braq_Control_network.gexf")


#Network for Urochloa CS + ADE (19)
Braq.csade <- subset_samples(input.final, Tree == "U. brizantha" & Substrate == "CS + ADE")
braq.csade <- phyloseq2meco(Braq.csade)
braq.csade$tax_table %<>% tidy_taxonomy
braq.csade$cal_abund()
t1 <- trans_network$new(dataset = braq.csade, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
braq.csade.csv = t1$res_network_attr
write.csv(braq.csade.csv, "./Figures/Braq_CSade_atributes.csv")
t1$save_network(filepath = "./Figures/Braq_CSade_network.gexf")


#Network for Urochloa ADE (20)
Braq.ade <- subset_samples(input.final, Tree == "U. brizantha" & Substrate == "ADE")
braq.ade <- phyloseq2meco(Braq.ade)
braq.ade$tax_table %<>% tidy_taxonomy
braq.ade$cal_abund()
t1 <- trans_network$new(dataset = braq.ade, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
braq.ade.csv = t1$res_network_attr
write.csv(braq.ade.csv, "./Figures/Braq_ADE_atributes.csv")
t1$save_network(filepath = "./Figures/Braq_ADE_network.gexf")







### Network comparison for U. brizantha
Bri
Bri.ag.new <- aggregate_rare(Bri, level = "Genus",
                                              detection = 1/100, prevalence = 1/100)
braq.r6 <- phyloseq2meco(Bri.ag.new)
braq.r6$tax_table %<>% tidy_taxonomy
braq.r6$cal_abund()
soil_amp_network <- list()
tmp <- clone(braq.r6)
tmp$sample_table %<>% subset(Substrate == "CS")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)


###
#Network for Bulk Control (21)
Bulk.ctrl <- subset_samples(input.final, Tree == "Bulk" & Substrate == "Control")
bulk.ctrl <- phyloseq2meco(Bulk.ctrl)
bulk.ctrl$tax_table %<>% tidy_taxonomy
bulk.ctrl$cal_abund()
t1 <- trans_network$new(dataset = bulk.ctrl, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
bulk.ctrl.csv = t1$res_network_attr
write.csv(bulk.ctrl.csv, "./Figures/Bulk_Ctrl_atributes.csv")
t1$save_network(filepath = "./Figures/Bulk_Control_network.gexf")


#Network for Urochloa CS (22)
Bulk.cs <- subset_samples(input.final, Tree == "Bulk" & Substrate == "CS")
bulk.cs <- phyloseq2meco(Bulk.cs)
bulk.cs$tax_table %<>% tidy_taxonomy
bulk.cs$cal_abund()
t1 <- trans_network$new(dataset = bulk.cs, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
bulk.cs.csv = t1$res_network_attr
write.csv(bulk.cs.csv, "./Figures/Bulk_cs_atributes.csv")
t1$save_network(filepath = "./Figures/Bulk_Control_network.gexf")


#Network for Urochloa CS + ADE (23)
Bulk.csade <- subset_samples(input.final, Tree == "Bulk" & Substrate == "CS + ADE")
bulk.csade <- phyloseq2meco(Bulk.csade)
bulk.csade$tax_table %<>% tidy_taxonomy
bulk.csade$cal_abund()
t1 <- trans_network$new(dataset = bulk.csade, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
bulk.csade.csv = t1$res_network_attr
write.csv(bulk.csade.csv, "./Figures/Bulk_CSade_atributes.csv")
t1$save_network(filepath = "./Figures/Bulk_CSade_network.gexf")


#Network for Urochloa ADE (24)
Bulk.ade <- subset_samples(input.final, Tree == "Bulk" & Substrate == "ADE")
bulk.ade <- phyloseq2meco(Bulk.ade)
bulk.ade$tax_table %<>% tidy_taxonomy
bulk.ade$cal_abund()
t1 <- trans_network$new(dataset = bulk.ade, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.001, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
bulk.ade.csv = t1$res_network_attr
write.csv(bulk.ade.csv, "./Figures/Bulk_ADE_atributes.csv")
t1$save_network(filepath = "./Figures/Bulk_ADE_network.gexf")


##########################################################################


#### Alpha diversity
summarize_phyloseq(input.final)
richness <- estimate_richness(input.final)
head(richness)
head(sample_data(input.final))
alpha.all = cbind(richness, as.data.frame(sample_data(input.final)))

df5 <- data_summary(alpha.all, varname="Observed",
                    groupnames=c("Substrate", "Tree", "Treatment"))
head(df5)
df5$Substrate <- factor(df5$Substrate,
                             levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS")
                        )

alpha.plot = 
  ggplot(data = df5, aes(x=Tree, y=Observed, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Observed, ymax=Observed+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Tree", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = -15, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Acid Phosphatase Activity") + guides(color = "none")
alpha.plot

##########
#Alpha diversity Cecropia
summarize_phyloseq(Bri)
richness <- estimate_richness(Bri)
head(richness)
head(sample_data(Bri))
alpha.all = cbind(richness, as.data.frame(sample_data(Bri)))
aaa = agricolae::kruskal(y = alpha.all$Observed, trt = alpha.all$Treatment)
aaa = agricolae::kruskal(y = alpha.all$InvSimpson, trt = alpha.all$Treatment)
library(plyr)
df5 <- data_summary(alpha.all, varname="Observed",
                    groupnames=c("Substrate", "Tree", "Treatment"))
head(df5)
df5$Substrate <- factor(df5$Substrate,
                        levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS")
)

alpha.plot = 
  ggplot(data = df5, aes(x=Substrate, y=Observed, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Observed, ymax=Observed+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "Number of different taxa") +
  theme(axis.text.x = element_text(angle = -15, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Observed Diversity") + guides(color = "none") +
  theme_classic()
alpha.plot

rstatix::dunn_test(alpha.all, Observed ~ Treatment, p.adjust.method = "fdr")


################################
#Enzymes for Cecropia
library(readxl)
MapTrat <- read_excel("EDITED/MapTrat.xlsx")
map_cec_full = dplyr::filter(MapTrat, Tree =="Cecropia sp.")
map_cec_full[-c(1:4)] <- lapply(map_cec_full[-c(1:4)], function(x) as.numeric(as.character(x)))
rstatix::dunn_test(map_cec_full, Phosphatase ~ Treatment, p.adjust.method = "fdr")
rstatix::dunn_test(map_cec_full, BetaGlucosidase ~ Treatment, p.adjust.method = "fdr")
rstatix::dunn_test(map_cec_full, Arilsulfatase ~ Treatment, p.adjust.method = "fdr")
#no significance for anyone
df.pho.cec <- data_summary(map_cec_full, varname="Phosphatase",
                    groupnames=c("Substrate", "Tree", "Treatment"))
df.beta.cec <- data_summary(map_cec_full, varname="BetaGlucosidase",
                           groupnames=c("Substrate", "Tree", "Treatment"))
df.ari.cec <- data_summary(map_cec_full, varname="Arilsulfatase",
                           groupnames=c("Substrate", "Tree", "Treatment"))
df.pho.cec$Substrate <- factor(df.pho.cec$Substrate,
                        levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS"))
df.beta.cec$Substrate <- factor(df.beta.cec$Substrate,
                               levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS"))
df.ari.cec$Substrate <- factor(df.pho.cec$Substrate,
                               levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS"))

pho.plot = 
  ggplot(data = df.pho.cec, aes(x = Substrate, y = Phosphatase, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Phosphatase, ymax=Phosphatase+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Acid Phosphatase") + guides(color = "none") +
  theme_classic()
pho.plot


beta.plot = 
  ggplot(data = df.beta.cec, aes(x = Substrate, y = BetaGlucosidase, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=BetaGlucosidase, ymax=BetaGlucosidase+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "β-Glucosidase") + guides(color = "none") +
  theme_classic()
beta.plot


ari.plot = 
  ggplot(data = df.ari.cec, aes(x = Substrate, y = Arilsulfatase, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Arilsulfatase, ymax=Arilsulfatase+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Arylsulfatase") + guides(color = "none") +
  theme_classic()
ari.plot

ggpubr::ggarrange(pho.plot, beta.plot, ari.plot, common.legend = T, ncol = 3, nrow = 1)


###############
#Enzymes for Schizolobium
map_schi_full = dplyr::filter(MapTrat, Tree =="S. amazonicum")
map_schi_full[-c(1:4)] <- lapply(map_schi_full[-c(1:4)], function(x) as.numeric(as.character(x)))
rstatix::dunn_test(map_schi_full, Phosphatase ~ Treatment, p.adjust.method = "fdr")
rstatix::dunn_test(map_schi_full, BetaGlucosidase ~ Treatment, p.adjust.method = "fdr")
rstatix::dunn_test(map_schi_full, Arilsulfatase ~ Treatment, p.adjust.method = "fdr")
#no significance for anyone
df.pho.schi <- data_summary(map_schi_full, varname="Phosphatase",
                            groupnames=c("Substrate", "Tree", "Treatment"))
df.beta.schi <- data_summary(map_schi_full, varname="BetaGlucosidase",
                             groupnames=c("Substrate", "Tree", "Treatment"))
df.ari.schi <- data_summary(map_schi_full, varname="Arilsulfatase",
                            groupnames=c("Substrate", "Tree", "Treatment"))
df.pho.schi$Substrate <- factor(df.pho.schi$Substrate,
                                levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS"))
df.beta.schi$Substrate <- factor(df.beta.schi$Substrate,
                                 levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS"))
df.ari.schi$Substrate <- factor(df.pho.schi$Substrate,
                                levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS"))

pho.plot = 
  ggplot(data = df.pho.schi, aes(x = Substrate, y = Phosphatase, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Phosphatase, ymax=Phosphatase+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Acid Phosphatase") + guides(color = "none") +
  theme_classic()
pho.plot


beta.plot = 
  ggplot(data = df.beta.schi, aes(x = Substrate, y = BetaGlucosidase, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=BetaGlucosidase, ymax=BetaGlucosidase+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "β-Glucosidase") + guides(color = "none") +
  theme_classic()
beta.plot


ari.plot = 
  ggplot(data = df.ari.schi, aes(x = Substrate, y = Arilsulfatase, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Arilsulfatase, ymax=Arilsulfatase+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Arylsulfatase") + guides(color = "none") +
  theme_classic()
ari.plot

ggpubr::ggarrange(pho.plot, beta.plot, ari.plot, common.legend = T, ncol = 3, nrow = 1)


#########
#Enzymes for Handroanthus
map_ipe_full = dplyr::filter(MapTrat, Tree =="H. avellaneadae")
map_ipe_full[-c(1:4)] <- lapply(map_ipe_full[-c(1:4)], function(x) as.numeric(as.character(x)))
rstatix::dunn_test(map_ipe_full, Phosphatase ~ Treatment, p.adjust.method = "fdr")
rstatix::dunn_test(map_ipe_full, BetaGlucosidase ~ Treatment, p.adjust.method = "fdr")
rstatix::dunn_test(map_ipe_full, Arilsulfatase ~ Treatment, p.adjust.method = "fdr")
#no significance for anyone
df.pho.ipe <- data_summary(map_ipe_full, varname="Phosphatase",
                           groupnames=c("Substrate", "Tree", "Treatment"))
df.beta.ipe <- data_summary(map_ipe_full, varname="BetaGlucosidase",
                            groupnames=c("Substrate", "Tree", "Treatment"))
df.ari.ipe <- data_summary(map_ipe_full, varname="Arilsulfatase",
                           groupnames=c("Substrate", "Tree", "Treatment"))
df.pho.ipe$Substrate <- factor(df.pho.ipe$Substrate,
                               levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS"))
df.beta.ipe$Substrate <- factor(df.beta.ipe$Substrate,
                                levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS"))
df.ari.ipe$Substrate <- factor(df.pho.ipe$Substrate,
                               levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS"))

pho.plot = 
  ggplot(data = df.pho.ipe, aes(x = Substrate, y = Phosphatase, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Phosphatase, ymax=Phosphatase+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Acid Phosphatase") + guides(color = "none") +
  theme_classic()
pho.plot


beta.plot = 
  ggplot(data = df.beta.ipe, aes(x = Substrate, y = BetaGlucosidase, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=BetaGlucosidase, ymax=BetaGlucosidase+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "β-Glucosidase") + guides(color = "none") +
  theme_classic()
beta.plot


ari.plot = 
  ggplot(data = df.ari.ipe, aes(x = Substrate, y = Arilsulfatase, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Arilsulfatase, ymax=Arilsulfatase+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Arylsulfatase") + guides(color = "none") +
  theme_classic()
ari.plot

ggpubr::ggarrange(pho.plot, beta.plot, ari.plot, common.legend = T, ncol = 3, nrow = 1)

#############
#Enzymes for Acacia
map_aca_full = dplyr::filter(MapTrat, Tree =="A. mangium")
map_aca_full[-c(1:4)] <- lapply(map_aca_full[-c(1:4)], function(x) as.numeric(as.character(x)))
rstatix::dunn_test(map_aca_full, Phosphatase ~ Treatment, p.adjust.method = "fdr")
rstatix::dunn_test(map_aca_full, BetaGlucosidase ~ Treatment, p.adjust.method = "fdr")
rstatix::dunn_test(map_aca_full, Arilsulfatase ~ Treatment, p.adjust.method = "fdr")
#no significance for anyone
df.pho.aca <- data_summary(map_aca_full, varname="Phosphatase",
                           groupnames=c("Substrate", "Tree", "Treatment"))
df.beta.aca <- data_summary(map_aca_full, varname="BetaGlucosidase",
                            groupnames=c("Substrate", "Tree", "Treatment"))
df.ari.aca <- data_summary(map_aca_full, varname="Arilsulfatase",
                           groupnames=c("Substrate", "Tree", "Treatment"))
df.pho.aca$Substrate <- factor(df.pho.aca$Substrate,
                               levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS"))
df.beta.aca$Substrate <- factor(df.beta.aca$Substrate,
                                levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS"))
df.ari.aca$Substrate <- factor(df.pho.aca$Substrate,
                               levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS"))

pho.plot = 
  ggplot(data = df.pho.aca, aes(x = Substrate, y = Phosphatase, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Phosphatase, ymax=Phosphatase+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Acid Phosphatase") + guides(color = "none") +
  theme_classic()
pho.plot


beta.plot = 
  ggplot(data = df.beta.aca, aes(x = Substrate, y = BetaGlucosidase, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=BetaGlucosidase, ymax=BetaGlucosidase+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "β-Glucosidase") + guides(color = "none") +
  theme_classic()
beta.plot


ari.plot = 
  ggplot(data = df.ari.aca, aes(x = Substrate, y = Arilsulfatase, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Arilsulfatase, ymax=Arilsulfatase+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Arylsulfatase") + guides(color = "none") +
  theme_classic()
ari.plot

ggpubr::ggarrange(pho.plot, beta.plot, ari.plot, common.legend = T, ncol = 3, nrow = 1)

##########
#Enzymes for Urochloa
map_bri_full = dplyr::filter(MapTrat, Tree =="U. brizantha" & Substrate != "Initial")
map_bri_full[-c(1:4)] <- lapply(map_bri_full[-c(1:4)], function(x) as.numeric(as.character(x)))
rstatix::dunn_test(map_bri_full, Phosphatase ~ Treatment, p.adjust.method = "fdr")
rstatix::dunn_test(map_bri_full, BetaGlucosidase ~ Treatment, p.adjust.method = "fdr")
rstatix::dunn_test(map_bri_full, Arilsulfatase ~ Treatment, p.adjust.method = "fdr")
#no significance for anyone
df.pho.bri <- data_summary(map_bri_full, varname="Phosphatase",
                           groupnames=c("Substrate", "Tree", "Treatment"))
df.beta.bri <- data_summary(map_bri_full, varname="BetaGlucosidase",
                            groupnames=c("Substrate", "Tree", "Treatment"))
df.ari.bri <- data_summary(map_bri_full, varname="Arilsulfatase",
                           groupnames=c("Substrate", "Tree", "Treatment"))
df.pho.bri$Substrate <- factor(df.pho.bri$Substrate,
                               levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS"))
df.beta.bri$Substrate <- factor(df.beta.bri$Substrate,
                                levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS"))
df.ari.bri$Substrate <- factor(df.pho.bri$Substrate,
                               levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS"))

pho.plot = 
  ggplot(data = df.pho.bri, aes(x = Substrate, y = Phosphatase, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Phosphatase, ymax=Phosphatase+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Acid Phosphatase") + guides(color = "none") +
  theme_classic()
pho.plot


beta.plot = 
  ggplot(data = df.beta.bri, aes(x = Substrate, y = BetaGlucosidase, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=BetaGlucosidase, ymax=BetaGlucosidase+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "β-Glucosidase") + guides(color = "none") +
  theme_classic()
beta.plot


ari.plot = 
  ggplot(data = df.ari.bri, aes(x = Substrate, y = Arilsulfatase, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Arilsulfatase, ymax=Arilsulfatase+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Arylsulfatase") + guides(color = "none") +
  theme_classic()
ari.plot

ggpubr::ggarrange(pho.plot, beta.plot, ari.plot, common.legend = T, ncol = 3, nrow = 1)


###############
#Enzymes for Bulk
map_bulk_full = dplyr::filter(MapTrat, Tree =="None")
map_bulk_full[-c(1:4)] <- lapply(map_bulk_full[-c(1:4)], function(x) as.numeric(as.character(x)))
rstatix::dunn_test(map_bulk_full, Phosphatase ~ Treatment, p.adjust.method = "fdr")
rstatix::dunn_test(map_bulk_full, BetaGlucosidase ~ Treatment, p.adjust.method = "fdr")
rstatix::dunn_test(map_bulk_full, Arilsulfatase ~ Treatment, p.adjust.method = "fdr")
#no significance for anyone
df.pho.bulk <- data_summary(map_bulk_full, varname="Phosphatase",
                            groupnames=c("Substrate", "Tree", "Treatment"))
df.beta.bulk <- data_summary(map_bulk_full, varname="BetaGlucosidase",
                             groupnames=c("Substrate", "Tree", "Treatment"))
df.ari.bulk <- data_summary(map_bulk_full, varname="Arilsulfatase",
                            groupnames=c("Substrate", "Tree", "Treatment"))
df.pho.bulk$Substrate <- factor(df.pho.bulk$Substrate,
                                levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS"))
df.beta.bulk$Substrate <- factor(df.beta.bulk$Substrate,
                                 levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS"))
df.ari.bulk$Substrate <- factor(df.pho.bulk$Substrate,
                                levels = c("Initial" , "Control"  , "CS", "ADE", "ADE + CS"))

pho.plot = 
  ggplot(data = df.pho.bulk, aes(x = Substrate, y = Phosphatase, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Phosphatase, ymax=Phosphatase+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Acid Phosphatase") + guides(color = "none") +
  theme_classic()
pho.plot


beta.plot = 
  ggplot(data = df.beta.bulk, aes(x = Substrate, y = BetaGlucosidase, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=BetaGlucosidase, ymax=BetaGlucosidase+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "β-Glucosidase") + guides(color = "none") +
  theme_classic()
beta.plot


ari.plot = 
  ggplot(data = df.ari.bulk, aes(x = Substrate, y = Arilsulfatase, fill = Substrate)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Arilsulfatase, ymax=Arilsulfatase+sd), width=.2,
                position=position_dodge(.9))+
  labs(x = "Substrate", y= "ug of PNP/g soil/h") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment', title = "Arylsulfatase") + guides(color = "none") +
  theme_classic()
ari.plot

ggpubr::ggarrange(pho.plot, beta.plot, ari.plot, common.legend = T, ncol = 3, nrow = 1)


#Phyla abundance U. brizantha
Bri.t6 <- phyloseq2meco(Bri)

t1 <- trans_abund$new(dataset = Bri.t6, taxrank = "Phylum", ntaxa = 6, groupmean = "Substrate")
t1$plot_donut(label = FALSE)
t1$plot_donut(label = TRUE)

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


Bri   = subset_samples(input1, Tree == "U. brizantha")


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


### Alpha diversity
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


#Phyla abundance U. brizantha
Bri.t6 <- phyloseq2meco(Bri)

t1 <- trans_abund$new(dataset = Bri.t6, taxrank = "Phylum", ntaxa = 6, groupmean = "Substrate")
t1$plot_donut(label = FALSE)
t1$plot_donut(label = TRUE)

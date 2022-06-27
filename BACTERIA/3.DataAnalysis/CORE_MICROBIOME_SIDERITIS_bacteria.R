###############    Core microbiota


library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(dplyr) # data handling



####### Core microbiota anlaysis

##load phyloseq object for bacteria, (VINE_16S_All_2019.2020..RDS)##
VINE_16S_All_2019.2020. <- readRDS("../../2.PhyloseqObjectPrep/VINE_16S_All_2019.2020..RDS")



##SIDERITIS cultivar##
##Subset DataSet to cultivar Sideritis from different terroir units of Aigialeia's viticultural zone##
# Remove samples that don't belong to the cultivar Sideritis #

##subset samples##
Bacteria1920_vine_S <- subset_samples(Bacteria_1920_VineMicrobiome, !(Cultivar=="Agiorgitiko"))
Bacteria1920_vine_S <- subset_samples(Bacteria1920_vine_S, !(Cultivar=="Vidiano"))
Bacteria1920_vine_S <- subset_samples(Bacteria1920_vine_S, !(Cultivar=="Roditis"))
#Prune taxa because we did "subset samples' before
Bacteria1920_vine_S <- prune_taxa(taxa_sums(Bacteria1920_vine_S)>0,Bacteria1920_vine_S)


# Calculate compositional version of the data#
# (relative abundances) #
Bacteria1920_vine_S.rel <- microbiome::transform(Bacteria1920_vine_S, "compositional")
print(Bacteria1920_vine_S.rel)

##agglomerate in Genus level##
Bacteria1920_vine_S.rel_Gglomed<- tax_glom(Bacteria1920_vine_S.rel, taxrank = "Genus")
print(Bacteria1920_vine_S.rel_Gglomed)


#Check for the core ASVs
core.taxa.standard <- core_members(Bacteria1920_vine_S.rel_Gglomed, detection = 0.001, prevalence = 60/100)

print(core.taxa.standard)

saveRDS(Bacteria1920_vine_S.rel_Gglomed, file = "Bacteria1920_vine_S.rel_Gglomed.RDS")

#### Core heatmaps ####
# Core with compositionals
#library(RColorBrewer)
#library(reshape)
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-5), log10(.2), length = 10), 3)


# Make the plot and define the R colour palette: "RdBu" is been used
p.core <- plot_core(Bacteria1920_vine_S.rel_Gglomed, 
                    plot.type = "heatmap", 
                    colours = rev(brewer.pal(5, "RdBu")),
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = .6) + 
  xlab("Detection Threshold (Relative Abundance (%))")


p.core <- p.core + theme_bw() + ylab("ASVs")

print(p.core) 

dev.off()

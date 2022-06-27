###############    Core microbiota


library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(dplyr) # data handling



####### Core microbiota anlaysis
##load phyloseq object for fungi, (VINE_ITS_All_2019.2020..RDS)##
VINE_ITS_All_2019.2020. <- readRDS("../../2.PhyloseqObjectPrep/VINE_ITS_All_2019.2020..RDS")



##RODITIS cultivar##
##Subset DataSet to cultivar Roditis from different terroir units of Aigialeia's viticultural zone##
# Remove samples that don't belong to the cultivar Roditis #
##subset samples##
Fungi1920_vine_R <- subset_samples(VINE_ITS_All_2019.2020., !(Cultivar=="Agiorgitiko"))
Fungi1920_vine_R <- subset_samples(Fungi1920_vine_R, !(Cultivar=="Vidiano"))
Fungi1920_vine_R <- subset_samples(Fungi1920_vine_R, !(Cultivar=="Sideritis"))
#Prune taxa because we did "subset samples' before
Fungi1920_vine_R <- prune_taxa(taxa_sums(Fungi1920_vine_R)>0,Fungi1920_vine_R)


##save phyloseq object ##
saveRDS(Fungi1920_vine_R, file = "Fungi1920_vine_R.RDS")

# Calculate compositional version of the data#
# (relative abundances) #
Fungi1920_vine_R.rel <- microbiome::transform(Fungi1920_vine_R, "compositional")
print(Fungi1920_vine_R.rel)

##agglomerate in Genus level##
Fungi1920_vine_R.rel_Gglomed<- tax_glom(Fungi1920_vine_R.rel, taxrank = "Genus")
print(Fungi1920_vine_R.rel_Gglomed)


#Check for the core ASVs
core.taxa.standard <- core_members(Fungi1920_vine_R.rel_Gglomed, detection = 0.001, prevalence = 93/100)

print(core.taxa.standard)

saveRDS(Fungi1920_vine_R.rel_Gglomed, file = "Fungi1920_vine_R.rel_Gglomed.RDS")

#### Core heatmaps ####
# Core with compositionals
#library(RColorBrewer)
#library(reshape)
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-5), log10(.2), length = 10), 3)


# Make the plot and define the R colour palette: "RdBu" is been used
p.core <- plot_core(Fungi1920_vine_R.rel_Gglomed, 
                    plot.type = "heatmap", 
                    colours = rev(brewer.pal(5, "RdBu")),
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = .93) + 
  xlab("Detection Threshold (Relative Abundance (%))")


p.core <- p.core + theme_bw() + ylab("ASVs")

print(p.core) 

dev.off()

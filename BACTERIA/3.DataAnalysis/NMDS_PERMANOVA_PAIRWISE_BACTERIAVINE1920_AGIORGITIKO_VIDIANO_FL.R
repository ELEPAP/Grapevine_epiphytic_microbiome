##NMDS Analysis##
##load phyloseq object for bacteria, (VINE_16S_All_2019.2020..RDS)##
##NMDS must be applied for varieties Vidiano and Agiorgitiko separately from  varieties Roditis and Sideritis  ##

VINE_16S_All_2019.2020. <- readRDS("../../2.PhyloseqObjectPrep/VINE_16S_All_2019.2020..RDS")


##Subset DataSet to varieties Agiorgitiko and Vidiano from different viticultural zones  and perform NMDS per Viticultural_Zones, Vintage and cultivar##
# Remove samples that don't belong to the cultivars Agiorgitiko and Vidiano #

##AGIORGITIKO & VIDIANO##

##subset samples##
Bacteria1920_vine_AV <- subset_samples(VINE_16S_All_2019.2020., !(Cultivar=="Roditis"))
Bacteria1920_vine_AV <- subset_samples(Bacteria1920_vine_AV, !(Cultivar=="Sideritis"))
Bacteria1920_vine_AV <- prune_taxa(taxa_sums(Bacteria1920_vine_AV)>0,Bacteria1920_vine_AV)

##save phyloseq object and continue with the DA analysis scripts##
saveRDS(Bacteria1920_vine_AV, file = "Bacteria1920_vine_AV.RDS")

##transform phyloseq object raw counts to relative abundance (100%)##
Bacteria1920_vine_AV_100 <- transform_sample_counts(Bacteria1920_vine_AV, function(OTU) 100*OTU/sum(OTU))

##calculate bray-curtis distance##
ord.nmds.bray.Bacteria1920_vine_AV_100 <- ordinate(Bacteria1920_vine_AV_100, method="NMDS", distance="bray")

##plot the NMDS for cultivars Agiorgitiko and Vidiano (colored per Viticultural_Zones  firstly, secondly per Vintage and thirdly per Cultivar)##
plot_ordination(Bacteria1920_vine_AV_100, ord.nmds.bray.Bacteria1920_vine_AV_100, color="Viticultural_Zones",label = "", title=paste("NMDS (stress ",round(ord.nmds.bray.Bacteria1920_vine_AV_100$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("blue3","red")) + theme_bw()

plot_ordination(Bacteria1920_vine_AV_100, ord.nmds.bray.Bacteria1920_vine_AV_100, color="Vintage",label = "", title=paste("NMDS (stress ",round(ord.nmds.bray.Bacteria1920_vine_AV_100$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("blue3","red")) + theme_bw()

plot_ordination(Bacteria1920_vine_AV_100, ord.nmds.bray.Bacteria1920_vine_AV_100, color="Cultivar",label = "", title=paste("NMDS (stress ",round(ord.nmds.bray.Bacteria1920_vine_AV_100$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("blue3","red")) + theme_bw()

##pairwise permanova to groups of NMDS###
##check if pairwise permanova between variables (Viticultural_Zones,Cultivar, Vintage) differ significally and use p-value to manuscript data##
##load package pairwiseAdonis##
library(pairwiseAdonis)

##Viticultural_Zones vs Vintage##
mycmpfactorBacteriaAV <- interaction(data.frame(Bacteria1920_vine_AV_100@sam_data)$Viticultural_Zones, data.frame(Bacteria1920_vine_AV_100@sam_data)$Vintage)

mympairwisepermBacteriaAV <- pairwise.adonis(Bacteria1920_vine_AV_100@otu_table, mycmpfactorBacteriaAV, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

##Viticultural_Zones vs Cultivar##
mycmpfactorBacteriaAV_1 <- interaction(data.frame(Bacteria1920_vine_AV_100@sam_data)$Viticultural_Zones, data.frame(Bacteria1920_vine_AV_100@sam_data)$Cultivar)

mympairwisepermBacteriaAV_1 <- pairwise.adonis(Bacteria1920_vine_AV_100@otu_table, mycmpfactorBacteriaAV_1, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

##Cultivar vs Vintage##
mycmpfactorBacteriaAV_2 <- interaction(data.frame(Bacteria1920_vine_AV_100@sam_data)$Cultivar, data.frame(Bacteria1920_vine_AV_100@sam_data)$Vintage)

mympairwisepermBacteriaAV_2 <- pairwise.adonis(Bacteria1920_vine_AV_100@otu_table, mycmpfactorBacteriaAV_2, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")
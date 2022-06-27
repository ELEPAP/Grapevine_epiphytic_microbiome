##NMDS Analysis##

##load phyloseq object for fungi, (VINE_ITS_All_2019.2020..RDS)##
##NMDS must be applied for varieties Vidiano and Agiorgitiko separately from  varieties Roditis and Sideritis  ##

VINE_ITS_All_2019.2020. <- readRDS("../../2.PhyloseqObjectPrep/VINE_ITS_All_2019.2020..RDS")


##Subset DataSet to varieties Agiorgitiko and Vidiano from different viticultural zones and perform NMDS per Viticultural_Zones, Vintage and Cultivar##
# Remove samples that don't belong to the cultivars Agiorgitiko and Vidiano #

##AGIORGITIKO & VIDIANO##

##subset samples##
Fungi1920_vine_AV <- subset_samples(VINE_ITS_All_2019.2020., !(variety=="Roditis"))
Fungi1920_vine_AV <- subset_samples(Fungi1920_vine_AV, !(variety=="Sideritis"))
Fungi1920_vine_AV <- prune_taxa(taxa_sums(Fungi1920_vine_AV)>0,Fungi1920_vine_AV)

##transform phyloseq object raw counts to relative abundance (100%)##
Fungi1920_vine_AV100 <- transform_sample_counts(Fungi1920_vine_AV, function(OTU) 100*OTU/sum(OTU))

##calculate bray-curtis distance##
ord.nmds.bray.Fungi1920_vine_AV100 <- ordinate(Fungi1920_vine_AV100, method="NMDS", distance="bray")

##plot the NMDS for varieties Agiorgitiko and Vidiano (colored per Biogeography condition firstly, secondly per Vintage and thirdly per cultivar)##
plot_ordination(Fungi1920_vine_AV100, ord.nmds.bray.Fungi1920_vine_AV100, color="Viticultural_Zones",label = "", title=paste("NMDS (stress ",round(ord.nmds.bray.Fungi1920_vine_AV100$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("blue3","red")) + theme_bw()

plot_ordination(Fungi1920_vine_AV100, ord.nmds.bray.Fungi1920_vine_AV100, color="Vintage",label = "", title=paste("NMDS (stress ",round(ord.nmds.bray.Fungi1920_vine_AV100$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("blue3","red")) + theme_bw()

plot_ordination(Fungi1920_vine_AV100, ord.nmds.bray.Fungi1920_vine_AV100, color="Cultivar",label = "", title=paste("NMDS (stress ",round(ord.nmds.bray.Fungi1920_vine_AV100$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("blue3","red")) + theme_bw()


##pairwise permanova to groups of NMDS###
##check if pairwise permanova between variables (Viticultural_Zones, Vintage,  Cultivar) differ significantly and use p-value to manuscript data##
##load package pairwiseAdonis##
library(pairwiseAdonis)

##Viticultural_Zones vs Vintage##
mycmpfactorFungiAV <- interaction(data.frame(Fungi1920_vine_AV100@sam_data)$Viticultural_Zones, data.frame(Fungi1920_vine_AV100@sam_data)$Vintage)

mympairwisepermFungiAV <- pairwise.adonis(Fungi1920_vine_AV100@otu_table, mycmpfactorFungiAV, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

##Viticultural_Zones vs Cultivar##
mycmpfactorFungiAV1 <- interaction(data.frame(Fungi1920_vine_AV100@sam_data)$Viticultural_Zones, data.frame(Fungi1920_vine_AV100@sam_data)$Cultivar)

mympairwisepermFungiAV1 <- pairwise.adonis(Fungi1920_vine_AV100@otu_table, mycmpfactorFungiAV1, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

##Vintage vs cultivar##
mycmpfactorFungiAV2 <- interaction(data.frame(Fungi1920_vine_AV100@sam_data)$Vintage, data.frame(Fungi1920_vine_AV100@sam_data)$Cultivar)

mympairwisepermFungiAV2 <- pairwise.adonis(Fungi1920_vine_AV100@otu_table, mycmpfactorFungiAV2, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

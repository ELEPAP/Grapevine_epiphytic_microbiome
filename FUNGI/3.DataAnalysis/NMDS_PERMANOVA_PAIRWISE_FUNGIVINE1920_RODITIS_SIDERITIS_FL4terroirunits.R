##NMDS Analysis##

##load phyloseq object for fungi, (VINE_ITS_All_2019.2020..RDS)##
##NMDS must be applied for varieties Vidiano and Agiorgitiko separately from  varieties Roditis and Sideritis  ##

VINE_ITS_All_2019.2020. <- readRDS("../../2.PhyloseqObjectPrep/VINE_ITS_All_2019.2020..RDS")


##Subset DataSet to cultivars Roditis and Sideritis from different terroir units of Aigialeia's viticultural zone and perform NMDS per terroir_unit, Vintage and cultivar##
# Remove samples that don't belong to the cultivars Roditis and Sideritis #

##RODITIS & SIDERITIS##

##subset samples##
Fungi1920_vine_RS <- subset_samples(VINE_ITS_All_2019.2020., !(variety=="Agiorgitiko"))
Fungi1920_vine_RS <- subset_samples(Fungi1920_vine_RS, !(variety=="Vidiano"))
Fungi1920_vine_RS <- prune_taxa(taxa_sums(Fungi1920_vine_RS)>0,Fungi1920_vine_RS)


##transform phyloseq object raw counts to relative abundance (100%)##
Fungi1920_vine_RS_100 <- transform_sample_counts(Fungi1920_vine_RS, function(OTU) 100*OTU/sum(OTU))

##calculate bray-curtis distance##
ord.nmds.bray.Fungi1920_vine_RS_100 <- ordinate(Fungi1920_vine_RS_100, method="NMDS", distance="bray")

##plot the NMDS for cultivars Roditis and Sideritis (colored per Cultivar firstly and secondly per Vintage)##
plot_ordination(Fungi1920_vine_RS_100, ord.nmds.bray.Fungi1920_vine_RS_100, color="Cultivar",label = "", title=paste("NMDS (stress ",round(ord.nmds.bray.Fungi1920_vine_RS_100$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("blue3","red")) + theme_bw()

plot_ordination(Fungi1920_vine_RS_100, ord.nmds.bray.Fungi1920_vine_RS_100, color="Vintage",label = "", title=paste("NMDS (stress ",round(ord.nmds.bray.Fungi1920_vine_RS_100$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("blue3","red")) + theme_bw()


##pairwise permanova to groups of NMDS###
##check if pairwise permanova between variables (cultivar, Vintage) differ significantly and use p-value to manuscript data##
##load package pairwiseAdonis##
library(pairwiseAdonis)

##Cultivar vs Vintage##
mycmpfactorFungiRS <- interaction(data.frame(Fungi1920_vine_RS_100@sam_data)$Cultivar, data.frame(Fungi1920_vine_RS_100@sam_data)$Vintage)

mympairwisepermFungiRS <- pairwise.adonis(Fungi1920_vine_RS_100@otu_table, mycmpfactorFungiRS, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")


##EACH CULTIVAR SEPARATELY##
##Subset DataSet to each cultivar and perform NMDS only per Terroir_unit##
##Roditis##
##subset samples##
Fungi1920_vine_R<- subset_samples(Fungi1920_vine_RS, !(Cultivar=="Sideritis"))
Fungi1920_vine_R <- prune_taxa(taxa_sums(Fungi1920_vine_R)>0,Fungi1920_vine_R)

##Roditis
##transform phyloseq object raw counts to relative abundance (100%)##
Fungi1920_vine_R_100 <- transform_sample_counts(Fungi1920_vine_R, function(OTU) 100*OTU/sum(OTU))

##Roditis
##calculate bray-curtis distance##
ord.nmds.bray.Fungi1920_vine_R <- ordinate(Fungi1920_vine_R_100, method="NMDS", distance="bray")

##plot the NMDS for cultivar Roditis (colored per Terroir_Unit)##
plot_ordination(Fungi1920_vine_R_100, ord.nmds.bray.Fungi1920_vine_R, color="Terroir_Unit",label = "", title=paste("NMDS (stress ",round(ord.nmds.bray.Fungi1920_vine_R$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("brown2","darkorchid1")) + theme_bw()

##plot the NMDS for cultivar Roditis (colored per Vintage)##
plot_ordination(Fungi1920_vine_R_100, ord.nmds.bray.Fungi1920_vine_R, color="Vintage",label = "", title=paste("NMDS (stress ",round(ord.nmds.bray.Fungi1920_vine_R$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("brown2","darkorchid1")) + theme_bw()


##Sideritis##
##subset samples##
Fungi1920_vine_S<- subset_samples(Fungi1920_vine_RS, !(Cultivar=="Roditis"))
Fungi1920_vine_S <- prune_taxa(taxa_sums(Fungi1920_vine_S)>0,Fungi1920_vine_S)

##Sideritis
##transform phyloseq object raw counts to relative abundance (100%)##
Fungi1920_vine_S_100 <- transform_sample_counts(Fungi1920_vine_S, function(OTU) 100*OTU/sum(OTU))

##Sideritis
##calculate bray-curtis distance##
ord.nmds.bray.Fungi1920_vine_S <- ordinate(Fungi1920_vine_S_100, method="NMDS", distance="bray")

##plot the NMDS for culticar Sideritis (colored per Terroir_Unit)##
plot_ordination(Fungi1920_vine_S_100, ord.nmds.bray.Fungi1920_vine_S, color="Terroir_Unit",label = "", title=paste("NMDS (stress ",round(ord.nmds.bray.Fungi1920_vine_S$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("brown2","darkorchid1")) + theme_bw()

##plot the NMDS for culticar Sideritis (colored per Vintage)##
plot_ordination(Fungi1920_vine_S_100, ord.nmds.bray.Fungi1920_vine_S, color="Vintage",label = "", title=paste("NMDS (stress ",round(ord.nmds.bray.Fungi1920_vine_S$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("brown2","darkorchid1")) + theme_bw()

##pairwise permanova to groups of NMDS###

##check if pairwise permanova between variables (Terroir_Unit, Vintage) differ significally and use p-value to manuscript data##

##Roditis##
##Terroir_Unit vs Vintage##
mycmpfactorRoditis <- interaction(data.frame(Fungi1920_vine_R_100@sam_data)$Terroir_Unit, data.frame(Fungi1920_vine_R_100@sam_data)$Vintage)

mympairwisepermRoditis <- pairwise.adonis(Fungi1920_vine_R_100@otu_table, mycmpfactorRoditis, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

##Sideritis##
##Terroir_Unit vs Vintage##
mycmpfactorSideritis <- interaction(data.frame(Fungi1920_vine_S_100@sam_data)$Terroir_Unit, data.frame(Fungi1920_vine_S_100@sam_data)$Vintage)

mympairwisepermSideritis <- pairwise.adonis(Fungi1920_vine_S_100@otu_table, mycmpfactorSideritis, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")
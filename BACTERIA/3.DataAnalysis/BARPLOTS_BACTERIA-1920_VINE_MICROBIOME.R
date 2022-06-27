##Bar Plots##

##load phyloseq object for bacteria, (VINE_16S_All_2019.2020..RDS)##
##Barplots must be constructed for cultivars Vidiano and Agiorgitiko separately from cultivars Roditis and Sideritis and for each group of cultivars every tissue will be separately (fruits and leaves)##

VINE_16S_All_2019.2020. <- readRDS("../../2.PhyloseqObjectPrep/VINE_16S_All_2019.2020..RDS")


##Subset DataSet to cultivars Agiorgitiko and Vidiano from different viticultural zones and construct Barplots separately for each plant tissue (fruits and leaves)##
# Remove samples that don't belong to the cultivars Agiorgitiko and Vidiano #

##AGIORGITIKO & VIDIANO##
##subset samples##
Bacteria1920_vine_AV <- subset_samples(Bacteria_1920_VineMicrobiome, !(Cultivar=="Roditis"))
Bacteria1920_vine_AV <- subset_samples(Bacteria1920_vine_AV, !(Cultivar=="Sideritis"))
Bacteria1920_vine_AV <- prune_taxa(taxa_sums(Bacteria1920_vine_AV)>0,Bacteria1920_vine_AV)

##AGIORGITIKO & VIDIANO - FRUITS ##
# subset samples-Keep Fruits#
Bacteria1920_vine_AV_F <- subset_samples(Bacteria1920_vine_AV, !(plant_part=="Leaves"))
Bacteria1920_vine_AV_F <- prune_taxa(taxa_sums(Bacteria1920_vine_AV_F)>0,Bacteria1920_vine_AV_F)

##check sample_data after separating Fruits##
sample_data(Bacteria1920_vine_AV_F)

##transform phyloseq object raw counts to relative abundance (100%)##
Bacteria1920_vine_AV_F_100 <- transform_sample_counts(Bacteria1920_vine_AV_F, function(OTU) 100*OTU/sum(OTU))

##agglomerate in Phylum level##
Bacteria1920_vine_AV_F_100GlomPhylum <-tax_glom(Bacteria1920_vine_AV_F_100, taxrank = "Phylum")

##save phyloseq object and continue with the statistical analysis scripts##
saveRDS(Bacteria1920_vine_AV_F_100GlomPhylum, file = "Bacteria1920_vine_AV_F_100GlomPhylum.RDS")

##subset Phylum to the desirable taxa for plotting,in our case Actinobacteriota, Bacteroidota, Chloroflexi, Cyanobacteria, Firmicutes, Myxococcota, Planctomycetota, Proteobacteria. Average > 1% calculated in Excel##
Bacteria1920_vine_AV_F_100GlomPhylumSelection <- prune_taxa(row.names(tax_table(Bacteria1920_vine_AV_F_100GlomPhylum)[which(tax_table(Bacteria1920_vine_AV_F_100GlomPhylum)[,"Phylum"]%in%c("Actinobacteriota","Bacteroidota","Chloroflexi","Cyanobacteria","Firmicutes","Myxococcota","Planctomycetota","Proteobacteria"))]), Bacteria1920_vine_AV_F_100GlomPhylum)

##Plot Bars for bacteria in selected taxa## 
plot_bar(Bacteria1920_vine_AV_F_100GlomPhylumSelection, x="comb_num", fill="Phylum", facet_grid ="Vintage") + geom_col()

##AGIORGITIKO & VIDIANO - LEAVES ##
# subset samples-Keep Leaves#
Bacteria1920_vine_AV_L <- subset_samples(Bacteria1920_vine_AV, !(plant_part=="Fruit"))
Bacteria1920_vine_AV_L <- prune_taxa(taxa_sums(Bacteria1920_vine_AV_L)>0,Bacteria1920_vine_AV_L)

##check sample_data after separating Fruits##
sample_data(Bacteria1920_vine_AV_L)

##transform phyloseq object raw counts to relative abundance (100%)##
Bacteria1920_vine_AV_L_100 <- transform_sample_counts(Bacteria1920_vine_AV_L, function(OTU) 100*OTU/sum(OTU))

##agglomerate in Phylum level##
Bacteria1920_vine_AV_L_100GlomPhylum <-tax_glom(Bacteria1920_vine_AV_L_100, taxrank = "Phylum")

##save phyloseq object and continue with the statistical analysis scripts##
saveRDS(Bacteria1920_vine_AV_L_100GlomPhylum, file = "Bacteria1920_vine_AV_L_100GlomPhylum.RDS")

##subset Phylum to the desirable taxa for plotting,in our case Actinobacteriota, Bacteroidota, Chloroflexi, Cyanobacteria, Firmicutes, Myxococcota, Planctomycetota, Proteobacteria. Average > 1% calculated in Excel##
Bacteria1920_vine_AV_L_100GlomPhylumSelection <- prune_taxa(row.names(tax_table(Bacteria1920_vine_AV_L_100GlomPhylum)[which(tax_table(Bacteria1920_vine_AV_L_100GlomPhylum)[,"Phylum"]%in%c("Actinobacteriota","Bacteroidota","Chloroflexi","Cyanobacteria","Firmicutes","Myxococcota","Planctomycetota","Proteobacteria"))]), Bacteria1920_vine_AV_L_100GlomPhylum)

##Plot Bars for bacteria in selected taxa## 
plot_bar(Bacteria1920_vine_AV_L_100GlomPhylumSelection, x="comb_num", fill="Phylum", facet_grid ="Vintage") + geom_col()



##RODITIS & SIDERITIS##
##Subset DataSet to cultivars Roditis and Sideritis from different terroir units of Aigialeia's viticultural zone and construct Barplots separately for each plant tissue(fruits and leaves)##
# Remove samples that don't belong to the cultivars Roditis and Sideritis #

##subset samples##
Bacteria1920_vine_RS <- subset_samples(Bacteria_1920_VineMicrobiome, !(Cultivar=="Agiorgitiko"))
Bacteria1920_vine_RS <- subset_samples(Bacteria1920_vine_RS, !(Cultivar=="Vidiano"))
Bacteria1920_vine_RS <- prune_taxa(taxa_sums(Bacteria1920_vine_RS)>0,Bacteria1920_vine_RS)


##save phyloseq object ##
saveRDS(Bacteria1920_vine_RS, file = "Bacteria1920_vine_RS.RDS")

##RODITIS & SIDERITIS- FRUITS ##
# subset samples-Keep Fruits#
Bacteria1920_vine_RS_F <- subset_samples(Bacteria1920_vine_RS_fl4, !(plant_part=="Leaves"))
Bacteria1920_vine_RS_F <- prune_taxa(taxa_sums(Bacteria1920_vine_RS_F)>0,Bacteria1920_vine_RS_F)

##check sample_data after separating Fruits##
sample_data(Bacteria1920_vine_RS_F)

##transform phyloseq object raw counts to relative abundance (100%)##
Bacteria1920_vine_RS_F_100 <- transform_sample_counts(Bacteria1920_vine_RS_F, function(OTU) 100*OTU/sum(OTU))

##agglomerate in Phylum level##
Bacteria1920_vine_RS_F_100GlomPhylum <-tax_glom(Bacteria1920_vine_RS_F_100, taxrank = "Phylum")

##save phyloseq object and continue with the statistical analysis scripts##
saveRDS(Bacteria1920_vine_RS_F_100GlomPhylum, file = "Bacteria1920_vine_RS_F_100GlomPhylum.RDS")

##subset Phylum to the desirable taxa for plotting,in our case Acidobacteriota, Actinobacteriota, Bacteroidota, Chloroflexi, Cyanobacteria, Firmicutes, Myxococcota, Proteobacteria. Average > 1% calculated in Excel##
Bacteria1920_vine_RS_F_100GlomPhylumSelection <- prune_taxa(row.names(tax_table(Bacteria1920_vine_RS_F_100GlomPhylum)[which(tax_table(Bacteria1920_vine_RS_F_100GlomPhylum)[,"Phylum"]%in%c("Acidobacteriota","Actinobacteriota","Bacteroidota","Chloroflexi","Cyanobacteria","Firmicutes","Myxococcota","Proteobacteria"))]), Bacteria1920_vine_RS_F_100GlomPhylum)

##Plot Bars for bacteria in selected taxa## 
plot_bar(Bacteria1920_vine_RS_F_100GlomPhylumSelection, x="comb_num", fill="Phylum", facet_grid ="Vintage") + geom_col()

##RODITIS & SIDERITIS- LEAVES ##
# subset samples-Keep Leaves#
Bacteria1920_vine_RS_L <- subset_samples(Bacteria1920_vine_RS_fl4, !(plant_part=="Fruit"))
Bacteria1920_vine_RS_L <- prune_taxa(taxa_sums(Bacteria1920_vine_RS_L)>0,Bacteria1920_vine_RS_L)

##check sample_data after separating Fruits##
sample_data(Bacteria1920_vine_RS_L)

##transform phyloseq object raw counts to relative abundance (100%)##
Bacteria1920_vine_RS_L_100 <- transform_sample_counts(Bacteria1920_vine_RS_L, function(OTU) 100*OTU/sum(OTU))

##agglomerate in Phylum level##
Bacteria1920_vine_RS_L_100GlomPhylum <-tax_glom(Bacteria1920_vine_RS_L_100, taxrank = "Phylum")

##save phyloseq object and continue with the statistical analysis scripts##
saveRDS(Bacteria1920_vine_RS_L_100GlomPhylum, file = "Bacteria1920_vine_RS_L_100GlomPhylum.RDS")

##subset Phylum to the desirable taxa for plotting,in our case Acidobacteriota,Actinobacteriota, Bacteroidota, Chloroflexi, Cyanobacteria, Firmicutes, Myxococcota, Proteobacteria. Average > 1% calculated in Excel##
Bacteria1920_vine_RS_L_100GlomPhylumSelection <- prune_taxa(row.names(tax_table(Bacteria1920_vine_RS_L_100GlomPhylum)[which(tax_table(Bacteria1920_vine_RS_L_100GlomPhylum)[,"Phylum"]%in%c("Acidobacteriota","Actinobacteriota","Bacteroidota","Chloroflexi","Cyanobacteria","Firmicutes","Myxococcota","Proteobacteria"))]), Bacteria1920_vine_RS_L_100GlomPhylum)

##Plot Bars for fungi in selected taxa## 
plot_bar(Bacteria1920_vine_RS_L_100GlomPhylumSelection, x="comb_num", fill="Phylum", facet_grid ="Vintage") + geom_col()
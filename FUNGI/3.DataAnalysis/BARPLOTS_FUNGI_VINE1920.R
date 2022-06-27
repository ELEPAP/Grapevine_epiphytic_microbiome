##Bar Plots##

##load phyloseq object for fungi, (VINE_ITS_All_2019.2020..RDS)##
##Barplots must be constructed for cultivars Vidiano and Agiorgitiko separately from cultivars Roditis and Sideritis and for each group of cultivars every tissue will be separately, too (fruits and leaves)##

VINE_ITS_All_2019.2020. <- readRDS("../../2.PhyloseqObjectPrep/VINE_ITS_All_2019.2020.RDS")



##Subset DataSet to cultivars Agiorgitiko and Vidiano from different viticultural zones and construct Barplots separately for each plant_part (fruits and leaves)##
# Remove samples that don't belong to the cultivars Agiorgitiko and Vidiano #

##AGIORGITIKO & VIDIANO##
##subset samples##
Fungi1920_vine_AV <- subset_samples(VINE_ITS_All_2019.2020., !(Cultivar=="Roditis"))
Fungi1920_vine_AV <- subset_samples(Fungi1920_vine_AV, !(Cultivar=="Sideritis"))
Fungi1920_vine_AV <- prune_taxa(taxa_sums(Fungi1920_vine_AV)>0,Fungi1920_vine_AV)

##AGIORGITIKO & VIDIANO - FRUITS ##
# subset samples-Keep Fruits#
Fungi1920_vine_AV_F <- subset_samples(Fungi1920_vine_AV, !(plant_part=="Leaves"))
Fungi1920_vine_AV_F <- prune_taxa(taxa_sums(Fungi1920_vine_AV_F)>0,Fungi1920_vine_AV_F)

##check sample_data after separating Fruits##
sample_data(Fungi1920_vine_AV_F)

##transform phyloseq object raw counts to relative abundance (100%)##
Fungi1920_vine_AV_F_100 <- transform_sample_counts(Fungi1920_vine_AV_F, function(OTU) 100*OTU/sum(OTU))

##agglomerate in Class level##
Fungi1920_vine_AV_F_100GlomClass <-tax_glom(Fungi1920_vine_AV_F_100, taxrank = "Class")

##save phyloseq object and continue with the statistical analysis scripts##
saveRDS(Fungi1920_vine_AV_F_100GlomClass, file = "Fungi1920_vine_AV_F_100GlomClass.RDS")

##subset Class to the desirable taxa for plotting,in our case Dothideomycetes, Eurotiomycetes, Leotiomycetes, Microbotryomycetes, Saccharomycetes, Sordariomycetes, Tremellomycetes, Ustilaginomycetes. Average > 1% calculated in Excel##
Fungi1920_vine_AV_F_100GlomClassSelection <- prune_taxa(row.names(tax_table(Fungi1920_vine_AV_F_100GlomClass)[which(tax_table(Fungi1920_vine_AV_F_100GlomClass)[,"Class"]%in%c("c__Dothideomycetes","c__Eurotiomycetes","c__Leotiomycetes","c__Microbotryomycetes","c__Saccharomycetes","c__Sordariomycetes","Tremellomycetes","Ustilaginomycetes."))]), Fungi1920_vine_AV_F_100GlomClass)

##Plot Bars for fungi in selected taxa## 
plot_bar(Fungi1920_vine_AV_F_100GlomClassSelection, x="comb_num", fill="Class", facet_grid ="Vintage") + geom_col()

##AGIORGITIKO & VIDIANO - LEAVES ##
# subset samples-Keep Leaves#
Fungi1920_vine_AV_L <- subset_samples(Fungi1920_vine_AV, !(plant_part=="Fruit"))
Fungi1920_vine_AV_L <- prune_taxa(taxa_sums(Fungi1920_vine_AV_L)>0,Fungi1920_vine_AV_L)

##check sample_data after separating Fruits##
sample_data(Fungi1920_vine_AV_L)

##transform phyloseq object raw counts to relative abundance (100%)##
Fungi1920_vine_AV_L_100 <- transform_sample_counts(Fungi1920_vine_AV_L, function(OTU) 100*OTU/sum(OTU))

##agglomerate in Class level##
Fungi1920_vine_AV_L_100GlomClass <-tax_glom(Fungi1920_vine_AV_L_100, taxrank = "Class")

##save phyloseq object and continue with the statistical analysis scripts##
saveRDS(Fungi1920_vine_AV_L_100GlomClass, file = "Fungi1920_vine_AV_L_100GlomClass.RDS")

##subset Class to the desirable taxa for plotting,in our case Dothideomycetes, Eurotiomycetes, Leotiomycetes, Malasseziomycetes, Saccharomycetes, Sordariomycetes, Tremellomycetes, Ustilaginomycetes. Average > 1% calculated in Excel##
Fungi1920_vine_AV_L_100GlomClassSelection <- prune_taxa(row.names(tax_table(Fungi1920_vine_AV_L_100GlomClass)[which(tax_table(Fungi1920_vine_AV_L_100GlomClass)[,"Class"]%in%c("c__Dothideomycetes","c__Eurotiomycetes","c__Leotiomycetes","c__Malasseziomycetes","c__Saccharomycetes","c__Sordariomycetes","Tremellomycetes","Ustilaginomycetes."))]), Fungi1920_vine_AV_L_100GlomClass)

##Plot Bars for fungi in selected taxa## 
plot_bar(Fungi1920_vine_AV_L_100GlomClassSelection, x="comb_num", fill="Class", facet_grid ="Vintage") + geom_col()
    

##RODITIS & SIDERITIS##
##Subset DataSet to cultivars Roditis and Sideritis from different terroir units of Aigialeia's viticultural zone and construct Barplots separately for each plant_part (fruits and leaves)##
# Remove samples that don't belong to the cultivars Roditis and Sideritis #

##subset samples##
Fungi1920_vine_RS <- subset_samples(VINE_ITS_All_2019.2020., !(variety=="Agiorgitiko"))
Fungi1920_vine_RS <- subset_samples(Fungi1920_vine_RS, !(variety=="Vidiano"))
Fungi1920_vine_RS <- prune_taxa(taxa_sums(Fungi1920_vine_RS)>0,Fungi1920_vine_RS)


##RODITIS & SIDERITIS- FRUITS ##
# subset samples-Keep Fruits#
Fungi1920_vine_RS_F <- subset_samples(Fungi1920_vine_RS, !(plant_part=="Leaves"))
Fungi1920_vine_RS_F <- prune_taxa(taxa_sums(Fungi1920_vine_RS_F)>0,Fungi1920_vine_RS_F)

##check sample_data after separating Fruits##
sample_data(Fungi1920_vine_RS_F)

##transform phyloseq object raw counts to relative abundance (100%)##
Fungi1920_vine_RS_F_100 <- transform_sample_counts(Fungi1920_vine_RS_F, function(OTU) 100*OTU/sum(OTU))

##agglomerate in Class level##
Fungi1920_vine_RS_F_100GlomClass <-tax_glom(Fungi1920_vine_RS_F_100, taxrank = "Class")

##save phyloseq object and continue with the statistical analysis scripts##
saveRDS(Fungi1920_vine_RS_F_100GlomClass, file = "Fungi1920_vine_RS_F_100GlomClass.RDS")

##subset Class to the desirable taxa for plotting,in our case Agaricomycetes, Cystobasidiomycetes, Dothideomycetes, Eurotiomycetes, Leotiomycetes, Microbotryomycetes, Saccharomycetes, Sordariomycetes, Tremellomycetes. Average > 1% calculated in Excel##
Fungi1920_vine_RS_F_100GlomClassSelection <- prune_taxa(row.names(tax_table(Fungi1920_vine_RS_F_100GlomClass)[which(tax_table(Fungi1920_vine_RS_F_100GlomClass)[,"Class"]%in%c("c__Agaricomycetes","c__Cystobasidiomycetes","c__Dothideomycetes","c__Eurotiomycetes","c__Leotiomycetes","c__Microbotryomycetes","c__Saccharomycetes","c__Sordariomycetes","Tremellomycetes"))]), Fungi1920_vine_RS_F_100GlomClass)

##Plot Bars for fungi in selected taxa## 
plot_bar(Fungi1920_vine_RS_F_100GlomClassSelection, x="comb_num", fill="Class", facet_grid ="Vintage") + geom_col()

##RODITIS & SIDERITIS- LEAVES ##
# subset samples-Keep Leaves#
Fungi1920_vine_RS_L <- subset_samples(Fungi1920_vine_RS, !(plant_part=="Fruit"))
Fungi1920_vine_RS_L <- prune_taxa(taxa_sums(Fungi1920_vine_RS_L)>0,Fungi1920_vine_RS_L)

##check sample_data after separating Fruits##
sample_data(Fungi1920_vine_RS_L)

##transform phyloseq object raw counts to relative abundance (100%)##
Fungi1920_vine_RS_L_100 <- transform_sample_counts(Fungi1920_vine_RS_L, function(OTU) 100*OTU/sum(OTU))

##agglomerate in Class level##
Fungi1920_vine_RS_L_100GlomClass <-tax_glom(Fungi1920_vine_RS_L_100, taxrank = "Class")

##save phyloseq object and continue with the statistical analysis scripts##
saveRDS(Fungi1920_vine_RS_L_100GlomClass, file = "Fungi1920_vine_RS_L_100GlomClass.RDS")

##subset Class to the desirable taxa for plotting,in our case Agaricomycetes, Cystobasidiomycetes, Dothideomycetes, Eurotiomycetes, Leotiomycetes, Microbotryomycetes, Saccharomycetes, Sordariomycetes, Tremellomycetes. Average > 1% calculated in Excel##
Fungi1920_vine_RS_L_100GlomClassSelection <- prune_taxa(row.names(tax_table(Fungi1920_vine_RS_L_100GlomClass)[which(tax_table(Fungi1920_vine_RS_L_100GlomClass)[,"Class"]%in%c("c__Agaricomycetes","c__Cystobasidiomycetes","c__Dothideomycetes","c__Eurotiomycetes","c__Leotiomycetes","c__Microbotryomycetes","c__Saccharomycetes","c__Sordariomycetes","Tremellomycetes"))]), Fungi1920_vine_RS_L_100GlomClass)

##Plot Bars for fungi in selected taxa## 
plot_bar(Fungi1920_vine_RS_L_100GlomClassSelection, x="comb_num", fill="Class", facet_grid ="Vintage") + geom_col()
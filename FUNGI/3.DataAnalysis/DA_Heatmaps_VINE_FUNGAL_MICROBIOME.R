##D.A Heatmaps##

##D.A Heatmaps must be applied to cultivars Agiorgitiko-Vidiano separately from cultivars Roditis-Sideritis## 

## AGIORGITIKO- VIDIANO ##
##D.A Heatmap Agiorgitiko-Vidiano##
##load phyloseq object for fungi, (VINE_ITS_All_2019.2020..RDS)##
VINE_ITS_All_2019.2020. <- readRDS("../../2.PhyloseqObjectPrep/VINE_ITS_All_2019.2020..RDS")

# Remove samples that don't belong to the cultivars Agiorgitiko and Vidiano #
##subset samples##
Fungi1920_vine_AV <- subset_samples(VINE_ITS_All_2019.2020., !(Cultivar=="Roditis"))
Fungi1920_vine_AV <- subset_samples(Fungi1920_vine_AV, !(Cultivar=="Sideritis"))
Fungi1920_vine_AV <- prune_taxa(taxa_sums(Fungi1920_vine_AV)>0,Fungi1920_vine_AV)


#rename phyloseq object###
psdt <- Fungi1920_vine_AV

##transform phyloseq object raw counts to relative abundance##
psdt_ra <- transform_sample_counts(psdt, function(x) x / sum(x))

##Define the number of OTUs to be used in the analysis, 200 ASVs in our case## 
mynumOTUs <- 200

## Prepare working object and select experimental variable of interest##
ps_htmp <- prune_taxa(names(taxa_sums(psdt_ra)[order(taxa_sums(psdt_ra), decreasing = T)][1:mynumOTUs]),psdt_ra)

## Prepare experimental variables table##
mydesign <- data.frame(sample_data(ps_htmp))

## Prepare the relative abundance table for the top 200 ASVs
dtcounts <- decostand(data.frame(otu_table(ps_htmp))[row.names(mydesign),], method = "range", MARGIN = 2)

## ASVs 
myotus <- colnames(dtcounts)

## Run Kruskal-Wallis, non-parametrix multiple comparison test##
mystatsout <- list()

for(myotu in myotus){
  if(sum(dtcounts[,myotu]) == 0){
    mystatsout[[myotu]]$krusk$statistics$p.chisq <- 1
  } else {
    mykrusk <- kruskal(dtcounts[,myotu], mydesign[rownames(dtcounts),]$Viticultural_Zones, group = T)
    mystatsout[[myotu]][["krusk"]] <- mykrusk
  }
}

## Extract p-values from the test results without p.adj##
mystatsoutkruskpvals <- unlist(lapply(myotus, function(x) mystatsout[[x]]$krusk$statistics$p.chisq))
names(mystatsoutkruskpvals) <- myotus

## Adjust p-values with "False Discovery Rate" method
mystatsoutkruskpvals.adj <- p.adjust(mystatsoutkruskpvals, method = "fdr")

# Check if there are any ASVs that significantly differ with p.adjust
for (i in c(1:length(mystatsoutkruskpvals.adj))) {
  if (mystatsoutkruskpvals.adj[i] < 0.05) {
    print(names(mystatsoutkruskpvals.adj[i]))
  }
}
###without p.adjust
for (i in c(1:length(mystatsoutkruskpvals))) {
  if (mystatsoutkruskpvals[i] < 0.05) {
    print(names(mystatsoutkruskpvals[i]))
  }
}
##p.adjust was preferred to continue if was applicable to the data## 
##Tables preparation of the most abundant micororganisms## 
# prepare the ASV names
mytax <- data.frame(tax_table(ps_htmp), stringsAsFactors = F)

# For ITS - Remove letter from taxonomy
for (i in c(1:nrow(mytax))) {
  for(j in c(1:ncol(mytax))) {
    mytax[i,j] <- gsub("[a-z]__","",mytax[i,j])
  }
}

mytxplot <- data.frame(OTU = row.names(mytax), 
                       txplt = paste(row.names(mytax), " ", mytax$Phylum,  ":", mytax$Genus,  sep = ""))

##Prepare stars, select among (mystatsoutkruskpvals,mystatsoutkruskpvals.adj)##
##without p-adjustment)##
mystarspheat <- array()
for(myotu_sig in mytxplot$OTU){
  mystarspheat[myotu_sig] <- stars.pval(get(paste("mystatsoutkruskpvals", sep = ""))[myotu_sig])
}

## OR with p-adjustment##
mystarspheat <- array()
for(myotu_sig in mytxplot$OTU){
  mystarspheat[myotu_sig] <- stars.pval(get(paste("mystatsoutkruskpvals.adj", sep = ""))[myotu_sig])
}

## Remove samples with missing values (NA)
mystarspheat1 <- mystarspheat[complete.cases(mystarspheat)]
mystarspheat1_tbl <- data.frame(star=mystarspheat1, row.names = names(mystarspheat1))
### Prepare final data and names and stars table
##transpose data##
dtcounts <- t(dtcounts) ##antimetathesi##


counts_nms <- merge(mytxplot,dtcounts, by.x = "OTU", by.y = "row.names", all = T)


counts_nms_stars <- merge(mystarspheat1_tbl,counts_nms, by.x = "row.names",by.y = "OTU", all = T)


counts_nms_stars$nms_star <- paste(counts_nms_stars$txplt,counts_nms_stars$star)


counts_fin <- counts_nms_stars[,4:(ncol(counts_nms_stars)-1)]
row.names(counts_fin) <- counts_nms_stars$nms_star

# Prepeare the mean RA % for each taxon in all samples##
myrelabund <- data.frame("mean RA %" = 100*rowMeans(counts_fin), check.names = F)

# Merge the taxon abundance per sample (counts_fin) with the mean RA
counts_fin_relabund <- merge(counts_fin,myrelabund, by = "row.names", all.x = T)
# Remove not-needed columns
counts_fin_relabund1 <- counts_fin_relabund[,-c(1,ncol(counts_fin_relabund))]

rownames(counts_fin_relabund1) <- counts_fin_relabund$Row.names

##Aggregate table according to Viticultural_Zones##
##Aggregated table
counts_fin_t <- t(counts_fin)

counts_fin_agg <- aggregate(counts_fin_t, by = list(mydesign$Viticultural_Zones), mean)
rownames(counts_fin_agg) <- counts_fin_agg$Group.1


counts_fin_agg_t <- t(counts_fin_agg[,-grep("Group.1",colnames(counts_fin_agg))])


counts_fin_agg_relabund <- merge(counts_fin_agg_t,myrelabund, by = "row.names", all.x = T)

# Prepare the table with only the treatments per column
counts_fin_agg_relabund1 <- counts_fin_agg_relabund[,grep(paste(levels(mydesign$Viticultural_Zones),collapse = "|"),colnames(counts_fin_agg_relabund))]
# Remove not-needed columns
counts_fin_agg_relabund1 <- counts_fin_agg_relabund[,-c(1,ncol(counts_fin_agg_relabund))]

row.names(counts_fin_agg_relabund1) <- paste(counts_fin_agg_relabund$Row.names)

##Plot Aggregated table##
mat_fr_plot_agg <- counts_fin_agg_relabund1[grep("\\*",row.names(counts_fin_agg_relabund1)),]

##you can write the table##
write.table(mat_fr_plot_agg, file = "mat_fr_plot_agg.txt")

##Prepare function variables##
htmp_mat <- mat_fr_plot_agg  ##Aggregated Table##

####HeatmapPlot######
htmp_plot <- function(htmp_mat, mycolors, mycuttreecols, mycuttreerows, mean_RA) { 
  pheatmap(htmp_mat, color = mycolors, cluster_rows = T,
           fontsize = 15,
           cutree_cols = mycuttreecols, cutree_rows = mycuttreerows,
           annotation_row = mean_RA, 
           border_color = "grey30",
           cellwidth = 15, cellheight = 15)
}
## Heatmap colors##
mycolors <- colorRampPalette(c("grey95","lightsteelblue4","darkred","darkorange2"))(n = 119)
#mycuttreecols 
mycutree_col <- cutree(hclust(vegan::vegdist(htmp_mat)), h = 0.8)
mycuttreecols <- max(mycutree_col)
#mycuttreerows
mycutree_row <- cutree(hclust(vegan::vegdist(htmp_mat)), h = 0.8)
mycuttreerows<- max(mycutree_row)
#mean_RA
mean_RA <- myrelabund

cairo_pdf(paste("DA_Fungi1920_vine_AV.pdf", sep = ""), width = 35, height =15)

htmp_plot(htmp_mat = htmp_mat,
          mycolors = mycolors,
          mycuttreecols = mycuttreecols,
          mycuttreerows = mycuttreerows,
          mean_RA = mean_RA)

dev.off()

##____________________________________________________________________##
## RODITIS - SIDERITIS ##
##D.A Heatmap Roditis-Sideritis##
# Remove samples that don't belong to the cultivars Roditis and Sideritis #

##subset samples##
Fungi1920_vine_RS <- subset_samples(VINE_ITS_All_2019.2020., !(variety=="Agiorgitiko"))
Fungi1920_vine_RS <- subset_samples(Fungi1920_vine_RS, !(variety=="Vidiano"))
Fungi1920_vine_RS <- prune_taxa(taxa_sums(Fungi1920_vine_RS)>0,Fungi1920_vine_RS)

#Rename 
psdt <- Fungi1920_vine_RS

##transform phyloseq object raw counts to relative abundance##
psdt_ra <- transform_sample_counts(psdt, function(x) x / sum(x))

##Define the number of OTUs to be used in the analysis, 200 ASVs in our case## 
mynumOTUs <- 200

## Prepare working object and select experimental variable of interest##
ps_htmp <- prune_taxa(names(taxa_sums(psdt_ra)[order(taxa_sums(psdt_ra), decreasing = T)][1:mynumOTUs]),psdt_ra)

## Prepare experimental variables table##
mydesign <- data.frame(sample_data(ps_htmp))

## Prepare the relative abundance table for the top 200 ASVs
dtcounts <- decostand(data.frame(otu_table(ps_htmp))[row.names(mydesign),], method = "range", MARGIN = 2)

## ASVs 
myotus <- colnames(dtcounts)

## Run Kruskal-Wallis, non-parametrix multiple comparison test##
mystatsout <- list()

for(myotu in myotus){
  if(sum(dtcounts[,myotu]) == 0){
    mystatsout[[myotu]]$krusk$statistics$p.chisq <- 1
  } else {
    mykrusk <- kruskal(dtcounts[,myotu], mydesign[rownames(dtcounts),]$Terroir_Unit, group = T)
    mystatsout[[myotu]][["krusk"]] <- mykrusk
  }
}

## Extract p-values from the test results without p.adj##
mystatsoutkruskpvals <- unlist(lapply(myotus, function(x) mystatsout[[x]]$krusk$statistics$p.chisq))
names(mystatsoutkruskpvals) <- myotus

## Adjust p-values with "False Discovery Rate" method
mystatsoutkruskpvals.adj <- p.adjust(mystatsoutkruskpvals, method = "fdr")

# Check if there are any ASVs that significantly differ with p.adjust
for (i in c(1:length(mystatsoutkruskpvals.adj))) {
  if (mystatsoutkruskpvals.adj[i] < 0.05) {
    print(names(mystatsoutkruskpvals.adj[i]))
  }
}
###without p.adjust
for (i in c(1:length(mystatsoutkruskpvals))) {
  if (mystatsoutkruskpvals[i] < 0.05) {
    print(names(mystatsoutkruskpvals[i]))
  }
}
##p.adjust was preferred to continue if was applicable to the data## 
##Tables preparation of the most abundant micororganisms## 
# prepare the ASV names
mytax <- data.frame(tax_table(ps_htmp), stringsAsFactors = F)

# For ITS - Remove letter from taxonomy
for (i in c(1:nrow(mytax))) {
  for(j in c(1:ncol(mytax))) {
    mytax[i,j] <- gsub("[a-z]__","",mytax[i,j])
  }
}

mytxplot <- data.frame(OTU = row.names(mytax), 
                       txplt = paste(row.names(mytax), " ", mytax$Phylum,  ":", mytax$Genus,  sep = ""))

##Prepare stars, select among (mystatsoutkruskpvals,mystatsoutkruskpvals.adj)##
##without p-adjustment)##
mystarspheat <- array()
for(myotu_sig in mytxplot$OTU){
  mystarspheat[myotu_sig] <- stars.pval(get(paste("mystatsoutkruskpvals", sep = ""))[myotu_sig])
}

## OR with p-adjustment##
mystarspheat <- array()
for(myotu_sig in mytxplot$OTU){
  mystarspheat[myotu_sig] <- stars.pval(get(paste("mystatsoutkruskpvals.adj", sep = ""))[myotu_sig])
}

## Remove samples with missing values (NA)
mystarspheat1 <- mystarspheat[complete.cases(mystarspheat)]
mystarspheat1_tbl <- data.frame(star=mystarspheat1, row.names = names(mystarspheat1))
### Prepare final data and names and stars table
##transpose data##
dtcounts <- t(dtcounts) ##antimetathesi##


counts_nms <- merge(mytxplot,dtcounts, by.x = "OTU", by.y = "row.names", all = T)


counts_nms_stars <- merge(mystarspheat1_tbl,counts_nms, by.x = "row.names",by.y = "OTU", all = T)


counts_nms_stars$nms_star <- paste(counts_nms_stars$txplt,counts_nms_stars$star)


counts_fin <- counts_nms_stars[,4:(ncol(counts_nms_stars)-1)]
row.names(counts_fin) <- counts_nms_stars$nms_star

# Prepeare the mean RA % for each taxon in all samples##
myrelabund <- data.frame("mean RA %" = 100*rowMeans(counts_fin), check.names = F)

# Merge the taxon abundance per sample (counts_fin) with the mean RA
counts_fin_relabund <- merge(counts_fin,myrelabund, by = "row.names", all.x = T)
# Remove not-needed columns
counts_fin_relabund1 <- counts_fin_relabund[,-c(1,ncol(counts_fin_relabund))]

rownames(counts_fin_relabund1) <- counts_fin_relabund$Row.names

##Aggregate table according to Terroir_Unit##
##Aggregated table
counts_fin_t <- t(counts_fin)

counts_fin_agg <- aggregate(counts_fin_t, by = list(mydesign$Terroir_Unit), mean)
rownames(counts_fin_agg) <- counts_fin_agg$Group.1


counts_fin_agg_t <- t(counts_fin_agg[,-grep("Group.1",colnames(counts_fin_agg))])


counts_fin_agg_relabund <- merge(counts_fin_agg_t,myrelabund, by = "row.names", all.x = T)

# Prepare the table with only the treatments per column
counts_fin_agg_relabund1 <- counts_fin_agg_relabund[,grep(paste(levels(mydesign$Terroir_Unit),collapse = "|"),colnames(counts_fin_agg_relabund))]
# Remove not-needed columns
counts_fin_agg_relabund1 <- counts_fin_agg_relabund[,-c(1,ncol(counts_fin_agg_relabund))]

row.names(counts_fin_agg_relabund1) <- paste(counts_fin_agg_relabund$Row.names)

##Plot Aggregated table##
mat_fr_plot_agg <- counts_fin_agg_relabund1[grep("\\*",row.names(counts_fin_agg_relabund1)),]

##you can write the table##
write.table(mat_fr_plot_agg, file = "mat_fr_plot_agg.txt")

##Prepare function variables##
htmp_mat <- mat_fr_plot_agg  ##Aggregated Table##

####HeatmapPlot######
htmp_plot <- function(htmp_mat, mycolors, mycuttreecols, mycuttreerows, mean_RA) { 
  pheatmap(htmp_mat, color = mycolors, cluster_rows = T,
           fontsize = 15,
           cutree_cols = mycuttreecols, cutree_rows = mycuttreerows,
           annotation_row = mean_RA, 
           border_color = "grey30",
           cellwidth = 15, cellheight = 15)
}
## Heatmap colors##
mycolors <- colorRampPalette(c("grey95","lightsteelblue4","darkred","darkorange2"))(n = 119)
#mycuttreecols 
mycutree_col <- cutree(hclust(vegan::vegdist(htmp_mat)), h = 0.8)
mycuttreecols <- max(mycutree_col)
#mycuttreerows
mycutree_row <- cutree(hclust(vegan::vegdist(htmp_mat)), h = 0.8)
mycuttreerows<- max(mycutree_row)
#mean_RA
mean_RA <- myrelabund

cairo_pdf(paste("DA_Fungi1920_vine_RS.pdf", sep = ""), width = 35, height =15)

htmp_plot(htmp_mat = htmp_mat,
          mycolors = mycolors,
          mycuttreecols = mycuttreecols,
          mycuttreerows = mycuttreerows,
          mean_RA = mean_RA)

dev.off()

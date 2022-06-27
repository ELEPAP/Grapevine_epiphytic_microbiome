##Quality Control, classification and phyloseq object construction##
##install the necessary packages##
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2")


##load the package##
library(dada2); packageVersion("dada2")

mymultthread <- 56

##input file path and storing in variables##
##change to the directory containing the fastq files after unzipping##
path1 <- "/mnt/4abcb8af-e2eb-47e2-86d6-eed71f6d5304/00_home_guests_do_not_move/fotisbs/Fotis2020/fotisbac1a/lOTUs_moth_basis/" 
path2 <- "/mnt/4abcb8af-e2eb-47e2-86d6-eed71f6d5304/00_home_guests_do_not_move/fotisbs/Fotis2020/fotisbac2a//lOTUs_moth_basis/"
path3 <- "/home/elepap/FOT1_2021/FOT1_2021_Bac/demuxd_raw/bac/"
path4 <- "/home/elepap/FOT2_2021/FOT2_2021_Bac/demuxd_raw/bac/"


##lists files in a path##
# lists files in the path
list.files(c(path1,path2,path3,path4))



# set the variables containing all the forward and the reverse paths to the files of interest with the list.files command#
fnFs <- sort(list.files(c(path1,path2,path3,path4), pattern="_R1.fastq", full.names = TRUE))

fnRs <- sort(list.files(c(path1,path2,path3,path4), pattern="_R2.fastq", full.names = TRUE))


# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq#
sample.names <- gsub("_lib_.+","",basename(fnFs))


# Plot the per-base qualities of all the samples in a single file
pdf("Initial Quality Files.pdf",onefile = T)
for (i in c(1:length(fnFs))) {
  print(i)
  plot1 <- plotQualityProfile(c(fnFs[i], fnRs[i]))
  print(plot1)
}
dev.off()


#### Sequence quality filtering and control, error modelling, and dereplication ----
# Set the file paths where the quality controlled sequences are going to be saved#
filtFs <- file.path("filtered2", paste(sample.names, "_F_filt.fastq.gz", sep = ""))
filtRs <- file.path("filtered2", paste(sample.names, "_R_filt.fastq.gz", sep = ""))


# Filter the sequences and save then in the folders provided above and get their statistics in a table#
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     trimLeft = 11, maxN=0,  
                     maxEE=c(2,2), truncQ=2, 
                     rm.phix=TRUE, compress=TRUE, 
                     multithread=TRUE, matchIDs=TRUE)

# View(out)
head(out)


# Learn error rates using a machine learning algorithm #
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)


# Visualize the estimated error rates in a single pdf
pdf("Estimated Error Rates.pdf",onefile = T)
plotErrF <- plotErrors(errF, nominalQ=TRUE)
print(plotErrF)
plotErrR <- plotErrors(errR, nominalQ=TRUE)
print(plotErrR)
dev.off()


# Dereplication of each one of the read pairs to unique sequences (collapsing of the identical sequences for each pair per sample)#
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)


# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


# sample composition inference after read correction #
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


# merge read pairs retaining the per sequence sample information
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


#### construct the sequence table, remove the chimeras, and create a summary ----
# construct sequence table (ASV table)
seqtab <- makeSequenceTable(mergers)

# View(seqtab)

#Number of samples
nrow(seqtab)
#Number of sequence Variants- ASVs
ncol(seqtab)

#Distribution of sequence lengths #
table(nchar(getSequences(seqtab)))


# Chimera removal (consensus, pooled or per sample)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

##dimensions of seqtab.nochim##
dim(seqtab.nochim)

##record the portion of good sequences out of the total prior the chimera removal##
sum(seqtab.nochim)/sum(seqtab)

# Track reads through the pipeline
getN <- function(x) {sum(getUniques(x))}

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#View(track)
head(track)

# Save the read quality control data (if you want to download the table at your computer, you need to go to the "More" option and select export)
write.table(track, file = "readQCBacteria.txt", col.names = NA, sep = "\t", quote = FALSE)


##taxonomically classify the sequences##
taxa <-assignTaxonomy(seqs = seqtab.nochim, minBoot = 80, refFasta = "silva_nr_v138_train_set.fa.gz", multithread=TRUE, tryRC = TRUE)


write.table(taxa, file = "taxa.txt", col.names = NA, sep = "\t", quote = F)

##proceed with analysis, load the phyloseq package##
library(phyloseq); packageVersion("phyloseq")

##load the ggplot2 package, for plotting functions##
library(ggplot2); packageVersion("ggplot2")

# the following globally sets the theme of the plots
theme_set(theme_bw())

# load the experimental design data##
samdf <- read.table(file = "DESIGN.txt", header = TRUE, row.names = 1, sep = "\t")

#View(samdf)

##construct the phyloseq object## 
VINE_16S_19.20 <- phyloseq(otu_table(seqtab.nochim, 
                                             taxa_are_rows=FALSE),
                                   sample_data(samdf),
                                   tax_table(taxa))

## Replace the taxon names (the sequences of each ASV with something easier to read and save the sequence info in the object)##
##load the Biostrings and stringr packages##
 
# First install Biostrings which is required for saving the sequences as fasta
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Biostrings")
library("Biostrings")

# The following retrieves the sequences from the taxa table in a convenient for storing format
sequences <- Biostrings::DNAStringSet(taxa_names(VINE_16S_19.20))

names(sequences) <- taxa_names(VINE_16S_19.20)

# The actual addition of the sequences is performed with the following command
ps <- merge_phyloseq(VINE_16S_19.20, sequences)

VINE_16S_All_2019.2020. <- ps

# the following command (uses the stringr::str_pad for padding leading zeros) performs the actual renaming
library(stringr)
taxa_names(VINE_16S_All_2019.2020.) <- paste("ASV",str_pad(1:length(taxa_names(VINE_16S_All_2019.2020.)),5, pad = "0"),sep = "")

rank_names(VINE_16S_All_2019.2020.)

VINE_16S_All_2019.2020.

##Create table,and check number of features for each phyla##
table(tax_table(VINE_16S_All_2019.2020.)[, "Kingdom"], exclude = NULL)
table(tax_table(VINE_16S_All_2019.2020.)[, "Phylum"], exclude = NULL) 
table(tax_table(VINE_16S_All_2019.2020.)[, "Class"], exclude = NULL) 
table(tax_table(VINE_16S_All_2019.2020.)[, "Order"], exclude = NULL)
table(tax_table(VINE_16S_All_2019.2020.)[, "Family"], exclude = NULL)
table(tax_table(VINE_16S_All_2019.2020.)[, "Genus"], exclude = NULL)

#Further, features with ambiguous annotation, low abundance or non target taxa removed from phyloseq object##
##For our analysis removed Eukaryota, Archaea, NA##
ps_bac_wd <- subset_taxa(VINE_16S_All_2019.2020., !is.na(Kingdom) & !Kingdom %in% (c("", "uncharacterized", "Eukaryota","Archaea")))

##also remove mitochondria and chloroplasts##

VINE_16S_All_2019.2020. <- subset_taxa(VINE_16S_All_2019.2020., !Order %in% c("Chloroplast"))

VINE_16S_All_2019.2020. <- subset_taxa(VINE_16S_All_2019.2020., !Family %in% c("Mitochondria"))

VINE_16S_All_2019.2020. <- prune_taxa(taxa_sums(VINE_16S_All_2019.2020.)>0, VINE_16S_All_2019.2020.)

##replace the NA names with the classified leftmost (higher level annotated) taxa##
VINE_16S_All_2019.2020. <- VINE_16S_All_2019.2020.

for(i in 1:nrow(tax_table(VINE_16S_All_2019.2020.))){
  for(j in 2:ncol(tax_table(VINE_16S_All_2019.2020.))){
    if(is.na(tax_table(VINE_16S_All_2019.2020.)[i,j])){
      tax_table(VINE_16S_All_2019.2020.)[i,j] <- tax_table(VINE_16S_All_2019.2020.)[i,j-1]
    }
  }
}


##save phyloseq object and continue with the statistical analysis scripts##
saveRDS(VINE_16S_All_2019.2020., file = "VINE_16S_All_2019.2020..RDS")
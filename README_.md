Different factors are operative in shaping grapevine microbiome across different geographical scales: biogeography, cultivar or vintage?
By Papadopoulou E. 1, Bekris F. 1, Vasileiadis S. 1, Papadopoulou K.K. 1, Karpouzas D.G. 1*
(* corr. author)
1 University of Thessaly, Department of Biochemistry and Biotechnology, Laboratory of Plant and Environmental Biotechnology, Viopolis – 41500 Larissa, Greece

The provided material includes the code used in the statistical analysis of the study.
For obtaining the code the users need to open a terminal and having the GitHub tools, git-clone or download the repository, and enter the base folder. E.g:

$ git clone https://github.com/ELEPAP/Grapevine_epiphytic_microbiome-
In the case of the computational methods, with the "Grapevine_epiphytic_microbiome-" folder as working directory, and assuming that the necessary software and R packages are installed, the used code can be executed as described in this Readme.md file. The necessary datasets for performing all sequencing based analysis can be downloaded implementing the code provided in the corresponding repository folders as explained below.

Description of the order of executed scripts.
Steps 1-3 concern the data retrieval from NCBI and preprocessing, while steps 4-6 concern the actual data analysis for total fungi and bacteria.

First, it is necessary to download the sequencing data. To do so, you need to enter the "0.DownloadData" subfolder of e.g. the "Fungi" and execute the "fetch_data.sh" bash script (this assumes that you are located at the working directory "Grapevine_epiphytic_microbiome-"). The script is based on the SRR accession numbers found in the 0.DownloadData folder. Once the download is done, you need to combine all forward reads to a single file and all reverse reads to another file as well.
for i in {1..3}
do
	cd Fungi/0.DownloadData/batch${i}
	sh -x fetch_data.sh
	cat *_1.fastq | gzip > forward.fastq.gz
	cat *_2.fastq | gzip > reverse.fastq.gz
	cd ../../../
	cd Bacteria/0.DownloadData/batch${i}
	sh -x fetch_data.sh
	cat *_1.fastq | gzip > forward.fastq.gz
	cat *_2.fastq | gzip > reverse.fastq.gz
	cd ../../../
done
Then you need to demultiplex the data according to our own demultiplexing method using our in-house script. This requires Flexbar v3.0.3 to be installed as described in the manuscript. A detailed description of our in-house multiplexing approach is provided in our [previous work] (https://github.com/SotiriosVasileiadis/mconsort_tbz_degr#16s). You need to enter the folder Fungi/1.Demultiplex and run the following commands (change the MY_PROCS variable to whatever number of logical processors you have available and want to devote). the following commands are going to save the demultiplexed files in the Fungi(or Bacteria)/1.Demultiplex/demux_out folder.
MY_WORKING_DIR_BASE=`pwd`
for i in {1..3}
do
  cd Fungi/1.Demultiplex
  MY_PROCS=56
  bash DemuxOwnBCsys_absPATH.sh demux_out${i} ${MY_WORKING_DIR_BASE}/Fungi/0.DownloadData/batch${i}/forward.fastq.gz ${MY_WORKING_DIR_BASE}/Fungi/0.DownloadData/batch${i}/reverse.fastq.gz fun${i}_map_file.txt ${MY_PROCS}
  cd demux_out${i}/analysis_ready
  gunzip *.gz # unzips files skipped by the Demux script
  cd ../../../../
  cd Bacteria/1.Demultiplex
  MY_PROCS=56
  bash DemuxOwnBCsys_absPATH.sh demux_ou${i} ${MY_WORKING_DIR_BASE}/Fungi/0.DownloadData/batch${i}/forward.fastq.gz ${MY_WORKING_DIR_BASE}/Fungi/0.DownloadData/batch${i}/reverse.fastq.gz bac${i}_map_file.txt ${MY_PROCS}
  cd demux_out${i}/analysis_ready
  gunzip *.gz # unzips files skipped by the Demux script
  cd ../../../../
done

cd Fungi/1.Demultiplex
mkdir -p demux_out/analysis_ready
cp demux_out[0-9]/analysis_ready/*.fastq demux_out/analysis_ready/
cd ../../

cd Bacteria/1.Demultiplex
mkdir -p demux_out/analysis_ready
cp demux_out[0-9]/analysis_ready/*.fastq demux_out/analysis_ready/
cd ../../
Following, the "PhyloseqObjectPrep_Fun_Vine1920.R" script of the Fungi and the "PhyloseqObjectPrep_Bac_Vine1920.R" of the Bacteria/2.PhyloseqObjectPerp folder is run in order to prepare the final phyloseq object to be used in the data analysis described below. Before runnin gthe script make sure that the necessary reference databases are found in the same folder.
cd Fungi/2.PhyloseqObjectPrep
# fetch the databases
wget https://files.plutof.ut.ee/public/orig/1D/B9/1DB95C8AC0A80108BECAF1162D761A8D379AF43E2A4295A3EF353DD1632B645B.gz
# run the R script
Rscript PhyloseqObjectPrep_Fun_Vine1920.R
cd ../../
cd Bacteria/2.PhyloseqObjectPrep
# fetch the databases
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz
tar vxf *.gz
# run the R script
Rscript PhyloseqObjectPrep_Bac_Vine1920.R
cd ../../

4a) Run the overall PERMANOVA tests.
cd Fungi/3.DataAnalysis/PERMANOVA
Rscript PERMANOVA.R
cd ../../../
cd Bacteria/3.DataAnalysis/PERMANOVA
Rscript PERMANOVA.R
cd ../../../

4b) Run the NMDS and PERMANOVA tests.
cd Fungi/3.DataAnalysis/NMDS_PERMANOVA_PAIRWISE_FUNGIVINE1920_AGIORGITIKO_VIDIANO
Rscript NMDS_PERMANOVA_PAIRWISE_FUNGIVINE1920_AGIORGITIKO_VIDIANO.R
cd Fungi/3.DataAnalysis/NMDS_PERMANOVA_PAIRWISE_FUNGIVINE1920_RODITIS_SIDERITIS
Rscript NMDS_PERMANOVA_PAIRWISE_FUNGIVINE1920_RODITIS_SIDERITIS.R
cd ../../../
cd Bacteria/3.DataAnalysis/NMDS_PERMANOVA_PAIRWISE_BACTERIAVINE1920_AGIORGITIKO_VIDIANO_FL
Rscript NMDS_PERMANOVA_PAIRWISE_BACTERIAVINE1920_AGIORGITIKO_VIDIANO_FL.R
cd Bacteria/3.DataAnalysis/NMDS_PERMANOVA_PAIRWISE_BACTERIAVINE1920_RODITIS_SIDERITIS_FL
Rscript NMDS_PERMANOVA_PAIRWISE_BACTERIAVINE1920_RODITIS_SIDERITIS_FL.R
cd ../../../

4c) Prepare the barplots.
cd Fungi/3.DataAnalysis/BARPLOTS_FUNGI_VINE1920
Rscript BARPLOTS_FUNGI_VINE1920.R
cd ../../../
cd Bacteria/3.DataAnalysis/BARPLOTS_BACTERIA-1920_VINE_MICROBIOME
Rscript BARPLOTS_BACTERIA-1920_VINE_MICROBIOME.R
cd ../../../

4d) Prepare the differential abundance (DA) heatmaps.
cd Fungi/3.DataAnalysis/DA_Heatmaps_VINE_FUNGAL_MICROBIOME
Rscript DA_Heatmaps_VINE_FUNGAL_MICROBIOME.R
cd ../../../
cd Bacteria/3.DataAnalysis/DA_Heatmaps_VINE_BACTERIA_MICROBIOME
Rscript DA_Heatmaps_VINE_BACTERIA_MICROBIOME.R
cd ../../../

4e) Run the core microbiome analysis.
cd Fungi/3.DataAnalysis/CORE_MICROBIOME_RODITIS_fungi
Rscript CORE_MICROBIOME_RODITIS_fungi.R
cd Fungi/3.DataAnalysis/CORE_MICROBIOME_SIDERITIS_fungi
Rscript CORE_MICROBIOME_SIDERITIS_fungi.R
cd ../../../
cd Fungi/3.DataAnalysis/CORE_MICROBIOME_RODITIS_bacteria
Rscript CORE_MICROBIOME_RODITIS_bacteria.R
cd Fungi/3.DataAnalysis/CORE_MICROBIOME_SIDERITIS_bacteria
Rscript CORE_MICROBIOME_SIDERITIS_bacteria.R
cd ../../../

Code Usage disclaimer
The following is the disclaimer that applies to all scripts, functions, one-liners, etc. This disclaimer supersedes any disclaimer included in any script, function, one-liner, etc.

You running this script/function means you will not blame the author(s) if this breaks your stuff. This script/function is provided AS IS without warranty of any kind. Author(s) disclaim all implied warranties including, without limitation, any implied warranties of merchantability or of fitness for a particular purpose. The entire risk arising out of the use or performance of the sample scripts and documentation remains with you. In no event shall author(s) be held liable for any damages whatsoever (including, without limitation, damages for loss of business profits, business interruption, loss of business information, or other pecuniary loss) arising out of the use of or inability to use the script or documentation. Neither this script/function, nor any part of it other than those parts that are explicitly copied from others, may be republished without author(s) express written permission. Author(s) retain the right to alter this disclaimer at any time. This disclaimer was copied from a version of the disclaimer published by other authors in https://ucunleashed.com/code-disclaimer and may be amended as needed in the future.
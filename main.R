##################################################################################
#                                                                                #
# PiCnIc/TRONCO Examples -- CRC Case Study COADREAD                              #
#                                                                                #
##################################################################################
# Copyright (c) 2015, Giulio Caravagna, Luca De Sano, Daniele Ramazzotti         #
# email: tronco@disco.unimib.it                                                  #
# All rights reserved. This program and the accompanying materials               #
# are made available under the terms of the GNU GPL v3.0                         #
# which accompanies this distribution                                            #
#                                                                                #
##################################################################################

# You might install TRONCO's version from our Github as the lates version, which is
# development, or master (stable)
if(!require(devtools)) install.packages('devtools')
install_github("BIMIB-DISCo/TRONCO", ref = 'development')
library(TRONCO)

# Set SINK to FALSE to avoid creating a log file
SINK = TRUE
if(SINK) sink(paste0(getwd(), "/PiCnIc-COADREAD-logfile.txt"), append=FALSE, split=TRUE)

#setwd('/Volumes/DATA/Work/Software/Github/TRONCO')
#library(devtools)
#document()
#setwd('/Volumes/DATA/Work/Software/Github/PiCnIc-COADREAD')

#Working directory
workdir = "TCGA-data/"
dir.create(workdir)

# Data files
datafile = 'TCGA-COADREAD-TRONCO.zip'
download.file('https://github.com/BIMIB-DISCo/datasets/raw/master/TCGA-COADREAD-TRONCO.zip',
	destfile=datafile,
	method='curl',
	extra='-L')
unzip(datafile, exdir = workdir)

# Name input files
clinical.file = paste0(workdir, "Clinical/crc_clinical_sheet.txt")
MAF.file = paste0(workdir, "Mutations/TCGA_CRC_Suppl_Table2_Mutations_20120719.csv")
GISTIC.file = paste0(workdir, "CNV/crc_gistic.txt")
clusters.file = paste0(workdir, "Clusters/TCGA-clusters.csv") 
mutex.msi.file = paste0(workdir, "Mutex/msi_results.txt") 
mutex.mss.file = paste0(workdir, "Mutex/mss_results.txt") 

###############################################################################
# PiCnIc steps are separately implemented in the files in folder scripts      #
#                                                                             #
# 1 - Cohort subtyping                                                        #
# 2 - Drivers selection                                                       #
# 3 - Detection of groups of exclusive alterations                            #
# 4 - Models inference, confidence assessmnent                                #
###############################################################################

###############################################################################
## Steps: preamble                                                            #
## Prepare COADREAD data, this involved minor manual curation of              #
## project data, which we did offline.                                        #
###############################################################################
source('scripts/TCGA-import.R', echo = T)         

###############################################################################
## Steps: 1/2                                                                 #
## Cohort subtyping (MSI-HIGH/MSS) and drivers selection are done together    #
## as subtype status is provided by TCGA via clinical annotations             #
## and the list of 33 driver genes is prepared by TCGA (manual curation and   #
## MutSigCV execution)                                                        #
###############################################################################
source('scripts/1-2.subtyping-drivers.R', echo = T)      

###############################################################################
## Step: 3                                                                    #
## Detection of groups of mutual exclusivity is done here. This is showing    #
## how to exploit the interface between TRONCO and the MUTEX tool, as well    #
## as how to include groups provided elsewhere (via prior knowledge and       #
## via the MEMO tool -- run by TCGA)                                          #
###############################################################################
source('scripts/3.groups-exclusivity.R', echo = T)

###############################################################################
## Step: 4                                                                    #
## Finalize data processing, build CAPRI's hypotheses from exclusivity groups #
## and run the algorithm                                                      #
###############################################################################
source('scripts/4.reconstruction.R', echo = T)

###############################################################################
## Step: 5                                                                    #
## Assess the statistical confidence of each model. We use various measures   #
###############################################################################
source('scripts/5.statistics.R', echo = T)

tronco.plot(MSI.models, 
	 pathways = pathway.list, 
	 fontsize = 15,
	 edge.cex = 1.5,
 	 legend.cex = .7,
	 scale.nodes = .6,
	 confidence = c('tp', 'pr', 'hg', 'npb'), # Display p-values 
	 pathways.color = pathways.color,
	 label.edge.size = 9,
	disconnected = F, 
	height.logic = .3,
	file = paste0(workdir, '/msi.pdf'))

tronco.plot(MSS.models, 
	 pathways = pathway.list, 
	 fontsize = 15,
	 edge.cex = 1.5,
 	 legend.cex = .7,
	 scale.nodes = .6,
	 confidence = c('tp', 'pr', 'hg', 'npb'), # Display p-values 
	 pathways.color = pathways.color,
	 label.edge.size = 9,
	disconnected = F, 
	height.logic = .3,
	file = paste0(workdir, '/mss.pdf'))

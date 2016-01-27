##################################################################################
#                                                                                #
# TRONCO Examples -- CRC Case Study COADREAD                                     #
#                                                                                #
##################################################################################
# Copyright (c) 2015, Giulio Caravagna, Luca De Sano, Daniele Ramazzotti         #
# email: tronco@disco.unimib.it                                                  #
# All rights reserved. This program and the accompanying materials               #
# are made available under the terms of the GNU GPL v3.0                         #
# which accompanies this distribution                                            #
#                                                                                #
##################################################################################

clinical.data = TCGA.map.clinical.data(
  file = clinical.file, 
  column.samples = 'patient', 
  column.map = 'tumor_stage')
head(clinical.data)

# Load MAF - use is.TCGA to match samples to patients
MAF = import.MAF(
	file = MAF.file, 
	is.TCGA = TRUE, 
	sep = ';')

# Add stage annotation - use match.TCGA.patients to match long/short barcodes
MAF = annotate.stages(MAF, clinical.data, match.TCGA.patients = TRUE)
show(MAF)

# Check for duplicated samples - we find them
TCGA.multiple.samples(MAF)

# Remove duplicated samples according to TCGA criteria, shorten barcodes and add stages
MAF = TCGA.remove.multiple.samples(MAF)
MAF = TCGA.shorten.barcodes(MAF)
MAF = annotate.stages(MAF, clinical.data)

# Load as a plain table
GISTIC = read.table(
  GISTIC.file, 
  check.names = FALSE,
  stringsAsFactors = FALSE, 
  header = TRUE)
  
# Have a look at this table to remove useless information
head(GISTIC[, 1:5])
GISTIC$Entrez_Gene_Id = NULL
rownames(GISTIC) = GISTIC$Hugo_Symbol
GISTIC$Hugo_Symbol = NULL

# Import all GISTIC data 
GISTIC = import.GISTIC( t(GISTIC) )
show(GISTIC)

# We want to use only high-confidence scores in GISTIC, renamed as Amplification/Deletion
GISTIC = delete.type(GISTIC, 'Heterozygous Loss') # low-level deletions
GISTIC = delete.type(GISTIC, 'Low-level Gain')    # low-level amplifications
GISTIC = rename.type(GISTIC, 'Homozygous Loss', 'Deletion')    
GISTIC = rename.type(GISTIC, 'High-level Gain', 'Amplification')    
GISTIC = annotate.stages(GISTIC, clinical.data)
show(GISTIC)

c("PIK3CA", "FAM123B") %in% as.genes(GISTIC) # Shall be FALSE

# we set intersect.genomes to FALSE to take the union of altered genes 
MAF.GISTIC = intersect.datasets(GISTIC, MAF, intersect.genomes = FALSE)

# We remove events which have no observations in the dataset, and annotate stages
MAF.GISTIC = trim(MAF.GISTIC)
MAF.GISTIC = annotate.stages(MAF.GISTIC, clinical.data)
show(MAF.GISTIC)

# Export these datasets as Rdata
save(MAF.GISTIC, file = paste0(workdir, 'MAF.GISTIC.Rdata'))
save(GISTIC, file = paste0(workdir, 'GISTIC.Rdata'))

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

#################################################################################
# Driver events.                                                                #
#                                                                               #
# 33 genes mapped to 5 pathways by TCGA, we declare them here                   #
# as we will load from the MAF and GISTIC files only the mutations annotated to #
# these genes, as this is increasing the speed of running all the script.       #
#                                                                               #
# This list of genes is compiled by manual curation and by running the MutSiGCV #
# tool -- see the TCGA paper for details                                        #
#################################################################################
Wnt = c("APC", "CTNNB1", "DKK1", "DKK2", "DKK3", "DKK4", "LRP5", "FZD10", "FAM123B", 
        "AXIN2", "TCF7L2", "FBXW7", "ARID1A", "SOX9")
RAS = c("ERBB2", "ERBB3", "NRAS", "KRAS", "BRAF")
PI3K = c("IGF2", "IRS2", "PIK3CA", "PIK3R1", "PTEN")
TGFb = c("TGFBR1", "TGFBR2", "ACVR1B", "ACVR2A", "SMAD2", "SMAD3", "SMAD4")
P53 = c("TP53", "ATM")

# Some variable which will be processed by TRONCO plotting functions
pathway.genes = c(Wnt, RAS, PI3K, TGFb, P53)
pathway.names = c('Wnt', 'RAS', 'PI3K', 'TGFb', 'P53')
pathway.list = list(Wnt = Wnt, RAS = RAS, PI3K = PI3K, TGFb = TGFb, P53 = P53)

#################################################################################
# COADREAD mutations.                                                           #
#                                                                               #
# We load all mutations annotated by the consortium in the MAF file that they   #
# provide. TRONCO provides a function to process such data format and transform #
# the listed mutations in a TRONCO object. In doing so, we also use the         #
# following options:                                                            #
#   - option is.TCGA, which checks if every patient is associated to a unique   #
#     patient -- if not, we resolve that with other TRONCO functions. Notice    #
#     that if this is the case the function will raise a warning, which is the  #
#     case for COADREAD samples.                                                #
#   - filter.fun, which is a custom function that one can define to filter out  #
#     some MAF entries that he does not want to process. The function is then   #
#     applied on each row of the MAF, via apply(), so it is supposed to return  #
#     a TRUE/FALSE flag for every . In this case, as we know that the MAF has   #
#     many more genes annotated, we just check if the annotated gene is one of  #
#     the list that we declared. An alternative option would have been to load  #
#     the full MAF and then use events.selection() function with option         #
#     filter.in.names=pathway.genes.                                            #
#################################################################################
MAF = import.MAF(
	file = MAF.file, 
	is.TCGA = TRUE, 
	sep = ';',
	filter.fun = function(x){ return(x['Hugo_Symbol'] %in% pathway.genes) } # filter
	)

# If one wants to check which information are stored in a TRONCO object a set of
# functions are available. Here every alteration inside TRONCO is called "event"
# and is constituted of a name and a type. These could be for example
#
# TP53 (name) "Mutation" (type)
# PTEN (name) "Deletion" (type)
# 3q (name) "Amplification" (type)
# 17p13.1 (name) "Deletion" (type)
#
# or any other type of event that we want to include in TRONCO. Clearly, these
# could be also custom; in this case all the information available in the MAF
# object are of the first type.
#
# For historical reasons we always refer to genes (even when an event refers
# to a cytoband), which means that we can "query" the MAF object as follows.
as.genes(MAF)                     # List of genes
as.genes(MAF, types = 'Mutation') # List of genes, redundant
ngenes(MAF)                       # Number of

# Events can be similary queried
as.events(MAF)                            # List of events
as.events(MAF, genes = c('KRAS', 'TP53')) # List of events, for two specific genes
as.events(MAF,keysToNames = T)            # Again list of events, but compacted
nevents(MAF)                              # Number of

# Types of alterations 
as.types(MAF)                             # Types of events
as.types(MAF, genes = c('KRAS', 'TP53'))  # Types of events, for two specific genes
ntypes(MAF)                               # Number of

# Samples and their alterations
as.samples(MAF)                                                 # Samples list
which.samples(MAF, gene = 'KRAS', type = 'Mutation')            # KRAS mutations
which.samples(MAF, gene = 'KRAS', type = 'Mutation', neg = T)   # KRAS wildtype
nsamples(MAF)                                                   # Number of

# We can annotated clinical stage data for a TRONCO object. We need the TCGA 
# map patient -> stage
clinical.data = TCGA.map.clinical.data(
  file = clinical.file, 
  column.samples = 'patient', 
  column.map = 'tumor_stage')
head(clinical.data)

# Add stage annotation. We use match.TCGA.patients to match long/short barcodes
# This annotation will assign to each sample its stage, and NA to those which 
# do not have that. Warnings will be raised for such samples.
MAF = annotate.stages(MAF, clinical.data, match.TCGA.patients = TRUE)

# This function gives a textual representation of a TRONCO object
view(MAF)

# Check for duplicated samples - we find them in this cohort
TCGA.multiple.samples(MAF)

# Remove duplicated samples according to TCGA criteria, shorten barcodes and add stages
MAF = TCGA.remove.multiple.samples(MAF)
MAF = TCGA.shorten.barcodes(MAF)
MAF = annotate.stages(MAF, clinical.data)

# We annotate the MAF object with a simple textual description which will become
# the figures title
MAF = annotate.description(x = MAF, label = "COADREAD MAF data for driver genes")

# You can see the MAF data loaded as follows. We here exploit the oncoprint visualization
# function of TRONCO which plots the classical grid of mutations/CNAs or other events in
# general. Here we use it in its simples form, with many defaults -- later we will
# show how to produce fancy plots
if(DOPLOTS) oncoprint(MAF)

# install.packages("devtools")
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("AnnotationDbi", "biomaRt", "Biostrings", "GenomicFeatures", "GenomicRanges", "Rsamtools"))
# install.packages(c("rmarkdown", "knitr"))
# devtools::install_github("griffithlab/GenVisR")
# library(GenVisR)
# MAF.dataframe = MAF = import.MAF(file = MAF.file, is.TCGA = TRUE, sep = ';', to.TRONCO = FALSE,
#                                  filter.fun = function(x){ return(x['Hugo_Symbol'] %in% pathway.genes) }) 
# 
# dev.new(noRStudioGD = TRUE)
# 
# MAF.dataframe = MAF.dataframe[which(MAF.dataframe$Variant_Classification != 'De_novo_Start_OutOfFrame'),] # Not supported mutation type
# 
# library(RColorBrewer)
# mut.colors = brewer.pal('Dark2', n = 7)
# names(mut.colors) = unique(MAF.dataframe$Variant_Classification)
# 
# waterfall(MAF.dataframe, mainGrid = T, mainDropMut = T, mainPalette = mut.colors)
# 
#   
#   
# library(reshape2)
# 
# names(MAF.dataframe)


#################################################################################
# COADREAD Copy Number Alterations - CNAs                                       #
#                                                                               #
# We load all CNAs annotated by the consortium in the GISTIC file that they     #
# provide. As for mutations, TRONCO provides a function to process such data    #
# format and transform it in a TRONCO object. In doing so, we also use one      #
# options:                                                                      #
#   - filter.genes = pathway.genes which will ensure that we load only data     #
#################################################################################
GISTIC = import.GISTIC(x = GISTIC.file, filter.genes = pathway.genes)
GISTIC = annotate.stages(GISTIC, clinical.data)
GISTIC = annotate.description(x = GISTIC, label = "COADREAD CNA data for driver genes")

# The imported CNAs contain 4 types of alterations, which can be seen as follows
as.types(GISTIC)

# or equivalently with show/oncoprint functions. The oncoprint is particulalry interesting
# as it shows that GISTIC data is much full of heterozygous losses and low-level gains
# which we do not want to process 
view(GISTIC)
if(DOPLOTS) oncoprint(GISTIC)

# We want to use only high-confidence scores in GISTIC, renamed as Amplification/Deletion
# and thus we use TRONCO's editing functions which allow to modify labels of alterations
# in a certain object, or delete certain alteration types
GISTIC = delete.type(GISTIC, 'Heterozygous Loss') # low-level deletions
GISTIC = delete.type(GISTIC, 'Low-level Gain')    # low-level amplifications
GISTIC = rename.type(GISTIC, 'Homozygous Loss', 'Deletion')    
GISTIC = rename.type(GISTIC, 'High-level Gain', 'Amplification')    
GISTIC = annotate.stages(GISTIC, clinical.data)
view(GISTIC)

# Shall be FALSE, as we know that these genes do not harbour CNAs
c("PIK3CA", "FAM123B") %in% as.genes(GISTIC) 

# Notice that GISTIC contains much more and different samples from MAF, and that not all
# samples in the MAF are also in the GISTIC object
nsamples(GISTIC)
nsamples(GISTIC) > nsamples(MAF)
as.samples(MAF) %in% as.samples(GISTIC) 

# So we have two TRONCO objects, one for MAF, one for CNAs. We want to intersect the two datasets
# to have only samples for which we have both CNA and mutation data available. We use a TRONCO
# function which intersects the samples set, but we ask that not to be strict about genomes via
# intersect.genomes = FALSE. We do that as some genes such as c("PIK3CA", "FAM123B") do not appear
# in both datasets.
MAF.GISTIC = intersect.datasets(GISTIC, MAF, intersect.genomes = FALSE)
view(MAF.GISTIC)

# We remove events which have no observations in the dataset, and annotate stages
MAF.GISTIC = trim(MAF.GISTIC)
MAF.GISTIC = annotate.stages(MAF.GISTIC, clinical.data)
MAF.GISTIC = annotate.description(x = MAF.GISTIC, label = "COADREAD MAF/CNA data for driver genes")
view(MAF.GISTIC)

# View and export these datasets as Rdata
if(DOPLOTS) oncoprint(MAF.GISTIC)
save(MAF.GISTIC, file = paste0(workdir, 'MAF.GISTIC.Rdata'))
save(GISTIC, file = paste0(workdir, 'GISTIC.Rdata'))
save(MAF, file = paste0(workdir, 'MAF.Rdata'))

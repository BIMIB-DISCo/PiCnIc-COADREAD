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

# Prepare folders
dir.create('./MSS')
dir.create('./MSI')
sub.dir = c('mutex', 'Rdata-lifted', 'Rdata-models')
sapply(paste0('./MSS/', sub.dir), dir.create)
sapply(paste0('./MSI/', sub.dir), dir.create)

# Load MAF.GISTIC file,  set  some fancy  colors  to get cute visualization
load(paste0(workdir, '/MAF.GISTIC.Rdata'))
MAF.GISTIC = change.color(MAF.GISTIC, 'Mutation', 'darkolivegreen3')
MAF.GISTIC = change.color(MAF.GISTIC, 'Amplification', 'coral')
MAF.GISTIC = change.color(MAF.GISTIC, 'Deletion', 'cornflowerblue')

# Load table data
file = read.delim(clusters.file, sep = ";")
head(file)

# Select  just certain annotations, remove blank lines
tab = file[, c("patient",  "MSI_status", "sequenced")]
tab = tab[1:276, ] 
rownames(tab) = tab$patient

# Filter out non-sequenced samples, and order them (for console visualization)
tab = tab[tab$sequenced == 1, ]
tab = tab[order(tab$MSI_status), ]
print(tab)

# Define the maps to split samples
map.MSS = tab[tab$MSI_status ==  "MSS", , drop = FALSE]
map.MSI.H = tab[tab$MSI_status ==  "MSI-H", , drop = FALSE]

# These are the samples that we actually use
MSS.samples = rownames(map.MSS)
MSI.H.samples = rownames(map.MSI.H)

# Split is done by using samples.selection with appropriate vectors as input
MSS = trim(samples.selection(MAF.GISTIC, MSS.samples))
MSI.H = trim(samples.selection(MAF.GISTIC, MSI.H.samples))
show(MSS)
show(MSI.H)

alteration.color = 'dimgray'
pathways.color = c('firebrick1', 'darkblue', 'darkgreen', 'darkmagenta', 'darkorange')

# Driver events - 33 genes mapped to 5 pathways by TCGA
Wnt = c("APC", "CTNNB1", "DKK1", "DKK2", "DKK3", "DKK4", "LRP5", "FZD10", "FAM123B", "AXIN2", "TCF7L2", "FBXW7", "ARID1A", "SOX9")
RAS = c("ERBB2", "ERBB3", "NRAS", "KRAS", "BRAF")
PI3K = c("IGF2", "IRS2", "PIK3CA", "PIK3R1", "PTEN")
TGFb = c("TGFBR1", "TGFBR2", "ACVR1B", "ACVR2A", "SMAD2", "SMAD3", "SMAD4")
P53 = c("TP53", "ATM")

# Some variable which will be processed by TRONCO plotting functions
pathway.genes = c(Wnt, RAS, PI3K, TGFb, P53)
pathway.names = c('Wnt', 'RAS', 'PI3K', 'TGFb', 'P53')
pathway.list = list(Wnt = Wnt, RAS = RAS, PI3K = PI3K, TGFb = TGFb, P53 = P53)

# MSS tumors  
MSS = trim(events.selection(MSS, filter.in.names = pathway.genes))
MSS = annotate.description(MSS, 'MSS subtype')

# MSI-HIGH tumors  
MSI.H = trim(events.selection(MSI.H, filter.in.names = pathway.genes))
MSI.H = annotate.description(MSI.H, 'MSI-HIGH subtype')

# We use TRONCO visualization function - oncoprint - to view these dataset 
w = oncoprint(MSS, 
	title = 'MSS tumors - with all driver genes',
	legend.cex = .5,			      		   # Legend size for events type
	gene.annot = pathway.list, 		     # List of mapping to pathways/groups
	gene.annot.color = pathways.color, # Mapping color
	sample.id = T) 					   	    	 # Sample names

w = oncoprint(MSI.H, 
	legend.cex = .5,
	gene.annot = pathway.list,
	gene.annot.color = pathways.color,
	sample.id = T)


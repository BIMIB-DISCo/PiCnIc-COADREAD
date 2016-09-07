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

# Prepare folders as we want to separate outputs
dir.create('./MSS')
dir.create('./MSI')
sub.dir = c('mutex', 'Rdata-lifted', 'Rdata-models')
sapply(paste0('./MSS/', sub.dir), dir.create)
sapply(paste0('./MSI/', sub.dir), dir.create)

# Load MAF.GISTIC file, set some fancy colors to get cute visualization later
load(paste0(workdir, '/MAF.GISTIC.Rdata'))
MAF.GISTIC = change.color(MAF.GISTIC, 'Mutation', 'darkolivegreen3')
MAF.GISTIC = change.color(MAF.GISTIC, 'Amplification', 'coral')
MAF.GISTIC = change.color(MAF.GISTIC, 'Deletion', 'cornflowerblue')
if(DOPLOTS) oncoprint(MAF.GISTIC)

#################################################################################
# Subytping.                                                                    #
#                                                                               #
# We exploit a clinical annotation provided by TCGA as a table with             #
# cluster-association data that we load, pre-process and then use to split      #
#################################################################################
file = read.delim(clusters.file, sep = ";")
head(file)

# Select  just certain annotations, remove blank lines at the end of the file
tab = file[, c("patient",  "MSI_status", "sequenced")]
tab = tab[1:276, ] 
rownames(tab) = tab$patient

# Filter out non-sequenced samples, and order them (for console visualization)
tab = tab[tab$sequenced == 1, ]
tab = tab[order(tab$MSI_status), ]
print(tab)

# Define the maps to split samples -- these are our clustering assignment
map.MSS = tab[tab$MSI_status ==  "MSS", , drop = FALSE]
map.MSI.H = tab[tab$MSI_status ==  "MSI-H", , drop = FALSE]

# These are the samples that we actually use
MSS.samples = rownames(map.MSS)
MSI.H.samples = rownames(map.MSI.H)

# TRONCO provides two functions to select a subset of samples. The one we use here
# subsets a TRONCO dataset according to a vector of names of samples. Once we 
# subset we also trim as some of the events might have 0 frequency in our sub-cohort.
# The other function which splits a dataset according to a map is ssplit(). In the 
# splitting process we get some warning as some samples for which we have clinical
# annotation are not found in the dataset -- for example TCGA-A6-3810. This 
# happens because had not both mutations and CNAs
#
# 'TCGA-A6-3810' %in% as.samples(GISTIC) # FALSE
# 'TCGA-A6-3810' %in% as.samples(MAF)    # TRUE
MSS = trim(samples.selection(MAF.GISTIC, MSS.samples))
MSI.H = trim(samples.selection(MAF.GISTIC, MSI.H.samples))
view(MSS)
view(MSI.H)

# MSS and MSI-HIGH subtypes
MSS = annotate.description(MSS, 'COADREAD - MSS subtype')
MSI.H = annotate.description(MSI.H, 'COADREAD - MSI-HIGH subtype')

# Now we start doing some fancier oncoprint visualization. We first set some colors
# also to annotate the pathways
alteration.color = 'dimgray'
pathways.color = c('firebrick1', 'darkblue', 'darkgreen', 'darkmagenta', 'darkorange')

if(DOPLOTS) oncoprint(MSS, 
	legend.cex = .5,			      		   # Legend size for events type
  cellwidth = 3,                     # Grid size
  cellheight = 10,
	gene.annot = pathway.list, 		     # List of mapping to pathways/groups
	gene.annot.color = pathways.color, # Mapping color
	sample.id = T) 					   	    	 # Sample names

if(DOPLOTS) oncoprint(MSI.H, 
          legend.cex = .5,			      		   # Legend size for events type
          cellwidth = 3,                     # Grid size
          cellheight = 10,
          gene.annot = pathway.list, 		     # List of mapping to pathways/groups
          gene.annot.color = pathways.color, # Mapping color
          sample.id = T) 					   	    	 # Sample names

# Oncoprint has a lot of options, you might want to try these visualizations as well
# 
# It can cluster samples/genes
# oncoprint(MSS, legend.cex = .5, cellwidth = 3, cellheight = 10, gene.annot = pathway.list,
#           gene.annot.color = pathways.color, sample.id = T, samples.cluster = T,
#           genes.cluster = T)
# 
# It can group genes with multiple alterations so to easily allow to spot which genes
# do not have, for instance, a CNA associated
# oncoprint(MSS, legend.cex = .5, cellwidth = 3, cellheight = 10, gene.annot = pathway.list,
#           gene.annot.color = pathways.color, sample.id = T, group.by.label = T) 					  
# 
# It can group genes with multiple alterations so to easily allow to spot which genes
# do not have, for instance, a CNA associated
# oncoprint(MSS, legend.cex = .5, cellwidth = 3, cellheight = 10, gene.annot = pathway.list,
#           gene.annot.color = pathways.color, sample.id = T, group.by.stage = T, excl.sort = F) 					  
# 
# It can highlight if some events - i.e. rows - have the same 'signature' namely occurr in the
# the same set of samples
# consolidate.data(MSI.H)
# oncoprint(MSI.H, legend.cex = .5, cellwidth = 3, cellheight = 10, gene.annot = pathway.list,
#           gene.annot.color = pathways.color, sample.id = T, annotate.consolidate.events = T) 					  
# 
# It can group samples according to a map, which helps to visualize clustering results
# tab = tab[as.samples(MAF.GISTIC), ]
# tab$patient = NULL
# tab$sequenced = NULL
# oncoprint(MAF.GISTIC, legend.cex = .5, cellwidth = 3, cellheight = 10, gene.annot = pathway.list,
#           gene.annot.color = pathways.color, sample.id = T, group.samples = tab) 					  
#
# It can show genes mapped to pathways
# oncoprint(as.pathway(MSS, pathway.genes = Wnt, pathway.name = 'Wnt'))
# oncoprint(as.pathway(MSS, pathway.genes = Wnt, pathway.name = 'Wnt', aggregate.pathway = FALSE))


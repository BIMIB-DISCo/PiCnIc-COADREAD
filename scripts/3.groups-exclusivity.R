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


# Export a file MSS.txt compliant to MUTEX input, which is one of the tools listed in
# PiCnIc table in our manuscript. We provide in TRONCO a wrapper to transform a dataset
# into the textual format compliant with MUTEX's specifications
export.mutex(MSS, 
	filename = 'MSS/mutex/MSS.txt', 
	label.mutation = 'Mutation',  
	label.amplification = 'Amplification', 
	label.deletion = 'Deletion'   
)

export.mutex(MSI.H, 
	filename = 'MSI/mutex/MSI.H.txt', 
	label.mutation = 'Mutation', 
	label.amplification = 'Amplification',
	label.deletion = 'Deletion'
	)

# MUTEX is a java tool, which is supposed to be executed outside this R
# environment. We refer to


# MUTEX output already provided. If you prefer to run MUTEX on your pc you need to
# execute it with default parameters and rename the output file as msi_result.txt
# and mss_result.txt

MSI.H.mutex = import.mutex.groups(mutex.msi.file)
MSS.mutex = import.mutex.groups(mutex.mss.file)

# Visualize groups via console
print(MSI.H.mutex)

# Then combine oncoprint via grid.arrange
require(gridExtra)
grid.arrange(
	oncoprint(
		events.selection(MSI.H, filter.in.names = MSI.H.mutex[[1]]), # Select only events for genes in group 1
				title = paste("MSI-H - Mutex group 1"),
				legend.cex = .3,
				font.row = 6,
				ann.hits = FALSE, # Avoid annotating the hits for these groups
				cellheight = 10,
				silent = T,       # Do not plot the oncoprint, just compute the GROB table
				gene.annot = pathway.list,
				gene.annot.color = pathways.color,
	)$gtable,
	oncoprint(
		events.selection(MSI.H, filter.in.names = MSI.H.mutex[[2]]), 
				title = paste("MSI-H - Mutex group 2"),
				legend.cex = .3,
				silent = T,
				font.row = 6,
				ann.hits = FALSE,
				cellheight = 10,
				gene.annot = pathway.list,
				gene.annot.color = pathways.color,
	)$gtable,
	oncoprint(
		events.selection(MSI.H, filter.in.names = MSI.H.mutex[[3]]), 
				title = paste("MSI-H - Mutex group 3"),
				legend.cex = .3,
				silent = T,
				font.row = 6,
				ann.hits = FALSE,
				cellheight = 10,
				gene.annot = pathway.list,
				gene.annot.color = pathways.color,
	)$gtable, 
	ncol=1 # Display all plots in a single column
)

# Apriori CRC knowledge
KNOWLEDGE.PRIOR.WNT = c('APC', 'CTNNB1')
KNOWLEDGE.PRIOR.RAF = c('KRAS', 'NRAS', 'BRAF')

# MEMO group estimated by TCGA  for the non-hypermutated tumors:
TCGA.MEMO = c('ERBB2', 'IGF2', 'PIK3CA', 'PTEN')

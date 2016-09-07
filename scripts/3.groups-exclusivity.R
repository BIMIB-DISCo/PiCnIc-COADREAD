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

############################################################################
# MUTEX is a Java tool, which is supposed to be executed outside this R    #
# environment, to detect groups of exclusive alteration; it is one of the  #
# tools mentioned in the PiCnIc manuscript and supported in TRONCO         #
# The original project is hosted at                                        #
#                                                                          #
# https://code.google.com/archive/p/mutex/                                 #
#                                                                          #
# and then mirrored at https://github.com/ozgunbabur/mutex                 #
#                                                                          #
# We require the user to install MUTEX on its own machine. At the time of  #
# this writing the tool is executed with the following commands            #
#                                                                          #
# java -jar target/mutex.jar MSS/mutex/MSS.txt                             #
# java -jar target/mutex.jar MSI/mutex/MSI.H.txt                           #
#                                                                          #
# where the input files are cread with TRONCO functions which provide a    #
# wrapper to transform a dataset into the textual format compliant to      #
# MUTEX's specifications                                                   #
#                                                                          #
# export.mutex(MSS,                                                        #
#             filename = 'MSS/mutex/MSS.txt',                              #
#             label.mutation = 'Mutation',                                 #
#             label.amplification = 'Amplification',                       #
#             label.deletion = 'Deletion'                                  #
# )                                                                        #
#                                                                          #
# export.mutex(MSI.H,                                                      #
#             filename = 'MSI/mutex/MSI.H.txt',                            #
#             label.mutation = 'Mutation',                                 #
#             label.amplification = 'Amplification',                       #
#             label.deletion = 'Deletion'                                  #
# )                                                                        #
#                                                                          #
# We executed the tool with its default parameters (and renamed the        #
# outputs as msi_result.txt and mss_result.txt). The execution time was    #
# approximately 6 and 3.5 hours, so we make available the MUTEX output in  #
# advance.                                                                 #
#                                                                          #
# Such results are indeed part of the files downloaded by the intial       #
# script. We provide in TRONCO functions to import results by MUTEX. Here  #
# we select only groups with score below 0.2 (default).                    #
############################################################################
MSI.H.mutex = import.mutex.groups(mutex.msi.file)
MSS.mutex = import.mutex.groups(mutex.mss.file)

# Visualize groups via console
print(MSI.H.mutex)

# Then combine oncoprint via grid.arrange utilities. As each oncoprint call
# returns a gtable GROB object we use $gtable as we want. In a similar way
# we can compute the same thing for MSS tumors.
#
# Notice that here we use a very useful function to subset data, which is
#
# events.selection(MSI.H, filter.in.names = ...)
#
# which selects all events which are associated to a gene with name in
# filter.in.names. This function has also other possible criteria to 
# select events, and will be used later.
if(DOPLOTS) grid.arrange(
	oncoprint(
		events.selection(MSI.H, filter.in.names = MSI.H.mutex[[1]]), # Select only events for genes in group 1
				title = paste("MSI-H - Mutex group 1"),
				legend.cex = .3,
				font.row = 6,
				ann.hits = FALSE, # Avoid annotating the hits for these groups
				cellheight = 10,
				silent = TRUE,       # Do not plot the oncoprint, just compute the GROB table
				gene.annot = pathway.list,
				gene.annot.color = pathways.color,
				gtable = TRUE
	)$gtable,
	oncoprint(
		events.selection(MSI.H, filter.in.names = MSI.H.mutex[[2]]), 
				title = paste("MSI-H - Mutex group 2"),
				legend.cex = .3,
				silent = TRUE,
				font.row = 6,
				ann.hits = FALSE,
				cellheight = 10,
				gene.annot = pathway.list,
				gene.annot.color = pathways.color,
				gtable = TRUE
	)$gtable,
	oncoprint(
		events.selection(MSI.H, filter.in.names = MSI.H.mutex[[3]]), 
				title = paste("MSI-H - Mutex group 3"),
				legend.cex = .3,
				silent = TRUE,
				font.row = 6,
				ann.hits = FALSE,
				cellheight = 10,
				gene.annot = pathway.list,
				gene.annot.color = pathways.color,
				gtable = TRUE
	)$gtable, 
	ncol=1 # Display all plots in a single column
)

#################################################################################
# In PiCnIc one can accommodate prior knowledge about other groups, which       #
# might have been fetched by the literature or other tools. In the COADREAD     #
# case we can exploit the fact that Wnt alterations in APC/CTNNB1 are often     #
# exclusive, as well as KRAS/RNAS/BRAF from RAF pathway. Also, TCGA has run     #
# the same analysis using another tool, called MEMO, and found significance     #
# trends of exclusivitiy among some PI3K genes and activators of that pathway.  #
# These are ERBB2, IGF2, PIK3CA and PTEN. Note that the MEMO group described by #
# TCGA contained mutations/CNAs and IGF2 overexpression. Since we are dealing   #
# with the problem of progression inference, and here we are forced to use      #
# altered states which are persistent across tumor evelution (mutations and     #
# CNAs), we are not using as a sampl IGF2 overexpression -- but its             #
# amplification yes, as it is a CNA.                                            #
#################################################################################
# Apriori CRC knowledge
KNOWLEDGE.PRIOR.WNT = c('APC', 'CTNNB1')
KNOWLEDGE.PRIOR.RAF = c('KRAS', 'NRAS', 'BRAF')

# MEMO group estimated by TCGA  for the non-hypermutated tumors:
TCGA.MEMO = c('ERBB2', 'IGF2', 'PIK3CA', 'PTEN')

if(DOPLOTS) grid.arrange(
  oncoprint(
    events.selection(MSI.H, filter.in.names = KNOWLEDGE.PRIOR.WNT),
    title = paste("MSI-H - Wnt APC/CTNNB1 exclusivity (knowledge prior)"),
    legend.cex = .3, font.row = 6, ann.hits = FALSE, cellheight = 10,
    cellwidth = 3, silent = T, gene.annot = pathway.list,
    gene.annot.color = pathways.color,
    gtable = TRUE
  )$gtable,
  oncoprint(
    events.selection(MSI.H, filter.in.names = KNOWLEDGE.PRIOR.RAF), 
    title = paste("MSI-H - RAF KRAS/NRAS/BRAF exclusivity (knowledge prior)"),
    legend.cex = .3, font.row = 6, ann.hits = FALSE, cellheight = 10,
    cellwidth = 3, silent = T, gene.annot = pathway.list,
    gene.annot.color = pathways.color,
    gtable = TRUE
  )$gtable,
  oncoprint(
    events.selection(MSI.H, filter.in.names = TCGA.MEMO), 
    title = paste("MSI-H - MEMO group (computational prior)"),
    legend.cex = .3, font.row = 6, ann.hits = FALSE, cellheight = 10,
    cellwidth = 3, silent = T, gene.annot = pathway.list,
    gene.annot.color = pathways.color,
    gtable = TRUE
  )$gtable, 
  ncol=1 # Display all plots in a single column
)

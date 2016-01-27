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

# Selected MSS tumors
tcga.f.ag.MSS.select = select(tcga.f.ag.MSS, 
	MIN.FREQ, 
	unique(c(TCGA.MEMO, 
		KNOWLEDGE.PRIOR.WNT, 
		KNOWLEDGE.PRIOR.RAF, 
		unlist(MSS.mutex)
		)
		))
tcga.f.ag.MSS.select = annotate.description(tcga.f.ag.MSS.select, 'TCGA Validation MSS colorectal tumors')

# Selected MSI tumors
tcga.f.ag.MSI.select = select(tcga.f.ag.MSI, 
	MIN.FREQ, 
	unique(c(TCGA.MEMO, 
		KNOWLEDGE.PRIOR.WNT, 
		KNOWLEDGE.PRIOR.RAF, 
		unlist(MSI.H.mutex)
		)
		))
tcga.f.ag.MSI.select = annotate.description(tcga.f.ag.MSI.select, 'TCGA Validation MSI-HIGH colorectal tumors')

# Consolidate data
cd = consolidate.data( tcga.f.ag.MSS, T)
cd = consolidate.data( tcga.f.ag.MSI, T)

# Reconstructions
tcga.f.ag.MSS.models = recon(x = tcga.f.ag.MSS.select, folder = 'VALIDATION', mutex = MSS.mutex)
tcga.f.ag.MSI.models = recon(x = tcga.f.ag.MSI.select, folder = 'VALIDATION', mutex = MSS.mutex)

# P-values for MSS training: temporal priority, probability raising and hypergeometric test
as.selective.advantage.relations(MSS.models) # edges which contribute to MLE
as.selective.advantage.relations(MSS.models, type = 'pf') # all edges 

# These should be matched to those for test. Matching depend on nodes' names, which depends for events in the node, or on
# Discarded events which were not included by consolidation etc.
# For this reason, syntactic match is not working unless for some simple cases, try: 
#
# merge(as.selective.advantage.relations(..), as.selective.advantage.relations(..), by = c('SELECTS', 'SELECTED'))
#
# In general, matching is hand-curated
as.selective.advantage.relations(tcga.f.ag.MSS.models) # edges which contribute to MLE
as.selective.advantage.relations(tcga.f.ag.MSS.models, type = 'pf') # all edges 


# P-values for MSI-HIGH
as.selective.advantage.relations(MSI.models) # edges which contribute to MLE
as.selective.advantage.relations(MSI.models, type = 'pf') # all edges 
as.selective.advantage.relations(tcga.f.ag.MSI.models) # edges which contribute to MLE
as.selective.advantage.relations(tcga.f.ag.MSI.models, type = 'pf') # all edges 

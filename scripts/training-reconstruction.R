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

select = function(x, min.freq, forced.genes) 
{
	# Collapse multiple events per gene in one unique event		
	x.sel = as.alterations(x)	
	
	# Get a list of those with minimum frequency > min.freq
	# but force inclusion of all events for genes in "forced.genes"
	x.sel = events.selection(x.sel, filter.freq = min.freq, filter.in.names = forced.genes)
		
	# Subset input - select all events for any gene in "x.sel"	
	x = events.selection(x, filter.in.names = as.genes(x.sel))
	return(x)	
}

# We set this as a variable
MIN.FREQ = 0.05

# Selected MSS dataset
MSS.select = select(MSS, 
	MIN.FREQ, 
	unique(  # All group genes
		c(TCGA.MEMO, 
		KNOWLEDGE.PRIOR.WNT, 
		KNOWLEDGE.PRIOR.RAF, 
		unlist(MSS.mutex))
		)
	)
MSS.select = annotate.description(MSS.select, 'TCGA MSS colorectal tumors')

# Selected MSI dataset
MSI.H.select = select(MSI.H, 
	MIN.FREQ, 
	unique(c(TCGA.MEMO, 
		KNOWLEDGE.PRIOR.WNT, 
		KNOWLEDGE.PRIOR.RAF, 
		unlist(MSI.H.mutex)
		)
		))
MSI.H.select = annotate.description(MSI.H.select, 'TCGA MSI-HIGH colorectal tumors')

oncoprint(MSS.select, gene.annot = pathway.list, gene.annot.color = pathways.color)

# Selected MSS is ok
cd = consolidate.data( MSS.select, T)

# In MSI we need to manipulate data instead
cd = consolidate.data( MSI.H.select, T)

# FBXW7 amplification and PTEN mutations are discarded
MSI.H.select = delete.event(MSI.H.select, gene='FBXW7', type ='Amplification')
MSI.H.select = delete.event(MSI.H.select, gene='PTEN', type ='Mutation')

recon = function(x, folder, mutex, ...) {

	lift = x
	
	# 1. create formulas from every MUTEX group
	if(!is.null(mutex))
		# TRONCO function to create  formulas from groups.
		# Notice that only 1 formula per group is created via dim.min 
		for (w in mutex) {
			lift = hypothesis.add.group(lift, 
				FUN = OR, # formula type is "soft exclusivity" (OR)
				group = w, # the group
				dim.min = length(w) # only 1 group has maximal length  
				) 
		}
	
	# Now we subset Wnt groups to genes actually present in lift, and we add again the formula.
	# Notice that formulas for groups with dimension smaller than 2 will not be created 
	KNOWLEDGE.PRIOR.WNT.subtype = KNOWLEDGE.PRIOR.WNT[
		KNOWLEDGE.PRIOR.WNT %in% as.genes(lift)
		]
	
	lift = hypothesis.add.group(lift, 
		FUN = OR, 
		group = KNOWLEDGE.PRIOR.WNT.subtype, 
		dim.min = length(KNOWLEDGE.PRIOR.WNT.subtype)) 

	# The same for the RAF group
	KNOWLEDGE.PRIOR.RAF.subtype = KNOWLEDGE.PRIOR.RAF[
		KNOWLEDGE.PRIOR.RAF %in% as.genes(lift)
		]
	
	lift = hypothesis.add.group(lift, 
		FUN = OR, 
		group = KNOWLEDGE.PRIOR.RAF.subtype, 
		dim.min = length(KNOWLEDGE.PRIOR.RAF.subtype)) 

	# And the MEMO group
	TCGA.MEMO.subtype = TCGA.MEMO[TCGA.MEMO %in% as.genes(lift)]
	
	lift = hypothesis.add.group(lift, 
		FUN = OR, 
		group = TCGA.MEMO.subtype, 
		dim.min = length(TCGA.MEMO.subtype)) 
	
	# Homologous in soft exclusivity
	lift = hypothesis.add.homologous(lift, FUN = OR)
	
	# Save to file the PDF of the lifted dataset, and its Rdata
	lift  = annotate.description(lift, as.description(x))
	oncoprint(lift, file = paste0(folder, '/Rdata-lifted/lifted.pdf'))		
	save(lift, file=paste0(folder, '/Rdata-lifted/lifted.Rdata'))	
	
	# CAPRI execution with seed set, default parameters
	model = tronco.capri(lift, boot.seed = 12345)

	# As takes long to bootstrap, we make a simple visualization first
	tronco.plot(model, 
		 pathways = pathway.list, 
		 confidence = c('tp', 'pr', 'hg'), # Display p-values 
		 ... )

	# Example non-parametric and statistical bootstrap
	model = tronco.bootstrap(model, nboot = 5)
	model = tronco.bootstrap(model, type = "statistical", nboot = 5)

	# Save the Rdata
	save(model, file=paste0(folder, '/Rdata-models/model.Rdata'))

	# Plot the TRONCO model
	tronco.plot(model, 
		 pathways = pathway.list, 
		 edge.cex = 1.5,
	 	 legend.cex = .5,
		 scale.nodes = .6,
		 confidence = c('tp', 'pr', 'hg'), # Display p-values 
		 pathways.color = pathways.color,
		disconnected = F, 
		height.logic = .3,
		file = paste0(folder, '/Rdata-models/model.pdf'),
		 ... )

	return(model)
}

# Work..
MSS.models = recon(x = MSS.select, folder = 'MSS', mutex = MSS.mutex)
MSI.models = recon(x = MSI.H.select, folder = 'MSI', mutex = MSI.H.mutex)

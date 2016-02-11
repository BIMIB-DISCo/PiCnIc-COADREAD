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

###############################################################################
# Selection                                                                   #
#                                                                             #
# We want to restrict to considering genes accorfing to their frequence. We   #
# do that for a reason: we used a TCGA curated list which ensures that those  #
# are drive, however some of them might have very low frequency (in terms of  #
# marginals) and we would like to limit their use for a model.                #
#                                                                             #
# However, some genes might be part of a group of exclusive alterations, and  #
# in that case it would be biased to account for their frequency, regardless  #
# of the group they belong to (this is also supported as many new tools try   #
# detect higly unbalanced exclusivity groups). Thus, we implement a double    #
# criterion which is in the following function.                               # 
###############################################################################
select = function(x, min.freq, forced.genes) 
{
	# First, for each gene we want to have its frequency, regardless the types
  # of alterations it harbours. So we want TP53 overall frequency by accounting
  # for mutations and CNAs.
  #
  # To get that, we collapse multiple events per gene in one unique event with
  # one of TRONCO functions
	x.sel = as.alterations(x)	
	
	# We want a list of genes which have such frequency above a threshold, which
	# is a parameter here. This list is also containig a set of genes regardless
	# of this frequency paramter - this special set of genes is a parameter as 
	# well. TRONCO events.selection supports both type of criteri.
	x.sel = events.selection(x.sel, 
	                         filter.freq = min.freq, 
	                         filter.in.names = forced.genes)
		
	# Subset input - select all events for any gene in "x.sel"	
	x = events.selection(x, filter.in.names = as.genes(x.sel))
	return(x)	
}

# We set this as a variable -- 5% minimum frequency
MIN.FREQ = 0.05

# Selected MSS dataset. The list of genes which are forced to be sleected are
# those which are part of an exclusivity group
MSS.select = select(MSS, 
	MIN.FREQ, 
	unique(  # All group genes
		c(TCGA.MEMO, 
		KNOWLEDGE.PRIOR.WNT, 
		KNOWLEDGE.PRIOR.RAF, 
		unlist(MSS.mutex))
		)
	)
MSS.select = annotate.description(MSS.select, 'TCGA COADREAD\n MSS\n colorectal tumors')

if(DOPLOTS) oncoprint(MSS.select, 
          legend.cex = .5,			      		   # Legend size for events type
          cellwidth = 3,                     # Grid size
          cellheight = 10,
          gene.annot = pathway.list, 		     # List of mapping to pathways/groups
          gene.annot.color = pathways.color, # Mapping color
          sample.id = T) 					   	    	 # Sample names

# Selected MSI dataset
MSI.H.select = select(MSI.H, 
	MIN.FREQ, 
	unique(c(TCGA.MEMO, 
		KNOWLEDGE.PRIOR.WNT, 
		KNOWLEDGE.PRIOR.RAF, 
		unlist(MSI.H.mutex)
		)
		))
MSI.H.select = annotate.description(MSI.H.select, 'TCGA COADREAD \nMSI-HIGH \ncolorectal tumors')

if(DOPLOTS) oncoprint(MSI.H.select, 
          legend.cex = .5,			      		   # Legend size for events type
          cellwidth = 3,                     # Grid size
          cellheight = 10,
          gene.annot = pathway.list, 		     # List of mapping to pathways/groups
          gene.annot.color = pathways.color, # Mapping color
          sample.id = T) 					   	    	 # Sample names

# Before the reconstruction we dataset can be consolidated, in the sense that we check
# for alterations which have the same signature across the samples. If we find them
# we could decide to either merge the events - maybe giving them a mnemonic name - 
# or delete one of them. 
#
# This is necessary as such alterations will be statistically undistinguishable.
# The function below checks also for alterations which are always/never present 
# in the samples, but it will not find those for COADREAD. If the function returns
# empty results, data is ok.

# Selected MSS is ok
cd = consolidate.data( MSS.select, T)

# In MSI we need to manipulate data instead
cd = consolidate.data( MSI.H.select, T)

# FBXW7 amplification and PTEN mutations are discarded
MSI.H.select = delete.event(MSI.H.select, gene='FBXW7', type ='Amplification')
MSI.H.select = delete.event(MSI.H.select, gene='PTEN', type ='Mutation')

# We make use of a function to transform groups in CAPRI's hypotheses,
# to add hypotheses for genes with multiple types of alterations and
# to run CAPRI
recon = function(x, folder, mutex, ...) {

	lift = x
	
	# 1. create formulas from every MUTEX group
	if(!is.null(mutex))
	  # This is TRONCO's hypothesis adding function from groups. It is the quickest way to include
	  # a formula, but not necessarily the one you want to use, as you might want something 
	  # custom, while this is pretty automatic -- and hence constrained.
	  #
	  # It takes a function FUN which is a connective from the language of n-ary AND/OR/XOR 
	  # operators, and builds a corresponding formula. We make an example, do not consider
	  # that we set dim.min = length(w).
	  #
	  # Example: consider a group of alterations in genes KRAS, TP53, APC. Assume that in our data we 
	  #   have alterations as follows: for KRAS, of both type Mutation and Deletion where no samples
	  #   harbour boths, and for TP53 and APC only of type Mutation. Let FUN = OR, also.
    #
	  #	  This function will create the following formula
	  #
	  #    OR( XOR(KRAS:Mutation, KRAS:Deletion), TP53:Mutation, APC:Mutation )
	  #
	  #   where XOR(KRAS:Mutation, KRAS:Deletion) is the subformula created as KRAS has two 
	  #   distinct types of hard-exclusive alterations. Such a formula could be equivalently
	  #   created with the hypothesis.add() function, in that case one could decide to use
	  #   only one of the KRAS alterations, for instance.
	  #
	  # This function is creating in principle as many formulas as many are the subsets of
	  # genes it can sample from the group. We do not want that in some cases, such as now,
	  # and we control that via dim.min.
    #	
		# When we set  the minimum group size (dim.min) parameter equal to the
	  # number of elements in the group, we get exactly that.
		for (w in mutex) {
			lift = hypothesis.add.group(lift, 
				FUN = OR,  # formula type is "soft exclusivity" (OR)
				group = w, # the group
				dim.min = length(w) # only 1 group has maximal length  
				) 
		}
	
	# For the sake of correctness, we check that for groups which were not explicitely 
	# computed from a cohort, we have the actual releveant genes/alteraitons available.
	#
	# Now we subset Wnt groups to genes actually present in lift, and we add again the formula.
	# Notice that formulas for groups with dimension smaller than 2 will not be created 
	KNOWLEDGE.PRIOR.WNT.subtype = KNOWLEDGE.PRIOR.WNT[
		KNOWLEDGE.PRIOR.WNT %in% as.genes(lift)
		]
	
	# This is as above
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
	
	# TRONCO also allows to automatically insert subformulas as explained in the
	# example above -- XOR(KRAS:Mutation, KRAS:Deletion). If one wants that he
	# can use this function. Here we add them in soft exclusivity via OR.
	lift = hypothesis.add.homologous(lift, FUN = OR)
	
	# Save to file the PDF of the lifted dataset, and its Rdata
	lift  = annotate.description(lift, as.description(x))
	if(DOPLOTS) oncoprint(lift, file = paste0(folder, '/Rdata-lifted/lifted.pdf'))		
	save(lift, file=paste0(folder, '/Rdata-lifted/lifted.Rdata'))	
	
	# CAPRI execution with seed set, default parameters:
	# - regularization --> AIC/BIC
	# -         pvalue --> 0.05
	# -     heuristics --> Hill Climbing
	# -          #boot --> 100
	model = tronco.capri(lift, boot.seed = 12345)

	# Save the Rdata
	save(model, file=paste0(folder, '/Rdata-models/model.Rdata'))

	# Plot the TRONCO model. TRONCO implements a visualization based on 
	# standard graph libraries, and offers various options to manipulate 
	# its basic plot. We here use some of them
	if(DOPLOTS) tronco.plot(model, 
		 pathways = pathway.list, # every node has a border annotated with pathway colors
		 edge.cex = 1.5,          # scale edge size
	 	 legend.cex = .5,         # scale legend size
		 scale.nodes = .6,        # scale node size
		 confidence = c('tp', 'pr', 'hg'), # display p-values for these statistics 
		 pathways.color = pathways.color,  # color for each pathway
		 disconnected = F,        # do not display nodes without incoming/outgoing edges
		 height.logic = .3,       # scale logical connectives
		 #file = paste0(folder, '/Rdata-models/model.pdf'), # save to file
		 ... )

	return(model)
}

# Work..
MSS.models = recon(x = MSS.select, folder = 'MSS', mutex = MSS.mutex)
MSI.models = recon(x = MSI.H.select, folder = 'MSI', mutex = MSI.H.mutex)

# Now the view function displays also some other information about the models
view(MSS.models)
view(MSI.models)

# We show the edges in each of these models through a TRONCO function which puts in a table each
# selective advantage relation, and provides also the marginal frequencies of the upstream/downstream 
# event. Also, p-values for the temporal priority, probability raising and hypergeometric testing
# are reported. Two table are returned, one for each regularizer as each model has its own set of
# edges.
as.selective.advantage.relations(MSS.models)
as.selective.advantage.relations(MSI.models)

# CAPRI has an optimization strategy which selects among all edges which satisfy Suppes' conditions
# only those which are considered best. If one is willing to see the full set of edges the above
# function can be used with a type parameter.
as.selective.advantage.relations(MSS.models, type = 'pf') # these edges are called 'prima facie' (pf)

# Notice that one can also visualize these relations in a graphical way. Try for instance
# tronco.plot(MSS.models, pf = TRUE)

# This function has some options to query only a subset of the relations, for instance one might be 
# interested in KRAS/PIK3CA ones
# as.selective.advantage.relations(MSS.models, 
#                                 events = as.events(MSS.models, genes = c('KRAS', 'PIK3CA')))
#
# Or those involving just Wnt genes
# as.selective.advantage.relations(MSS.models, 
#                                 events = as.events(MSS.models, genes = Wnt))
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
dir.create('./VALIDATION')
sub.dir = c('mutex', 'Rdata-lifted', 'Rdata-models')
sapply(paste0('./VALIDATION/', sub.dir), dir.create)

# Load mutations (binary matrix), process names
new.data = load(paste0(workdir, "/Mutations/TCGA_fun.Rdata"))
head(tcga.f.ag)
rownames(tcga.f.ag) = tcga.f.ag$gene
tcga.f.ag$gene = NULL

tcga.f.ag = t(tcga.f.ag)
rownames(tcga.f.ag) = gsub('_', '-', rownames(tcga.f.ag))

# Import in TRONCO format
tcga.f.ag = import.genotypes(tcga.f.ag, 'Mutation')
show(tcga.f.ag)

# Mutation rate to discriminate MSS and MSI-HIGH
plot(sort(rowSums(as.genotypes(tcga.f.ag))), col = 'darkgreen', lty = 3, xlab = 'Tumor sample', ylab = 'Mutations (number)')
title('Validation dataset - Mutations per sample')
abline(h = 500, col = 'red', lty = 'dashed')
abline(v = 204, col = 'black', lty = 'dashed')
text(x = 230, y = 10, 'MSI-H (Val.)', cex = .8)
text(x = 170, y = 10, 'MSS (Val.)', cex = .8)
text(x = 30, y = 600, '500 mutations (cutoff)', cex = .5)

# MSS tumors
tcga.f.ag.MSS = rownames(as.genotypes(tcga.f.ag))[
	which(rowSums(as.genotypes(tcga.f.ag)) < 500)
	] 

# MSI-HIGH tumors
tcga.f.ag.MSI = rownames(as.genotypes(tcga.f.ag))[
	which(rowSums(as.genotypes(tcga.f.ag)) > 500)
	] 	

# TRONCO objects for these tumors
tcga.f.ag.MSS = trim(samples.selection(tcga.f.ag, tcga.f.ag.MSS))
tcga.f.ag.MSI = trim(samples.selection(tcga.f.ag, tcga.f.ag.MSI))
tcga.f.ag.MSS = TCGA.shorten.barcodes(tcga.f.ag.MSS)
tcga.f.ag.MSI = TCGA.shorten.barcodes(tcga.f.ag.MSI)
show(tcga.f.ag.MSS)
show(tcga.f.ag.MSI)

# Select just genes mapped to pathways
tcga.f.ag.MSS = trim(events.selection(tcga.f.ag.MSS, filter.in.names = pathway.genes))
tcga.f.ag.MSI = trim(events.selection(tcga.f.ag.MSI, filter.in.names = pathway.genes))
tcga.f.ag.MSS = change.color(tcga.f.ag.MSS, 'Mutation', 'darkolivegreen3')
tcga.f.ag.MSI = change.color(tcga.f.ag.MSI, 'Mutation', 'darkolivegreen3')

# Load CNA data
load(paste0(workdir, 'GISTIC.Rdata'))
show(GISTIC)
GISTIC = change.color(GISTIC, 'Amplification', 'coral')
GISTIC = change.color(GISTIC, 'Deletion', 'cornflowerblue')
GISTIC = events.selection(GISTIC, filter.in.names = pathway.genes)

# Join datasets with both mutations and CNAs: MSS
tcga.f.ag.MSS = 
	ebind(
		samples.selection(tcga.f.ag.MSS, 
			intersect(
				as.samples(tcga.f.ag.MSS), 
				as.samples(GISTIC))
			),
		samples.selection(GISTIC, 
			intersect(
				as.samples(tcga.f.ag.MSS), 
				as.samples(GISTIC))
			)
		)
tcga.f.ag.MSS = trim(tcga.f.ag.MSS)
tcga.f.ag.MSS = annotate.stages(tcga.f.ag.MSS, as.stages(GISTIC))		
oncoprint(tcga.f.ag.MSS)

# Join datasets with both mutations and CNAs: MSI-HIGH
tcga.f.ag.MSI = 
ebind(
	samples.selection(tcga.f.ag.MSI, 
		intersect(
			as.samples(tcga.f.ag.MSI), 
			as.samples(GISTIC))
	),
	samples.selection(GISTIC, 
		intersect(
			as.samples(tcga.f.ag.MSI), 
			as.samples(GISTIC))
	)
)
tcga.f.ag.MSI = trim(tcga.f.ag.MSI)
tcga.f.ag.MSI = annotate.stages(tcga.f.ag.MSI, as.stages(GISTIC))		
oncoprint(tcga.f.ag.MSI)

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


##################################################################################
# Post-reconstruction Non-parametric/Statistical Bootstrap Scores                #
#   - non-parametric: we re-sample a dataset of equale size of the one used for  #
#     model reconstruction, we re-run CAPRI, and then count how many times we    #
#     re-infer the same edge as compared to how many times for each edge;        #
#   - statistical: we hold original data fixed, and we re-run the statistical    #
#     test for temporal priority and probability raising by initializing a       # 
#     random number generator with different seeds.                              #
# Results of these procedures are matrices which report the score for each pair  #
# of edges in the model.                                                         #
##################################################################################
bootstrap.stat.conf = function(model, folder, ...)
{
  # This is set to 100 for our paper, here we set a much lower value as computation time
  # might span over many hours with NUM.BOOT.ITER = 100
  NUM.BOOT.ITER = 1
  
  # Example non-parametric and statistical bootstrap
  model = tronco.bootstrap(model, nboot = NUM.BOOT.ITER)
  model = tronco.bootstrap(model, type = "statistical", nboot = NUM.BOOT.ITER)
  
  # As takes long to bootstrap, we make a simple visualization first
  tronco.plot(model, 
              pathways = pathway.list, 
              edge.cex = 1.5,          
              legend.cex = .5,         
              scale.nodes = .6,        
              confidence = c('npb', 'sb'), # display bootstrap scores 
              pathways.color = pathways.color,  
              disconnected = F,        
              height.logic = .3,       
              file = paste0(folder, '/Rdata-models/model-bootstrap.pdf'), 
              ... 
              )

  # We also print to file the bootstrap scores -- we use the pheatmap package
  require(pheatmap)
  pheatmap(keysToNames(model, as.confidence(model, conf = 'npb')$npb$bic) * 100, 
            main = paste(folder, 'COADREAD \n non-parametric bootstrap scores for BIC model'),
            fontsize_row = 6,
            fontsize_col = 6,
            display_numbers = T,
            number_format = "%d"
    )
  dev.copy2pdf(file = paste0(folder, '/Rdata-models/scores-npb-bootstrap-bic.pdf'))

  pheatmap(keysToNames(model, as.confidence(model, conf = 'npb')$npb$aic) * 100, 
            main = paste(folder, 'COADREAD \n non-parametric bootstrap scores for AIC model'),
            fontsize_row = 6,
            fontsize_col = 6,
            display_numbers = T,
            number_format = "%d"
    )
  dev.copy2pdf(file = paste0(folder, '/Rdata-models/scores-npb-bootstrap-aic.pdf'))
  
  pheatmap(keysToNames(model, as.confidence(model, conf = 'sb')$sb$bic) * 100, 
           main = paste(folder, 'COADREAD \n statistical bootstrap scores for BIC model'),
           fontsize_row = 6,
           fontsize_col = 6,
           display_numbers = T,
           number_format = "%d"
  )
  dev.copy2pdf(file = paste0(folder, '/Rdata-models/scores-sb-bootstrap-bic.pdf'))

  pheatmap(keysToNames(model, as.confidence(model, conf = 'sb')$sb$aic) * 100, 
           main = paste(folder, 'COADREAD \n statistical bootstrap scores for AIC model'),
           fontsize_row = 6,
           fontsize_col = 6,
           display_numbers = T,
           number_format = "%d"
  )
  dev.copy2pdf(file = paste0(folder, '/Rdata-models/scores-sb-bootstrap-aic.pdf'))
  
  return(model)
}

dev.new(noRStudioGD = T)

MSS.models = bootstrap.stat.conf(MSS.models, folder = 'MSS')
MSI.models = bootstrap.stat.conf(MSI.models, folder = 'MSI')

# We print to screen a table with the bootstrap scores that we just computed.
# This function has pretty much the same options of as.selective.advantage
as.bootstrap.scores(MSS.models)
as.bootstrap.scores(MSI.models)

##################################################################################
# k-fold Cross-validation                                                        #
#   - eloss: entropy loss for each model computed from 10 repetitions of 10-fold #
#     cross-validation;                                                          #
##################################################################################
mss.eloss.bic = tronco.kfold.eloss(MSS.models, 'bic')$value
mss.eloss.aic = tronco.kfold.eloss(MSS.models, 'aic')$value
msi.eloss.bic = tronco.kfold.eloss(MSI.models, 'bic')$value
msi.eloss.aic = tronco.kfold.eloss(MSI.models, 'aic')$value

# We make an example violin plot for MSS tumors
library(vioplot)
par(mfrow=c(1,2))
vioplot(mss.eloss.bic, mss.eloss.aic, col = 'red', lty = 1, rectCol="gray",
  colMed = 'black', names = c('BIC', 'AIC'), pchMed = 15, horizontal = T)
title(main = 'Entropy loss \n MSS COADREAD tumors')
vioplot(msi.eloss.bic, msi.eloss.aic, col = 'red', lty = 1, rectCol="gray",
        colMed = 'black', names = c('BIC', 'AIC'), pchMed = 15, horizontal = T)
title(main = 'Entropy loss \n MSI-HIGH COADREAD tumors')


x = c(1:5)
y = x+4
boxplot(list(x,y), main = 'ssss', notch = T)

t=as.summary.statistics(MSS.models)$aic
rownames(t) = apply(t, 1, function(x){ return(paste(x[1], "-->" , x[2])) })
t$SELECTS = NULL
t$SELECTED = NULL
str(t)
tableplot(t,
          title="Statistics for COADREAD MSS tumors -- TRONCO/PiCnIc",
          sortCol = 'NONPAR.BOOT', 
          table.label = T,
          nBins = nrow(t), 
          sample = F
          )

#devtools::install_github("gaborcsardi/crayon")
#library(crayon)
#cat(red("Hello", "world!\n"))

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

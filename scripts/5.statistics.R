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
  NUM.BOOT.ITER = 10
  
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
  if(!require(pheatmap)) install.packages('pheatmap', dependencies = T)
  pheatmap(keysToNames(model, as.confidence(model, conf = 'npb')$npb$bic) * 100, 
            main = paste(folder, 'COADREAD \n non-parametric bootstrap scores for BIC model'),
            fontsize_row = 6,
            fontsize_col = 6,
            display_numbers = T,
            number_format = "%d"
    )
  #dev.copy2pdf(file = paste0(folder, '/Rdata-models/scores-npb-bootstrap-bic.pdf'))

  pheatmap(keysToNames(model, as.confidence(model, conf = 'npb')$npb$aic) * 100, 
            main = paste(folder, 'COADREAD \n non-parametric bootstrap scores for AIC model'),
            fontsize_row = 6,
            fontsize_col = 6,
            display_numbers = T,
            number_format = "%d"
    )
  #dev.copy2pdf(file = paste0(folder, '/Rdata-models/scores-npb-bootstrap-aic.pdf'))
  
  pheatmap(keysToNames(model, as.confidence(model, conf = 'sb')$sb$bic) * 100, 
           main = paste(folder, 'COADREAD \n statistical bootstrap scores for BIC model'),
           fontsize_row = 6,
           fontsize_col = 6,
           display_numbers = T,
           number_format = "%d"
  )
  #dev.copy2pdf(file = paste0(folder, '/Rdata-models/scores-sb-bootstrap-bic.pdf'))

  pheatmap(keysToNames(model, as.confidence(model, conf = 'sb')$sb$aic) * 100, 
           main = paste(folder, 'COADREAD \n statistical bootstrap scores for AIC model'),
           fontsize_row = 6,
           fontsize_col = 6,
           display_numbers = T,
           number_format = "%d"
  )
  #dev.copy2pdf(file = paste0(folder, '/Rdata-models/scores-sb-bootstrap-aic.pdf'))
  
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
MSS.models = tronco.kfold.eloss(MSS.models)
MSI.models = tronco.kfold.eloss(MSI.models)

# One can query this statistics with this function
as.kfold.eloss(MSS.models)
as.kfold.eloss(MSI.models)

# We make an example violin plot
if(!require(vioplot)) install.packages('vioplot')

vioplot(MSS.models$kfold$bic$eloss, MSS.models$kfold$aic$eloss, col = 'red', lty = 1, rectCol="gray",
  colMed = 'black', names = c('BIC', 'AIC'), pchMed = 15, horizontal = T)
title(main = 'Entropy loss \n MSS COADREAD tumors')
#dev.copy2pdf(file = 'MSS/MSS-kfold-eloss.pdf')

vioplot(MSI.models$kfold$bic$eloss, MSI.models$kfold$aic$eloss, col = 'red', lty = 1, rectCol="gray",
        colMed = 'black', names = c('BIC', 'AIC'), pchMed = 15, horizontal = T)
title(main = 'Entropy loss \n MSI-HIGH COADREAD tumors')
#dev.copy2pdf(file = 'MSI/MSI-kfold-eloss.pdf')

##################################################################################
# k-fold Cross-validation                                                        #
#   - prederr: prediction error for each parent set X                            #
##################################################################################
MSS.models = tronco.kfold.prederr(MSS.models)
MSI.models = tronco.kfold.prederr(MSI.models)

# As above, we can query this statistics as well
as.kfold.prederr(MSS.models)
as.kfold.prederr(MSI.models)

as.selective.advantage.relations(MSS.models)$bic
as.kfold.prederr(MSS.models)$bic

##################################################################################
# k-fold Cross-validation                                                        #
#   - posterr: posterior classification error for each edge                      #
##################################################################################
MSS.models = tronco.kfold.posterr(MSS.models)
MSI.models = tronco.kfold.posterr(MSI.models)

# As above, we can query this statistics as well
as.kfold.posterr(MSS.models)
as.kfold.posterr(MSI.models)

# For instance, one can visualize a table with all edge statistics by merging all
# the table that we have produced with these functions
Reduce(
  function(...) merge(..., all = T), 
  list(
    as.selective.advantage.relations(MSS.models)$aic,
    as.bootstrap.scores(MSS.models)$aic,
    as.kfold.prederr(MSS.models)$aic,
    as.kfold.posterr(MSS.models)$aic
  )
)
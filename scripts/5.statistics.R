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
bootsrap.stat.conf = function(model, folder, ...)
{
  # This is set to 100 for our paper, here we set a much lower value as computation time
  # might span over many hours with NUM.BOOT.ITER = 100
  NUM.BOOT.ITER = 1
  
  # Example non-parametric and statistical bootstrap
  model = tronco.bootstrap(model, nboot = 2)
  model = tronco.bootstrap(model, type = "statistical", nboot = 2)
  
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

    return(model)
}

MSS.models = bootsrap.stat.conf(MSS.models, folder = 'MSS',)



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

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
  NUM.BOOT.ITER = 100
  
  # Example non-parametric and statistical bootstrap -- set cores loading
  model = tronco.bootstrap(model, nboot = NUM.BOOT.ITER, cores.ratio = .5)
  model = tronco.bootstrap(model, type = "statistical", nboot = NUM.BOOT.ITER,  cores.ratio = .5)
  
  # As takes long to bootstrap, we make a simple visualization first
  if(DOPLOTS) tronco.plot(model, 
              pathways = pathway.list, 
              edge.cex = 1.5,          
              legend.cex = .5,         
              scale.nodes = .6,        
              confidence = c('npb', 'sb'), # display bootstrap scores 
              pathways.color = pathways.color,  
              disconnected = F,        
              height.logic = .3,       
       #       file = paste0(folder, '/Rdata-models/model-bootstrap.pdf'), 
              ... 
              )

  # We also print to file the bootstrap scores -- we use the pheatmap package
  if(DOPLOTS) pheatmap(keysToNames(model, as.confidence(model, conf = 'npb')$npb$capri_bic) * 100, 
            main = paste(folder, 'COADREAD \n non-parametric bootstrap scores for BIC model'),
            fontsize_row = 6,
            fontsize_col = 6,
            display_numbers = T,
            number_format = "%d"
    )
  #dev.copy2pdf(file = paste0(folder, '/Rdata-models/scores-npb-bootstrap-bic.pdf'))

  if(DOPLOTS) pheatmap(keysToNames(model, as.confidence(model, conf = 'npb')$npb$capri_aic) * 100, 
            main = paste(folder, 'COADREAD \n non-parametric bootstrap scores for AIC model'),
            fontsize_row = 6,
            fontsize_col = 6,
            display_numbers = T,
            number_format = "%d"
    )
  #dev.copy2pdf(file = paste0(folder, '/Rdata-models/scores-npb-bootstrap-aic.pdf'))
  
  if(DOPLOTS) pheatmap(keysToNames(model, as.confidence(model, conf = 'sb')$sb$capri_bic) * 100, 
           main = paste(folder, 'COADREAD \n statistical bootstrap scores for BIC model'),
           fontsize_row = 6,
           fontsize_col = 6,
           display_numbers = T,
           number_format = "%d"
  )
  #dev.copy2pdf(file = paste0(folder, '/Rdata-models/scores-sb-bootstrap-bic.pdf'))

  if(DOPLOTS)  pheatmap(keysToNames(model, as.confidence(model, conf = 'sb')$sb$capri_aic) * 100, 
           main = paste(folder, 'COADREAD \n statistical bootstrap scores for AIC model'),
           fontsize_row = 6,
           fontsize_col = 6,
           display_numbers = T,
           number_format = "%d"
  )
  #dev.copy2pdf(file = paste0(folder, '/Rdata-models/scores-sb-bootstrap-aic.pdf'))
  
  return(model)
}

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

# If you want also the values computed at each fold just use
# as.kfold.eloss(MSS.models, values = T)

# We make an example violin plot
if(DOPLOTS) {
  vioplot(MSS.models$kfold$capri_bic$eloss, MSS.models$kfold$aic$eloss, col = 'red', lty = 1, rectCol="gray",
  colMed = 'black', names = c('BIC', 'AIC'), pchMed = 15, horizontal = T)
  title(main = 'Entropy loss \n MSS COADREAD tumors')
#dev.copy2pdf(file = 'MSS/MSS-kfold-eloss.pdf')

vioplot(MSI.models$kfold$capri_bic$eloss, MSI.models$kfold$capri_aic$eloss, col = 'red', lty = 1, rectCol="gray",
        colMed = 'black', names = c('BIC', 'AIC'), pchMed = 15, horizontal = T)
  title(main = 'Entropy loss \n MSI-HIGH COADREAD tumors')
#dev.copy2pdf(file = 'MSI/MSI-kfold-eloss.pdf')
}
##################################################################################
# k-fold Cross-validation                                                        #
#   - prederr: prediction error for each parent set X                            #
##################################################################################
MSS.models = tronco.kfold.prederr(MSS.models)
MSI.models = tronco.kfold.prederr(MSI.models)

# As above, we can query this statistics as well
as.kfold.prederr(MSS.models)
as.kfold.prederr(MSI.models)

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
tabular = function(obj, M){
  tab = Reduce(
    function(...) merge(..., all = T), 
      list(
      as.selective.advantage.relations(obj)[[M]],
      as.bootstrap.scores(obj)[[M]],
      as.kfold.prederr(obj)[[M]],
      as.kfold.posterr(obj)[[M]]
    )
  )
  # merge reverses first with second column
  tab = tab[, c(2,1,3:ncol(tab))]
  tab = tab[order(tab$NONPAR.BOOT, na.last = TRUE, decreasing = TRUE), ]
  
  return(tab)
}

tabular(MSS.models, 'capri_bic')
tabular(MSS.models, 'capri_aic')
tabular(MSI.models, 'capri_bic')
tabular(MSI.models, 'capri_bic')

# We create an Excel file with these tables
excel.file = "PicNiC-COADREAD.statistics.xlsx"

excel.wbook = createWorkbook()

sheet.mss.bic <- createSheet( wb = excel.wbook, sheetName="MSS-bic")
sheet.mss.aic <- createSheet( wb = excel.wbook, sheetName="MSS-aic")
sheet.msi.bic <- createSheet( wb = excel.wbook, sheetName="MSI-HIGH-bic")
sheet.msi.aic <- createSheet( wb = excel.wbook, sheetName="MSI-HIGH-aic")

addDataFrame(x = tabular(MSS.models, 'capri_bic'), sheet = sheet.mss.bic, showNA = T, characterNA = 'NA')
addDataFrame(x = tabular(MSS.models, 'capri_aic'), sheet = sheet.mss.aic, showNA = T, characterNA = 'NA')
addDataFrame(x = tabular(MSI.models, 'capri_bic'), sheet = sheet.msi.bic, showNA = T, characterNA = 'NA')
addDataFrame(x = tabular(MSI.models, 'capri_aic'), sheet = sheet.msi.aic, showNA = T, characterNA = 'NA')

saveWorkbook(excel.wbook, excel.file)

#### We conclude by saving two Rdata files
save(MSS.models, file='MSS/Rdata-models/model.Rdata')
save(MSI.models, file='MSI/Rdata-models/model.Rdata')
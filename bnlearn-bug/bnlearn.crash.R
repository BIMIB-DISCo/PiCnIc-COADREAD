load('net.Rdata')
head(net$data)
head(net$net)

library(bnlearn)
graphviz.plot(net$net)

# sessionInfo()
##########################
# R version 3.2.2 (2015-08-14)
# Platform: x86_64-apple-darwin13.4.0 (64-bit)
# Running under: OS X 10.9.5 (Mavericks)

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] bnlearn_3.9

# loaded via a namespace (and not attached):
# [1] tools_3.2.2         parallel_3.2.2      grid_3.2.2          Rgraphviz_2.12.0    BiocGenerics_0.14.0 stats4_3.2.2        graph_1.46.0       


# OK - some warning as some nodes have marginal 0
bn.cv(net$data, bn = net$net, loss = 'logl')
warnings()

# OK - no warnings
bn.cv(net$data, 
      bn = net$net, 
      loss = 'pred', 
      loss.args = list(target = 'Mutation PIK3CA')
      )

# OK - warnings about parents of nodes that are not present in the original data
bn.cv(net$data, 
      bn = net$net, 
      loss = 'pred-lw', 
      loss.args = list(target = 'Mutation PIK3CA', from = 'Mutation KRAS')
      )

##############################################
# crash - zero marginals?
head(net$net$arcs)
bn.cv(net$data, 
      bn = net$net, 
      loss = 'pred-lw', 
      loss.args = list(target = 'Mutation ERBB2', from = 'Amplification ERBB2'),
      fit = "bayes", 
      fit.args = list(iss = 1)
)

# *** caught segfault ***
# address 0x7fa2bd0bbdc0, cause 'memory not mapped'

# Traceback:
# 1: .Call("mappred", node = node, fitted = fitted, data = data, n = as.integer(n),     from = from, debug = debug)
# 2: map.prediction(node, fitted, data, n = n, from = from)
# 3: classification.error(node = extra.args$target, fitted = fitted, prior = extra.args$prior, n = extra.args$n, from = extra.args$from, data = data, loss = loss, debug = debug)
# 4: loss.function(fitted = net, data = data[test, ], loss = loss,     extra.args = loss.args, debug = debug)
# 5: FUN(X[[i]], ...)
# 6: lapply(kcv, bn.cv.structure, data = data, bn = bn, loss = loss,     loss.args = loss.args, fit = fit, fit.args = fit.args, debug = debug)
# 7: crossvalidation(data = data, bn = bn, loss = loss, k = k, m = m,     algorithm.args = algorithm.args, loss.args = loss.args, fit = fit, fit.args = fit.args, method = method, cluster = cluster,     debug = debug)
# 8: bn.cv(net$data, bn = net$net, loss = "pred-lw", loss.args = list(target = "Mutation ERBB2",     from = "Amplification ERBB2"))



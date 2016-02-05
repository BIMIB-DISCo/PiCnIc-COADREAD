
load('net.Rdata')
head(net$data)
head(net$net)

library(bnlearn)
graphviz.plot(net$net)

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

# crash - low marginals
bn.cv(net$data, 
      bn = net$net, 
      loss = 'pred-lw', 
      loss.args = list(target = 'Mutation ERBB2', from = 'Amplification ERBB2')
)

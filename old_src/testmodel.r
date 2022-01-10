############## TEST FILE

troncomodel <- loadRData(file='input/model_boostrapALL.rda')

tronco.plot(
  troncomodel,
  pathways = pathway.list,
  edge.cex = 1.5,
  legend.cex = .4,
  scale.nodes = .6,
  confidence = c('tp', 'pr', 'hg', 'sb', 'npb'),
  pathways.color = pathways.color,
  disconnected = F,
  height.logic = .3,
  title = "LUAD all subtypes - final model"
)
par(.pardefault)

troncomodel <- loadRData(file='input/model_boostrapTRU.rda')
tronco.plot(
  troncomodel,
  pathways = pathway.list,
  edge.cex = 1.5,
  legend.cex = .4,
  scale.nodes = .6,
  confidence = c('tp', 'pr', 'hg', 'sb', 'npb'),
  pathways.color = pathways.color,
  disconnected = F,
  height.logic = .3,
  title = "LUAD TRU subtype - final model"
)
par(.pardefault)



troncomodel <- loadRData(file='input/model_boostrapPI.rda')
tronco.plot(
  troncomodel,
  pathways = pathway.list,
  edge.cex = 1.5,
  legend.cex = .4,
  scale.nodes = .6,
  confidence = c('tp', 'pr', 'hg', 'npb', 'sb'),
  pathways.color = pathways.color,
  disconnected = F,
  height.logic = .3,
  title = "LUAD PI subtype - final model"
)
par(.pardefault)


troncomodel <- loadRData(file='input/model_boostrapPP.rda')
tronco.plot(
  troncomodel,
  pathways = pathway.list,
  edge.cex = 1.5,
  legend.cex = .4,
  scale.nodes = .6,
  confidence = c('tp', 'pr', 'hg', 'npb', 'sb'),
  pathways.color = pathways.color,
  disconnected = F,
  height.logic = .3,
  title = "LUAD PP subtype - final model"
)
par(.pardefault)


############## END TEST FILE
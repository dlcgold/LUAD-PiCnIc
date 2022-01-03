

troncomodel <- loadRData(file='input/model_boostrapALL.rda')

tronco.plot(
  troncomodel,
  pathways = pathway.list,
  edge.cex = 1.5,
  legend.cex = .4,
  scale.nodes = .6,
  confidence = c('tp', 'pr', 'hg'),
  pathways.color = pathways.color,
  disconnected = F,
  height.logic = .3,
  title = "LUAD all subtypes - capri"
)


troncomodel <- loadRData(file='input/model_boostrapTRU.rda')
tronco.plot(
  troncomodel,
  pathways = pathway.list,
  edge.cex = 1.5,
  legend.cex = .4,
  scale.nodes = .6,
  confidence = c('tp', 'pr', 'hg'),
  pathways.color = pathways.color,
  disconnected = F,
  height.logic = .3,
  title = "LUAD TRU subtype - capri"
)



troncomodel <- loadRData(file='input/model_boostrapPI.rda')
tronco.plot(
  troncomodel,
  pathways = pathway.list,
  edge.cex = 1.5,
  legend.cex = .4,
  scale.nodes = .6,
  confidence = c('tp', 'pr', 'hg'),
  pathways.color = pathways.color,
  disconnected = F,
  height.logic = .3,
  title = "LUAD PI subtype - capri"
)

troncomodel <- loadRData(file='input/model_boostrapPP.rda')
tronco.plot(
  troncomodel,
  pathways = pathway.list,
  edge.cex = 1.5,
  legend.cex = .4,
  scale.nodes = .6,
  confidence = c('tp', 'pr', 'hg'),
  pathways.color = pathways.color,
  disconnected = F,
  height.logic = .3,
  title = "LUAD PP subtype - capri"
)



install_github('hdng/clonevol')
install.packages('gridBase')
install.packages('gridExtra')
install.packages('ggplot2')
install.packages('igraph')
install.packages('packcircles')
install_github('hdng/trees')

troncomodel <- loadRData(file='input/model_boostrapTRU.rda')
statistics(troncomodel, "TRU subtype", "TRU")



install_github("phillipnicol/OncoBN")
pkgbuild::check_build_tools(debug = TRUE)

devtools::load_all()

options(buildtools.check = function(action) TRUE )
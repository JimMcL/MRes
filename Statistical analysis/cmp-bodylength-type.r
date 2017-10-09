# Comparison of body size in mimics and models
LoadFns("morpho-distance")


plotSizeVsType <- function(d = read.csv(file.path(MORPHO_DIR, "Table s3.csv"))) {
  JSetParsTemporarily(mar = c(4.5, 3, 0.2, 0.2))
  #types <- c(MTP_MODEL, MTP_SPIDER, MTP_INSECT)
  types <- MTP_NAMES
  densities <- lapply(types, function(type) density(d[d$type == type,'bodylength'], na.rm = TRUE))
  
  JPlotDensities(densities, col = typeToCol(types), lty = typeToLty(types), xlab = 'Body length (mm)', yaxt = 'n', ylab = '')
  title(ylab = 'Probability density', line = 1)
  MLegend(types, table(d$type)[types], lwd = 2, lty = typeToLty(types))
}

plotSpidersSizeVsType <- function(d = read.csv(file.path(MORPHO_DIR, "Table s2.csv"))) {
  JSetParsTemporarily(mar = c(4.5, 3, 0.2, 0.2))
  d <- d[d$order == 'Araneae',]
  types <- unique(d$type)
  densities <- lapply(types, function(type) density(d[d$type == type,'bodylength'], na.rm = TRUE))
  
  JPlotDensities(densities, col = typeToCol(types), lty = typeToLty(types), xlab = 'Body length (mm)', yaxt = 'n', ylab = '')
  title(ylab = 'Probability density', line = 1)
  MLegend(types, table(d$type)[types], lwd = 2, lty = typeToLty(types))
}

#plotSizeVsType(read.csv(file.path(MORPHO_DIR, "Table s2.csv")))
plotSizeVsType()

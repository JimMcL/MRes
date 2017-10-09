# Functions for comparing lateral to dorsal mimetic accuracy

.plotLineAndConf <- function(l, colour, style, column, range) {
  # Line
  abline(l, col = colour, lty = style)
  
  # Confidence interval
  newx <- seq(range[1], range[2], length.out=100)
  newdata <- data.frame(newx)
  names(newdata) <- column
  preds <- predict(l, newdata = newdata, interval = 'confidence', level = .95)
  polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = JTransparentColour(colour, 20), border = NA)
  lines(newx, preds[, 3], lty = 3, col = colour)
  lines(newx, preds[, 2], lty = 3, col = colour)
}


plotlateralVsDorsal <- function(d, idCol = 'species', dorsalCol = 'dorsalAccuracyMorphometric', lateralCol = 'lateralAccuracyMorphometric', main = NULL) {
  JSetParsTemporarily(mar = c(2, 2, 2, 0.2))
  # Only interested in mimics (ants are equal by definition)
  types <- c(MTP_SPIDER, MTP_INSECT)
  d <- d[d$type %in% types,]
  f <- reformulate(dorsalCol, lateralCol)
  xlim <- extendrange(d[,dorsalCol], f = 0.1)
  plot(f, data = d, main = main,
       col = typeToCol(d$type), pch = typeToPch(d$type),
       xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', xlim = xlim)
  title(xlab = 'Dorsal accuracy', ylab = 'Lateral accuracy', line = 0.5)
  abline(a = 0, b = 1, lty = 2)
  #text(d[,dorsalCol], d[,lateralCol], d$scientificName, pos = 1, cex = .8, col = typeToCol(d$type))
  
  ls <- lm(f, data = d[d$type == MTP_SPIDER,])
  .plotLineAndConf(ls, typeToCol(MTP_SPIDER), typeToLty(MTP_SPIDER), dorsalCol, xlim)
  li <- lm(f, data = d[d$type == MTP_INSECT,])
  .plotLineAndConf(li, typeToCol(MTP_INSECT), typeToLty(MTP_INSECT), dorsalCol, xlim)
  
  MLegend(types, legendPos = 'bottomright', lwd = 2)
  #print(summary(ls))
  #print(summary(li))
}

CombinedDorsalLateralPlot <- function(dir) {
  d <- read.csv(file.path(dir, "Table s3.csv"))
  PlotToPng(file.path(dir, 'dorsal-lateral.png'), function () {
    par(mfrow = c(1, 2))
    plotlateralVsDorsal(d, main = 'Morphometrics')
    plotlateralVsDorsal(d, dorsalCol = 'dorsalAccuracyMachineLearning', lateralCol = 'lateralAccuracyMachineLearning', main = 'Machine learning')
  }, width = 800, height = 400)
}

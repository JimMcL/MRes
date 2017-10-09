# Functions for testing for relaxed selection based on specimen size

LoadFns("motion-stats")
LoadFns("graphics")

# Plots mimetic accuracy for individuals or species vs size
# 
# data - Coe$fac
PlotAccuracyVsSize <- function(data, types = c(MTP_SPIDER, MTP_INSECT), formula = accuracy ~ bodylength, label = NULL, cex = NULL, ...) {
  # Extend x-axis to 0
  xlim <- range(c(0, data[data$type %in% types,]$bodylength))
  lwd = 2
  r <- MPlotScatterRegression(data, formula, types, label, 
                              xlab = 'Body length (mm)', ylab = '', 
                              yaxt = 'n', xlim = xlim, cex = cex, 
                              lwd = lwd,
                              legendPos = NULL, ...)
  MLabelAccuracyAxis(axis = 'y', cex = cex)
  legend("topright", c('All mimics', typeToLabel(types)), 
         lwd = lwd, col = c('black', typeToCol(types)), lty = c(4, typeToLty(types)), 
         cex = ifelse(is.null(cex), 1, cex), inset = JInset())
  r
}

# Plots mimetic accuracy for individuals or species vs group-centred size,
# with single regression
# 
# data - Coe$fac
PlotMimicsAccuracyVsCentredSize <- function(data, label = NULL, cex = NULL, ...) {
  mimics <- data[data$type %in% MTP_MIMICS,]
  
  # Calculate group-centred body sizes
  groupSizeMeans <- ifelse(mimics$type == MTP_SPIDER, mean(mimics$bodylength[mimics$type == MTP_SPIDER], na.rm = TRUE), mean(mimics$bodylength[mimics$type == MTP_INSECT], na.rm = TRUE))
  mimics$centredBodyLength <- mimics$bodylength - groupSizeMeans
  # Plot
  plot(accuracy ~ centredBodyLength, data = mimics, 
       col = typeToCol(mimics$type), pch = typeToPch(mimics$type), cex = cex, 
       yaxt = 'n', xlab = 'Centred body length (mm)', ylab = '')#, ...)
  MLabelAccuracyAxis(axis = 'y', cex = cex)
  if (!is.null(label))
    mtext(label, line = -1.5, adj = 0.99)
  
  # Regression
  l <- lm(accuracy ~ centredBodyLength, data = mimics)
  abline(l, lwd = 2)
  
  MLegend(legendPos = 'bottomright', types = MTP_MIMICS, pch = typeToPch(MTP_MIMICS))
  
  list(mimics = l)
}

CombinedRelaxedSelectionPlots <- function(dorsalData, lateralData) {
  .reportLm <- function(lm, type, angle) {
    sl <- summary(lm)
    m <- sl$coefficients[2,1]
    cat(sprintf("%s' %s accuracy %s with size, %s\n", typeToLabel(type), angle, 
                ifelse(m > 0, 'increases', 'decreases'), .MlmInfo(lm)))
  }

  .pr <- function(data, angle, label) {
    #lml <- PlotAccuracyVsSize(data$species$Coe$fac, label = label)
    lml <- PlotMimicsAccuracyVsCentredSize(data$species$Coe$fac, label = label)
    for(i in 1:length(lml)) 
      .reportLm(lml[[i]], names(lml)[i], angle)
  }
  
  PlotToPng(file.path(MORPHO_DIR, paste0('relaxed-selection.png')), function () {
    JSetParsTemporarily(mfrow = c(2, 1), mar = c(4.5, 3.5, .1, .1))
    .pr(dorsalData, 'dorsal', "a) Dorsal shapes")
    .pr(lateralData, 'lateral', "b) Lateral shapes")
  }, width = 450, height = 700)
}

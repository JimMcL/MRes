# Some generally useful functions


MLegend <- function(types, typeCounts = NULL, legendPos = "topright", insetFn = JInset, fill = FALSE, ...) {
  ltypes <- AsPlottableMimicTypeFactor(types)
  idx <- order(ltypes)
  if (!is.null(typeCounts)) {
    labels <- sapply(1:length(ltypes), function(i) sprintf('%s (n = %d)', typeToLabel(ltypes[i]), typeCounts[i]))
  } else {
    labels <- typeToLabel(ltypes)
  }
  bg = NA
  if (fill)
    bg <- typeToCol(ltypes[idx])
  legend(legendPos, legend = labels[idx], inset = insetFn(), col = typeToCol(ltypes[idx]), pt.bg = bg, ...)
}

# Draws a scatter plot, with appropriate symbology, and adds linear regression lines for the different types.
# returns a list of the regressions
MPlotScatterRegression <- function(data, 
                                   formula = accuracy ~ bodylength, 
                                   types = c(MTP_SPIDER, MTP_INSECT), 
                                   label = NULL, 
                                   xaxt = 's', yaxt = 's', xlab = '', ylab = '', 
                                   cex = NULL, lwd = 2, 
                                   legendPos = 'topright', extraFn = NULL, ...) {
  
  # Only interested in subset of types
  mimics <- data[data$type %in% types,]
  
  plot(formula, data = mimics, 
       col = typeToCol(mimics$type), pch = typeToPch(mimics$type), cex = cex, 
       xaxt = xaxt, yaxt = yaxt, xlab = '', ylab = '', ...)
  title(xlab = xlab, line = ifelse(xaxt == 'n', 0.5, 2.5), cex.lab = cex)
  title(ylab = ylab, line = ifelse(yaxt == 'n', 0.5, 2), cex.lab = cex)
  if (!is.null(label))
    mtext(label, line = -1.5, adj = 0.01)
  if (!is.null(extraFn))
    extraFn(mimics)

  l <- lapply(types, function(tp) {
    l <- lm(formula, data = mimics[mimics$type == tp,])
    abline(l, lwd = lwd, col = typeToCol(tp), lty = typeToLty(tp))
    l
  })
  names(l) <- types

  if (!is.null(legendPos))
    legend(legendPos, typeToLabel(types), lwd = lwd, col = typeToCol(types), lty = typeToLty(types), cex = ifelse(is.null(cex), 1, cex), inset = JInset())

  l
}

MLabelAccuracyAxis <- function(label = 'Mimetic accuracy', poorLabel = 'Poor mimics', goodLabel = 'Good mimics', axis = 'x', cex = NULL) {
  xlab = ''
  ylab = ''
  if (axis == 'x') {
    xlab = label
    side = 1
    pAdj <- 0.5
  } else {
    ylab = label
    side = 2
    pAdj = -0.5
  }
  title(xlab = xlab, ylab = ylab, line = 2, cex.lab = cex)
  mtext(poorLabel, side, adj = 0.01, padj = pAdj, cex = cex)
  mtext(goodLabel, side, adj = .99, padj = pAdj, cex = cex)
  
}

.MlmInfo <- function(lm) {
  sl <- summary(lm)
  m <- sl$coefficients[2,1]
  n <- nrow(lm$model)
  f <- sl$fstatistic[1]
  df1 <- sl$fstatistic[2]
  df2 <- sl$fstatistic[3]
  p.value <- sl$coefficients[2,4]
  r.squared <- sl$adj.r.squared
  a.r.squared <- sl$adj.r.squared
  sprintf("slope = %g, n = %d, F(%d, %d) = %g, p = %g, adjusted r2 = %g %s", 
          signif(m, 2), n, 
          df1, df2, round(f, 3),
          signif(p.value, 3), 
          signif(a.r.squared, 3), 
          PValToSigStr(p.value))
}

# Connects classes of points in a scatterplot by drawing a line from each point to the class centroid 
JStars <- function(x, y, classes, col = par('fg'), labels = FALSE, labels.col = par('fg')) {
  .star <- function(x, y, col, c) {
    centroidX <- mean(x, na.rm = TRUE)
    centroidY <- mean(y, na.rm = TRUE)
    np <- length(x)
    segments(rep(centroidX, np), rep(centroidY, np), x, y, col = col)
    if (labels)
      text(centroidX, centroidY, c, col = labels.col)
  }
  
  if(length(col) == 1)
    col <- rep(col, length(x))
  
  for (c in unique(classes)) {
    idx <- which(classes == c)
    if (length(idx) > 1)
      .star(x[idx], y[idx], col[idx], c)
  }
}

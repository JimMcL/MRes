# Functions for analysing the relationship between behavioural and static mimicry, 
# i.e. motion and body shape.

.MnReadMotionMorpho <- function(motionDir, motionCols, morphoDir, morphoCols) {
  morpho <- read.csv(file.path(morphoDir, 'Table s3.csv'))
  motion <- read.csv(file.path(motionDir, 'Table s5.csv'))
  
  # Only keep rows which are common to both  
  both <- merge(motion[,motionCols], morpho[,morphoCols])
  # Just to be nice
  both <- both[order(both$scientificName),]
  both
}

MnCmpMotionMorpho <- function(motionDir, morphoDir, dir) {
  both <- .MnReadMotionMorpho(motionDir, c('accuracy', 'scientificName'), 
                            morphoDir, c('type', 'scientificName',
                                         'dorsalAccuracyMorphometric', 'lateralAccuracyMorphometric',
                                         'dorsalAccuracyMachineLearning', 'lateralAccuracyMachineLearning'))

  .plot <- function(x, y, label, ...) {
    .pf <- function(mm) {text(mm[,x], mm[,y], mm$scientificName, cex = 0.8, pos = 1, font = 3)}
    #.pf <- function(mm) {}
    f <- reformulate(x, response = y)
    xlim <- extendrange(both[both$type == MTP_SPIDER,x], f = .16) + 0.07
    ylim <- extendrange(both[both$type == MTP_SPIDER,y], f = .04) - 0.01
    l <- MPlotScatterRegression(both, formula = f, types = MTP_SPIDER,
                                xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', 
                                label = label, legendPos = NULL, extraFn = .pf, 
                                xlim = xlim, ylim = ylim, ...)[[1]]
    MLabelAccuracyAxis('Morphological accuracy')
    MLabelAccuracyAxis(axis = 'y', 'Behavioural accuracy')
    cat(sprintf("%s: %s\n", label, .MlmInfo(l)))
  }
  
  cat(sprintf('\nBehavioural accuracy vs morphological accuracy\n'))
  PlotToPng(file.path(dir, 'motion-morpho-corr.png'), function() {
    JSetParsTemporarily(mfrow = c(1, 2), mar = c(3, 3, .1, .1))
    .plot('dorsalAccuracyMorphometric', 'accuracy', 'a) Dorsal body shape')
    .plot('lateralAccuracyMorphometric', 'accuracy', 'b) Lateral body shape')
  }, width = 800, height = 360)
  
  # PlotToPng(file.path(outDir, 'motion-morpho-gv-corr.png'), function() {
  #   JSetParsTemporarily(mfrow = c(1, 2), mar = c(.1, .1, .1, .1))
  #   .plot('dorsalAccuracyMachineLearning', 'accuracy', 'Dorsal machine learning ')
  #   .plot('lateralAccuracyMachineLearning', 'accuracy', 'Lateral machine learning')
  # }, width = 800, height = 400)
}

MnCmpSpeedMorpho <- function(motionDir, morphoDir, dir) {
  both <- .MnReadMotionMorpho(motionDir, c('mean.speed.bodylength', 'maximum.speed.bodylength', 'scientificName'), 
                            morphoDir, c('type', 'scientificName',
                                         'dorsalAccuracyMorphometric', 'lateralAccuracyMorphometric',
                                         'dorsalAccuracyMachineLearning', 'lateralAccuracyMachineLearning'))
  EILICA <- 'Eilica sp1'
  normal <- both[!(both$scientificName %in% c(EILICA)),]
  .plot <- function(xCol, yCol, label, rightLabels = c(), ...) { 
    x <- both[,c(xCol, yCol, 'type', 'scientificName')]
    names(x) <- c('x', 'y', 'type', 'scientificName')
    nx <- normal[normal$type == c(MTP_SPIDER), c(xCol, yCol, 'type')]
    names(nx) <- c('x', 'y', 'type')
    .normalLn <- function(mm) { 
      text(mm$x, mm$y, mm$scientificName, cex = 0.8, pos = ifelse(mm$scientificName %in% rightLabels, 4, 1), font = 3)
      l <- lm(y ~ x, data = nx); 
      #abline(l, col = MyPallete$lightGrey); 
      cat(sprintf('Excluding Eilica: %s\n', .MlmInfo(l))) 
    }
    xlim = extendrange(both[both$type == MTP_SPIDER, xCol], f = 0.15) - .3
    ylim = extendrange(both[both$type == MTP_SPIDER, yCol], f = 0.02)
    ll <- MPlotScatterRegression(x, y ~ x, types = MTP_SPIDER, extraFn = .normalLn, 
                                 yaxt = 'n', legendPos = NULL, xlim = xlim, ylim = ylim, ...)
    mtext(label, line = -1.5, adj = 0.97)
    MLabelAccuracyAxis(axis = 'y', 'Morphological accuracy')
    title(xlab = 'Walking speed (body lengths/sec)', line = 2)
    cat(sprintf('Including Eilica: %s\n', .MlmInfo(ll[[1]])))
  }
 
  cat(sprintf('\nMorphological accuracy vs speed\n'))
  PlotToPng(file.path(dir, 'speed-morpho-corr.png'), function() {
    JSetParsTemporarily(mfrow = c(1, 2), mar = c(3, 3, .1, .1))
    .plot('mean.speed.bodylength', 'dorsalAccuracyMorphometric', 'a) Dorsal body shape', rightLabels = c('Rhombonotus gracilis', 'Myrmarachne luctuosa', 'Iridonyssus kohouti', 'Myrmarachne macleayana'))
    .plot('mean.speed.bodylength', 'lateralAccuracyMorphometric', 'b) Lateral body shape', rightLabels = c('Judalana lutea', 'Myrmarachne luctuosa', 'Myrmarachne erythrocephala'))
  }, width = 800, height = 360)
  #MpPlotAccuracyDensities(both$type, both$mean.speed.bodylength, subset = both$type == MTP_SPIDER)
  # .plot('mean.speed.bodylength', 'dorsalAccuracyMachineLearning', 'Dorsal machine learning')
  # .plot('mean.speed.bodylength', 'lateralAccuracyMachineLearning', 'Lateral machine learning')
  nws <- normal[normal$type == MTP_SPIDER, 'mean.speed.bodylength']
  cat(sprintf('Summary of mean walking speed for mimics excluding Eilica: mean = %g, sd = %g\n', signif(mean(nws), 2), signif(sd(nws), 2)))
  cat(sprintf('Eilica mean walking speed = %g\n', signif(both[both$scientificName == EILICA, 'mean.speed.bodylength'], 2)))
}

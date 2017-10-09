# Compare imperfection results for the same photos/species using morphometric analysis and google vision
# 

LoadFns("morpho-distance")
LoadFns("sample-dbs")
LoadFns("graphics")

.readAndMerge <- function(morphoFile, gvFile, commonCol, morphoAngle = NA) {
  morpho <- read.csv(file.path(MORPHO_DIR, morphoFile))
  morpho <- morpho[morpho$outline == 'actual',]
  if (!is.na(morphoAngle))
    morpho <- morpho[morpho$angle == morphoAngle,]
  gv <- read.csv(file.path(MORPHO_DIR, gvFile))
  
  commonVals <- intersect(morpho[,commonCol], gv[,commonCol])
  df <- data.frame(commonVals, 
                   morpho[match(commonVals, morpho[,commonCol]),'imperfection'], 
                   gv[match(commonVals, gv[,commonCol]), 'imperfection'],
                   morpho[match(commonVals, morpho[,commonCol]),'mimicType'])
  names(df) <- c(commonCol, 'morpho', 'gv', 'type')
  df
}

# Data are not normally distributed since they are bounded (by 0 & 1)
# Hence lm is not appropriate
.logitRegression <- function(x) {
  x <- xy.coords(x)
  l <- glm(y ~ x, data = x, family = binomial(link = 'logit'))
  intercept <- l$coefficients[1]
  slope <- l$coefficients[2]
  # Get the backtransformed values across the whole range
  line <- data.frame(x=seq(min(x$x),max(x$x), by=0.01))
  line$y <- exp(intercept + slope * line$x) / (1 + exp(intercept + slope * line$x))
  
  list(line = line, r = l)
}

.plotCorrelation <- function(both) {
  # Data are not normally distributed since they are bounded (by 0 & 1)
  # Hence lm is not appropriate
  l <- glm(gv ~ morpho, family = binomial(link = 'logit'), data = both)
  intercept <- l$coefficients[1]
  slope <- l$coefficients[2]
  # Get the backtransformed values across the whole range
  newY <- data.frame(x=seq(min(both$morpho),max(both$morpho), by=0.01))
  newY$predictedY <- exp(intercept + slope * newY$x) / (1 + exp(intercept + slope * newY$x))
  
  # Plot it
  plot(gv ~ morpho, data=both, 
       main = "Morphometrics vs machine learning", xlab = 'Morphometrics', ylab = 'Machine learning', 
       pch=16, col = typeToCol(both$type)) #Plot underlying data
  lines(newY$x, newY$predictedY, lwd=2, col = MyPallete$darkRed) #add the curve to the plot
  #legend("topright", c("Logistic regression"), inset = c(.05, .15), col = c(MyPallete$darkRed), lwd = 2)
  MLegend(unique(both$type), sapply(unique(both$type), function(tp) sum(both$type == tp)), pch = 16)
  
  print(summary(l))
}

PhotoCorrelation <- function() {
  morpho <- read.csv(file.path(MORPHO_DIR, 'imperfection-per-photo-morpho.csv'))
  morpho <- morpho[morpho$outline == 'actual',]
  gv <- read.csv(file.path(MORPHO_DIR, 'imperfection-per-photo-gv.csv'))
  
  commonIds <- intersect(morpho$id, gv$id)
  both <- data.frame(id = commonIds, 
                     morpho = morpho[match(commonIds, morpho$id),'imperfection'], 
                     gv = gv[match(commonIds, gv$id), 'imperfection'])
  
  ###################################
  photos <- DbQueryPhotos("ptype=Photo&imageable_type=Specimen")
  specimens <- DbQuerySpecimensForPhotos(photos)
  specimens$mimicType <- mimicType(specimens) 
  specimenIds <- photos[match(both$id, photos$id),]$imageableid
  types <- AsPlottableMimicTypeFactor(specimens[match(specimenIds, specimens$id),]$mimicType)
  ##################################
  
  # Linear - wrong
  # plot(gv ~ morpho, data = both)
  # l <- lm(gv ~ morpho, data = both)
  # print(summary(l))
  # abline(l, col = MyPallete$darkRed)
  
  
  # Data are not normally distributed since they are bounded (by 0 & 1)
  # Hence lm is not appropriate
  l <- glm(gv ~ morpho, family = binomial(link = 'logit'), data = both)
  intercept <- l$coefficients[1]
  slope <- l$coefficients[2]
  # Get the backtransformed values across the whole range
  newY <- data.frame(x=seq(min(both$morpho),max(both$morpho), by=0.01))
  newY$predictedY <- exp(intercept + slope * newY$x) / (1 + exp(intercept + slope * newY$x))
  
  # Plot it
  plot(gv ~ morpho, data=both, 
       main = "Morphometrics vs machine learning", xlab = 'Morphometrics', ylab = 'Machine learning', 
       pch=16, col = typeToCol(types)) #Plot underlying data
  lines(newY$x, newY$predictedY, lwd=2, col = MyPallete$darkRed) #add the curve to the plot
  #legend("topright", c("Logistic regression"), inset = c(.05, .15), col = c(MyPallete$darkRed), lwd = 2)
  MLegend(typeList, sapply(unique(types), function(tp) sum(types == tp)), pch = 16)
  
  print(summary(l))
}

SpeciesCorrelation <- function() {
  bothd <- .readAndMerge('imperfection-per-species-angle-morpho.csv', 'imperfection-per-species-gv.csv', 'species', 'Dorsal')
  .plotCorrelation(bothd)
  bothl <- .readAndMerge('imperfection-per-species-angle-morpho.csv', 'imperfection-per-species-gv.csv', 'species', 'Lateral')
  .plotCorrelation(bothl)
}

GVAngleCorrelation <- function() {
  gv <- read.csv(file.path(MORPHO_DIR, 'imperfection-per-species-angle-gv.csv'))
  gvd <- gv[gv$angle == 'Dorsal',]
  gvl <- gv[gv$angle == 'Lateral' | gv$angle == 'Lateral right side',]
  gvl <- aggregate(gvl$imperfection, by = list(gvl$mimicType, gvl$species, gvl$bodylength), mean)
  names(gvl) <- c('mimicType', 'species', 'bodylength', 'imperfection')
  both <- merge(gvd[,cols], gvl[,cols], by = c('mimicType', 'species', 'bodylength'))
  names(both)[4] = 'dorsal'
  names(both)[5] = 'lateral'
  
  plot(dorsal ~ lateral, data = both, col = typeToCol(both$mimicType), pch = 16)
  MLegend(unique(both$mimicType), sapply(unique(both$mimicType), function(tp) sum(both$mimicType == tp)))
  
  ll <- .logitRegression(both[,c('dorsal', 'lateral')])
  lines(ll$line$x, ll$line$y, lwd=2, col = MyPallete$darkRed) #add the curve to the plot
}


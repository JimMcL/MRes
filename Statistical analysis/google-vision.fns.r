# Functions for calculating mimetic accuracy for ant mimics using the 
# google cloud vision photo classification service.
#
# The functions in this file assume the existence of a CSV file which contains 
# the result of asking google cloud vision to classify specimen photos.
# The google cloud vision query is a separate step because it costs real money, 
# so I only want to do it once.

# The location of the pre-built CSV containing the google vision query results
GVDataCSV <- JUSecureFile('photo_labels.csv', subDir = 'Classes/MRes Thesis/GoogleVisionAccuracy')

#host <- '130.56.244.191'

# Returns a data frame containing mimetic accuracy per photo.
# Mimetic accuracy is simply the probability that the subject 
# of the photo is an ant as determioned by google cloud vision.
GVPhotoMimeticAccuracy <- function(csvFile = GVDataCSV) {
  labels <- read.csv(csvFile, header = T)
  # Get a list of all tested photos so that we can differentiate between 
  # photos which aren't at all ant-like, and those which have not been tested
  result <- data.frame(id = unique(labels$id))
  # Default accuracy is 0
  result$accuracy <- 0
  
  # When known, accuracy is just the probability it's an ant
  probability <- labels[labels$label == 'ant', c('id', 'probability')]
  result[match(probability$id, result$id),'accuracy'] <- probability$probability

  result  
}

# Adds mimetic accuracy based on google visison for photos, individuals and species to the specified dataset 
GVIncorporateAccuracy <- function(data) {
  
  # Get google vision accuracy for all known photos
  acc <- GVPhotoMimeticAccuracy()
  
  # Handle complication - morpho is applied to outlines, GV to source photos.
  # Assign GV accuracy to outlines so it can be easily compared to morpho.
  # Each outline has a source value of "Photo <id>", so match up on them
  outlineIds <- match(paste('Photo', acc$id), data$photo$fac$source)
  accIds <- which(!is.na(outlineIds))
  outlineIds <- outlineIds[accIds]
  # Now acc[accIds,] <-> data$photo$fac$source[outlineIds,]
  
  # Extract out accuracy for just the photos in data
  data$photo$fac$gvAccuracy[outlineIds] <- acc$accuracy[accIds]
  
  # Exclude google vision accuracy on pickled specimens as google does a particularly bad job
  data$photo$fac$gvAccuracy[grep('pickle', data$photo$fac$state, ignore.case = TRUE)] <- NA

  # Average photos to get accuracy for individuals
  ind <- aggregate(data$photo$fac$gvAccuracy, by = list(data$photo$fac$imageableid), mean, na.rm = TRUE)
  names(ind) <- c('imageableid', 'gvAccuracy')
  data$individual$Coe$fac <- merge(data$individual$Coe$fac, ind, all.x = TRUE)
  
  # Average individuals to get accuracy for species
  sp <- aggregate(data$individual$Coe$fac$gvAccuracy, by = list(data$individual$Coe$fac$scientificName, data$individual$Coe$fac$type), mean, na.rm = TRUE)
  names(sp) <- c('species', 'type', 'gvAccuracy')
  # DON'T USE MERGE because it changes the order of the rows!!!
  #data$species$Coe$fac <- merge(data$species$Coe$fac, sp, all.x = TRUE)
  for(i in 1:nrow(data$species$Coe$fac)) {
    r <- data$species$Coe$fac[i,]
    data$species$Coe$fac[i, 'gvAccuracy'] <- sp[sp$type == r$type & sp$species == r$species, 'gvAccuracy']
  }
  
  data
}

# Plot a correlation between computer learning vs morphometric estimates of mimetic accuracy
MorphoGvCorrelation <- function(dataFac, label = NULL) {
  sf <- data.frame(sAccuracy = scale(dataFac$accuracy), sgvAccuracy = scale(dataFac$gvAccuracy), 
                   scientificName = dataFac$scientificName, type = dataFac$type)
  .plotExtras <- function(mm) {
    list[x, y] <- JInset(1, 1)
    segments(c(-x, 0), c(0, -y), c(+x, 0), c(0, y))
    #text(mm$sAccuracy, mm$sgvAccuracy, mm$scientificName, cex = .5, col = typeToCol(mm$type))
    angle <- dataFac$angle[1]
    # For dorsal, draw insect line without GV accuracy 0 species
    # if (angle == 'Dorsal') {
    #   mi <- mm[mm$type == MTP_INSECT,]
    #   #outliers <- mi$sgvAccuracy <= (min(mi$sgvAccuracy) + .0001)
    #   outliers <- grepl('Vespoidea', mi$scientificName)
    #   aa <- mi[outliers,]
    #   text(aa$sAccuracy, aa$sgvAccuracy, aa$scientificName, cex = .8, pos = 4, col = typeToCol(aa$type))
    # 
    #   no <- mi[!outliers,]
    #   l <- lm(sgvAccuracy ~ sAccuracy, data = no)
    #   abline(l, lty = 3, lwd = 2, col = typeToCol(no$type[1]))
    #   cat(sprintf("%s insects excluding outliers: %s\n", angle, .MlmInfo(l)))
    # }
  }
  r <- MPlotScatterRegression(sf, formula = sgvAccuracy ~ sAccuracy, xlab = '', ylab = '', legendPos = 'bottomright', label = label, xaxt = 'n', yaxt = 'n', asp = 1, extraFn = .plotExtras)
  MLabelAccuracyAxis(axis = 'x', label = 'Mimetic accuracy - morphometric analysis')
  MLabelAccuracyAxis(axis = 'y', label = 'Mimetic accuracy - machine learning')
  # For debugging, draw all scientific names
  # mim <- sf[sf$type %in% MTP_INSECT,]
  # text(mim$sAccuracy, mim$sgvAccuracy, mim$scientificName, cex = 0.7)
  r
}

CombinedMorphoGvCorrelation <- function(dorsalData, lateralData) {
  .reportLm <- function(lm, type, angle) {
    cat(sprintf("%s' %s morpho/GV correlation, %s\n", typeToLabel(type), angle, .MlmInfo(lm)))
  }
  
  .pc <- function(data, angle, label) {
    lml <- MorphoGvCorrelation(data$species$Coe$fac, label = label)
    for(i in 1:length(lml)) 
      .reportLm(lml[[i]], names(lml)[i], angle)
  }
  
  PlotToPng(file.path(MORPHO_DIR, paste0('gv-corr.png')), function () {
    JSetParsTemporarily(mfrow = c(2, 1), mar = c(3.5, 3.5, .1, .1))
    .pc(dorsalData, 'dorsal', "a) Dorsal shapes")
    .pc(lateralData, 'lateral', "b) Lateral shapes")
  }, width = 450, height = 700)
  
}

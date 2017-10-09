#!Rscript
#
# Exports mimetic accuracy for a photo using Google vision API.
# Note that this is not required for my masters research, but is used to provide CSV files for use by others.
#
# Reads in a CSV file containing the labels assigned to each photo by Google vision,
# then converts that to a score for each photo. 
# Currently, imperfection is just 1 - probability that photo is of an ant.

LoadFns("morpho-distance")
LoadFns("sample-dbs")
LoadFns("google-vision")

GVPlotDensities <- function(densities, title, cols) {
  JSetParsTemporarily(mar=c(2.5, 2.5, 3, 1))
  
  xlim <- range(lapply(densities, function(d) d$x))
  ylim <- range(lapply(densities, function(d) d$y))
  plot(NULL, xlim = xlim, ylim = ylim, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = title)
  mtext('Imperfection', 1, 1)
  mtext('Density', 2, 1)
  i <- 1
  for(d in densities) {
    lines(d, col = cols[i], lwd = 2)
    i <- i + 1
  }
}

GVPlotVariationDensity <- function(species, title) {
  # Average imperfection for all angles
  avg <- aggregate(species$imperfection, by = list(species$species, species$mimicType), mean)
  avg <- na.omit(avg)
  names(avg) <- c('species', 'mimicType', 'imperfection')
  
  types <- unique(avg$mimicType)
  cols = c(MyPallete$darkOrange, MyPallete$midGreen, MyPallete$midBlue, MyPallete$darkRed)
  densities <- lapply(types, function(type) density(avg[avg$mimicType == type,'imperfection'], na.rm = TRUE))
  GVPlotDensities(densities, title, cols)
  leg <- sapply(types, function(t) sprintf("%s (n = %d)", capsentence(t), sum(avg$mimicType == t)))
  legend("topright", leg, col = cols, lwd = 2)
}

# Returns info about individuals/angle/mimic type by averaging photo imperfection
GVAggregateIndividuals <- function(photos) {
  individuals <- aggregate(photos$imperfection, by = list(photos$specimen.id, photos$angle, photos$mimicType, photos$species, photos$bodylength), mean)
  names(individuals) <- c('id', 'angle', 'mimicType', 'species', 'bodylength', 'imperfection')
  individuals
}

# Exports mimetic imperfection as assessed using Google could vision API
GVExportImperfection <- function(dir) {
  imp <- GVCalcPhotoMimeticImperfection()
  names(imp) <- c('id', 'imperfection')
  # Assume undefined imperfection is very bad
  imp$imperfection[is.na(imp$imperfection)] <- 1
  # Get data about all specimen photos
  photos <- DbQueryPhotos("ptype=Photo&imageable_type=Specimen")
  # Combine with mimetic imperfection
  photos <- merge(photos, imp, all.x = TRUE)
  photos$imperfection[is.na(photos$imperfection)] <- 1
  
  # Get specimen info
  specimens <- DbQuerySpecimens()
  
  # Add specimen info to photo info
  photos <- merge(x = photos, y = specimens, by.x = 'imageableid', by.y = 'id')
  photos <- plyr::rename(photos, c("imageableid" = "specimen.id", 
                                   "description.x" = "photo.description", 
                                   "description.y" = "description",
                                   "id" = "photo.id"))
  
  # Add mimic type
  photos$mimicType <- mimicType(photos)
  
  # Assess variation between photos of the same individual
  #individualSd <- aggregate(photos$imperfection, by = list(photos$specimen.id, photos$angle, photos$mimicType, photos$species), sd)
  #names(individualSd) <- c('specimen.id', 'angle', 'mimicType', 'species', 'sd')
  #individualSd <- na.omit(individualSd)
  #plot(sd ~ as.factor(mimicType), data = individualSd)
  
  # Take means of all photos, excluding preserved specimens, for individuals
  individuals <- GVAggregateIndividuals(photos[photos$state == 'Alive' | photos$state == 'Dry',])
  
  # Average for individuals across all angles
  indAvgAngle <- aggregate(individuals$imperfection, by = list(individuals$id, individuals$mimicType, individuals$species, individuals$bodylength), mean)
  names(indAvgAngle) <- c('id', 'mimicType', 'species', 'bodylength', 'imperfection')
  
  # Take means of all individuals for a species
  species <- aggregate(data.frame(individuals$imperfection, individuals$bodylength), by = list(individuals$angle, individuals$mimicType, individuals$species), mean)
  names(species) <- c('angle', 'mimicType', 'species', 'imperfection', 'bodylength')

  # Get means for species across all angles
  spAvgAngle <- aggregate(data.frame(individuals$imperfection, individuals$bodylength), by = list(individuals$mimicType, individuals$species), mean)
  names(spAvgAngle) <- c('mimicType', 'species', 'imperfection', 'bodylength')
  
  #plot(imperfection ~ mimicType, data = species[species$angle == 'Lateral right side' | species$angle == 'Lateral',], main = 'Lateral')
  #plot(imperfection ~ mimicType, data = species[species$angle == 'Dorsal',], main = 'Dorsal')
  # GVPlotVariationDensity(species[species$angle == 'Lateral' | species$angle == 'Lateral right side',], 'Lateral')
  # GVPlotVariationDensity(species[species$angle == 'Dorsal',], 'Dorsal')
  # GVPlotVariationDensity(species, 'All angles')

  # PlotToPng(file.path(OUT_DIR, 'mimetic variation dorsal google.png'),
  #           function () GVPlotVariationDensity(species[species$angle == 'Dorsal',], 'Dorsal - Google vision'),
  #           width = 600, height = 600)
  # PlotToPng(file.path(OUT_DIR, 'mimetic variation lateral google.png'),
  #           function () GVPlotVariationDensity(species[species$angle == 'Lateral' | species$angle == 'Lateral right side',], 'Lateral - Google vision'),
  #           width = 600, height = 600)
  # PlotToPng(file.path(OUT_DIR, 'mimetic variation google.png'),
  #           function () GVPlotVariationDensity(species, 'Mimetic Imperfection as assessed by machine learning'),
  #           width = 600, height = 600)
  
  # Save imperfection per photo
  write.csv(imp, file.path(dir, 'imperfection-per-photo-gv.csv'), row.names = FALSE)
  write.csv(individuals, file.path(dir, 'imperfection-per-individual-angle-gv.csv'), row.names = FALSE)
  write.csv(indAvgAngle, file.path(dir, 'imperfection-per-individual-gv.csv'), row.names = FALSE)
  write.csv(species, file.path(dir, 'imperfection-per-species-angle-gv.csv'), row.names = FALSE)
  write.csv(spAvgAngle, file.path(dir, 'imperfection-per-species-gv.csv'), row.names = FALSE)
}

stop("Not yet implemented to use accuracy rather than imperfection")

# Where to export to - the sampleIt app so that it is accessible publically
csvExportDir <- 'C:/Jim/uni/apps/sampleit/public/ac'

# Export imperfection based on google vision
GVExportImperfection(csvExportDir)
GVExportImperfection(MORPHO_DIR)

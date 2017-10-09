# Functions to perform outline morphometric analysis on body outline images 
# Input is a set of photos and meta data, output is a fourier characterisation of
# photos, individuals, species and mimic types.
# These calculations are slow, so results are saved and reused 
# if no photos have changed since the last calculation.

suppressMessages(library(Momocs))
#library(plyr)
#LoadFns("morpho-distance")

# Loads a set of photos, converts them to outlines, subsamples and smooths them.
#
# photos - data frame which specifies the set of photos to be loaded. 
#          Must have a "file" column. Used as the "fac" for the Out object.
# sampleSize - each outline is subsampled to this number of points
# 
# value: momocs Out object
.MALoadPhotos <- function(photos, sampleSize = 1600) {
  # Convert images to outlines
  coords <- import_jpg(photos$file, verbose = FALSE)
  # Give the coords slightly more meaningful names
  names(coords) <- photos$id

  # Subsample or interpolate points to a standard number of points.
  coords <- lapply(coords, function(m) {
    if (nrow(m) < sampleSize) {
      coo_interpolate(m, sampleSize)
    } else {
      coo_sample(m, sampleSize)
    }
  })
  
  o <- Out(coords, fac = photos)
  
  # Close the outline and smooth
  coo_close(o) %>% coo_smooth(5)
}

# Performs a morphometric analysis on a set of photos
# Inputs:
# photos - data frame with a row for each photo to be processed. Columns: file is location of jpg file
# 
# Value - list with elements "photo", "individual", "species", "type".
#         photo is a momocs Coe object, remaining elements are lists with members "Coe" and "shp".
#         All elements contain both strict and lax outlines (strictOrLax column).
MorphoAnalysisForPhotos <- function(photos, startTime = NULL) {
  .st <- function(msg) {
    ShowTime(msg, startTime)
  }
    # Load photos and convert to subsampled, smoothed outlines
  outlines <- .MALoadPhotos(photos)
  .st(sprintf("Loaded %d outlines in", length(outlines)))
  
  # Align them
  pt <- proc.time()
  aligned <- fgProcrustes(outlines)
  .st(sprintf("Procrustes alignment (%g secs/shape):", round((proc.time() - pt)[3] / length(aligned$coo), 1)))

  minCoords <- min(sapply(outlines$coo, length)) / 2
  # Run the elliptical fourier analysis
  fr <- efourier(aligned, norm = T, nb.h = minCoords %/% 2)
  .st("Fourier analysis")
  # Average multiple photo outlines to individual outlines.
  # Note that specimen id is called imageableid
  individual <- mshapes(fr, 'individualId')
  # Average individuals to species
  species <- mshapes(individual$Coe, 'speciesId')
  # Average species to mimic types
  types <- mshapes(species$Coe, 'typeId')
  .st("Shape averaging")
  
  # Fix up mean values in averaged groups
  .agg <- function(specific, group, idCol, valueCol = 'bodylength', fun = mean, na.rm = TRUE) {
    r <- c()
    for(i in 1:nrow(group)) {
      r[i] <- fun(specific[specific[,idCol] == group[i,idCol], valueCol], na.rm = na.rm)
    }
    r
  }
  species$Coe$fac$bodylength <- .agg(individual$Coe$fac, species$Coe$fac, 'speciesId', 'bodylength')
  types$Coe$fac$bodylength <- .agg(species$Coe$fac, types$Coe$fac, 'speciesId', 'bodylength')
  

  list(photo = fr, individual = individual, species = species, type = types)
}

MpLDACalcAccuracy <- function(coe, priorWeights = c(1, 1)) {
  # Start with a PCA to reduce the number of dimensions and eliminate constant features
  p <- PCA(coe)
  # Now perform the LDA, attempt to discriminate ants from non-ants
  l <- LDA(p, factor(ifelse(p$fac$type == MTP_MODEL, 'ant', 'non-ant')), verbose = FALSE, priors = priorWeights / length(priorWeights))
  
  # Just to simplify life...
  #l$accuracy <- l$mod.pred$posterior[,'ant']
  l$accuracy <- scale(-l$mod.pred$x[,1])
  
  # To plot: MpPlotAccuracyDensities(data$species$Coe$fac, l$accuracy)
  l
}

# Reads in a precalculated morphometric analysis, checks if it matches the photos, and if not, recalculates it.
# All photos are assumed to have been taken from the same angle (photos$angle).
# If photos are not specified, precalculated results are unconditionally returned (and angle must be specified).
#
# photos - data frame describing the photos to be analysed. 
#          If NULL, the cached analysis is returned (angle must be specified in this case).
# angle - viewing angle of photos to be analysed.
# see MorphoAnalysisForPhotos for return value
GetMorphoForPhotos <- function(photos, angle = photos[1,]$angle, verbose = TRUE, force = FALSE) {

  # Try to read in the result of the last analysis
  analFile <- JUDFile(paste0(angle, '.rds'), 'morpho', createDir = TRUE)
  result <- tryCatch(readRDS(analFile), 
                     error = function(c) NULL,
                     warning = function(c) NULL)

  # Did we manage to read in an up-to-date file? Assume photos can't change unless their ids change
  #cat(sprintf("Null result? %s, null photos? %s, idsets equal? %s, verbose? %s\n", is.null(result), is.null(photos), setequal(result$photo$id, photos$id), verbose))
  if (force || is.null(result) || (!is.null(photos) && !setequal(result$photo$id, photos$id))) {
    if (verbose)
      cat(sprintf("Recalculating morphometrics for %s photos\n", angle))
    # Something has changed. Perform morphometric analysis on all the photos
    startTime <- if(verbose) proc.time() else NULL

    result <- MorphoAnalysisForPhotos(photos, startTime)
    saveRDS(result, analFile)
    
    ShowTime("Total calculation time", startTime)
  }
  
  result
}

# see MorphoAnalysisForPhotos for return value
MGetData <- function(angle, checkUpToDate = !JCLIOption('g', 'dontcheck')) {
  cat(sprintf("Angle %s, checkUpToDate %s\n", angle, checkUpToDate))
  photos <- NULL
  if(checkUpToDate == 'force' || checkUpToDate) {
    photos <- DbGetOutlines(angle)
  }
  GetMorphoForPhotos(photos, angle = angle, force = checkUpToDate == 'force')
}


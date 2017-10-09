# Functions for analysing distances in morphological space

MORPHO_DIR <- file.path(JUCFile("Thesis"), 'Morphometrics')

# Returns variance and centroid of a set of points
getGroupInfo <- function(points, groupType) {
  # Centroid of a set of points is just the mean
  centroid <- apply(points, 2, mean)
  
  # Calculate distance between all pairs of points
  dists <- data.frame(dist = as.vector(dist(points)), type = groupType)
  
  list(centroid = centroid, dist = dists)
}

# Gathers sets of points into groups with equal fac[,column] values
# points - set of multidimensional points such as the x value from fourier or PCA analysis
# fac - factors for each row in points
# column - name of column within fac which is used to group rows
# type - name of the type of grouping, e.g. 'within individuals', 'within species'...
AnalyseGroups <- function(points, fac, column, type) {
  
  #cat(sprintf("AnalyseGroups %s, %s\n", column, type))
  values <- unique(fac[,column])
  # cat(sprintf("Analyse %s: %d points, %d fac, %d unique values\n", type, nrow(points), nrow(fac), length(values)))
  
  groupPoints <- matrix(nrow = 0, ncol = ncol(points))
  groupFac <- data.frame()
  groupDists <- data.frame()
  
  for(val in values) {
    idx <- which(fac[,column] == val)
    #cat(sprintf(" - val %s, indices %s\n", val, paste(idx, collapse = ", ")))
    if (length(idx) == 1) {
      vpoints <- points[idx,]
      vfac <- fac[idx,]
    } else {
      list[vpoints, dists] <- getGroupInfo(points[idx,], type)
      groupDists <- rbind(groupDists, dists)
      # Use the first row of group facs to represent the group
      vfac <- fac[idx[1],]
    }
    
    # Save centroid as position of this group
    groupPoints <- rbind(groupPoints, vpoints)
    
    # Save fac value for this group
    groupFac <- rbind(groupFac, vfac)
  }
  
  list(points = groupPoints, fac = groupFac, dists = groupDists)
}

# Calculates distances between groups of points
AnalyseVariance <- function(points, fac, columns, types) {
  result <- data.frame()
  for(i in 1:length(columns)) {
    list[points, fac, dists] <- AnalyseGroups(points, fac, columns[i], types[i])
    result <- rbind(result, dists)
  }
  result
}

GenerateRandomFacData <- function(numObs, speciesPerGenus, numSpecies, photosPerIndividual) {
  
  individualsPerSpecies <- numObs %/% (numSpecies * photosPerIndividual)
  
  r <- data.frame()
  for(sp in 1:numSpecies) {
    for(ind in 1:individualsPerSpecies) {
      for(ph in 1:photosPerIndividual){
        r <- rbind(r, list(kingdom = 'animal', genus = round(sp / speciesPerGenus), speciesId = sp, individualId = 10 * sp + ind))
      }
    }
  }
  r
}

# Generates density plots for distance between photos of individuals, individuals within a species, and between species
PlotMorphoDistance <- function(points, fac, title) {
  oldPars <- par(no.readonly = T)
  on.exit(par(oldPars))
  
  columns <- c("individualId", "speciesId", "genus", "kingdom")
  labels <- c("Within individuals", "Within species", "Within genera", "Between genera")
  av <- AnalyseVariance(points, fac, columns, labels)
  colours <- c(MyPallete$darkBlue, MyPallete$darkGreen, MyPallete$darkRed, MyPallete$darkOrange)
  
  dense <- function(idx) {
    rows <- av$type == labels[idx]
    if(sum(rows) > 1)
      density(av[rows, 'dist'])
    else
      NULL
  }
  indices <- 1:length(columns)
  legend <- sapply(indices, function(idx) sprintf("%s (n = %d)", labels[idx], nrow(av[av$type == labels[idx],])))

  ymax <- max(sapply(indices, function(idx) {
    d <- dense(idx)
    ifelse(is.null(d), 0, max(d$y))
  }))
  plot(NULL,
       xlab = "Distance",
       xlim = c(0, ceiling(max(av$dist))), 
       ylim = c(0, ymax),
       main = title)
  styles <- (indices + 1) %/% 2
  weight <- 2
  for(idx in indices) {
    lines(dense(idx), lwd = weight, lty = styles[idx], col = colours[idx])
  }
  legend("topright", legend=legend, col = colours, lty = styles, lwd = weight)
}

# Attempt to quantify the distribution of mimetic variance
# modelXy - points repreenting models in morphospace
# mimicXy - points representing mimics in morphospace
# modelIdxForMimic - index into modelXy for each row in mimicXy
# 
# value - set of distances from each row of mimicXy to the appropriate row of modelXy, 
# i.e. value[i] = dist(mimicXy[i,], modelXy[modelIdxForMimic[i],])
CalcMimeticVariation <- function(modelXy, mimicXy, modelIdxForMimic) {
  r <- numeric(length = nrow(mimicXy))
  for (i in 1:nrow(mimicXy)) {
    r[i] <- dist(rbind(mimicXy[i,], modelXy[modelIdxForMimic[i],]))
  }
  r
}

PlotMimeticVariation <- function(modelXy, mimicXy, modelIdxForMimic, add = FALSE, ...) {
  den <- density(CalcMimeticVariation(modelXy, mimicXy, modelIdxForMimic))
  # This is my own hack to attempt to correct for non-normal distribution - I think I am converting Rayleigh distribution to normal distribution
  #den$y <- den$y / (den$x * 2 * pi)
  if (add)
    lines(den, ...)
  else
    plot(den, ...)
}

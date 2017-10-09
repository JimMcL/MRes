#!Rscript
# 
# MRes Morphometric analysis

library(colorspace)
LoadFns("morpho-calculation")
LoadFns("sample-dbs")
LoadFns("morpho-distance")
LoadFns("morpho-export")
LoadFns("morpho-relaxed-selection")
LoadFns("motion")
LoadFns("google-vision")
LoadFns("cmp-lateral-dorsal")

# Queries SampleIt for photos.
# value - dataframe containing information about the photos, including name of downloaded file
# sampleSize - each outline is subsampled to this number of points
LoadPhotos <- function(query, sampleSize = 1600) {
  # Get the list of matching photos
  photos <- read.csv(paste0("http://localhost/photos.csv?", query), 
                    stringsAsFactors = F, strip.white=TRUE)
  # Download the photos
  file = JDownload(photos$url, verbose = F)

  ## Read specimen info for each photo
  specimenIds <- photos[photos$imageabletype == 'Specimen',]$imageableid
  facQuery <- sprintf("id=[%s]", paste(unique(specimenIds), collapse=','))
  specimens <- read.csv(paste0("http://localhost/specimens.csv?", facQuery), 
                        stringsAsFactors = T, strip.white=TRUE)
  fac <- specimens[match(photos$imageableid, specimens$id),]
  # Remove columns with duplicate names
  fac2 <- fac[,!colnames(fac) %in% c('id', 'description')]
  
  ## Create associated info
  sOrL <- strictOrLax(photos)
  mt <- mimicType(fac)
  p = cbind(photos, file, 
            type = mt, 
            label = sprintf("%s, %s", photos$imageableid, photos$id), 
            individualId = individualIds(photos, sOrL, mt),
            speciesId = speciesIds(fac, sOrL, mt),
            strictOrLax = sOrL,
            fac2,
            stringsAsFactors = F)
  coords <- import_jpg(p$file)
  # We don't really care about the file names, so replace names (which are file names) with photo ids
  names(coords) <- photos$id

  # Subsample or interpolate points to a standard number of points.
  coords <- lapply(coords, function(m) {
    if (nrow(m) < sampleSize) {
      coo_interpolate(m, sampleSize)
    } else {
      coo_sample(m, sampleSize)
    }
    })
  
  o <- Out(coords, fac = p)
  
  # Close the outline and smooth
  coo_close(o) %>% coo_smooth(5)
}

lightenColour <- function(col, lightness.factor) {
  # Convert to RGB
  rgb <- col2rgb(col)
  # Convert to HSL
  hsl <- rgb2hsv(rgb)
  # Adjust
  hsl['v',] <- pmin(1, hsl['v',] + (1 - hsl['v',]) * lightness.factor)
  hsl['s',] <- pmax(0, hsl['s',] * (1 - lightness.factor))
  apply(hsl, 2, function(col) do.call(hsv, as.list(col)))
}

# labels - column name for labelling points, or FALSE
# labels.col - colour for group labels. NA means use same colour as points. NULL means no labels
# legend - if if true, draw a legend
# bareAxes - plot only axes with morphospace shapes
plotMorphoPca <- function(pca, 
                    cex = 1.5, 
                    points = TRUE,
                    labels = "label", 
                    labels.cex = .8, 
                    labels.col = NA, 
                    ellipses = TRUE,
                    ellipsesax = TRUE,
                    legend = TRUE, 
                    legendInsetFn = JInset,
                    legendYIntersp = 1,
                    bareAxes = FALSE, 
                    extraPlottingFn = NULL,
                    zoom = 1) {
  
  # make shapes intense for bare axes,faded when plotting data
  list[col.shp, border.shp] <- lightenColour(c(MyPallete$lightBlue, '#000000'), ifelse(bareAxes, 0, 0.6))
  
  # Strict outlines are triangles, lax are circles
  pch <- typeToPch(pca$fac$type)

  plot(pca, "type", 
       points = points && !bareAxes, 
       pch = pch, 
       cex = cex, 
       col = typeToCol(levels(pca$type)),
       center.origin = FALSE, 
       zoom = zoom,
       ellipses = ellipses && !bareAxes, 
       ellipsesax = ellipses && !bareAxes && ellipsesax,
       labelsgroups = !bareAxes && !is.null(labels.col) && is.na(labels.col), 
       cex.labelsgroups = labels.cex,
       rect.labelsgroups = FALSE,
       morphospace = TRUE,
       size.shp = 2,
       pos.shp = 'full_axes',
       border.shp = border.shp,
       col.shp = col.shp,
       lwd.shp = 2,
       labelspoints = FALSE,
       eigen = FALSE,
       rug = FALSE,
       title = '',
       axisnames = FALSE, axisvar = FALSE)
  
  if (!bareAxes) {
    if (labels != FALSE) {
      # Don't use labelspoints in plot.PCA so we have more control over text, e.g. relative position
      text(pca$x[,c(1,2)], labels = pca$fac[, labels], pos = 1, cex = labels.cex, col = typeToCol(pca$fac$type))
    }
    
    # If label colour was specified, label groups explicitly, as plot.pca doesn't allow group label colour to be specified
    if (!is.null(labels.col) && !is.na(labels.col))
      .labelsgroups(pca$x[,c(1,2)], pca$fac[, 'type'], labels.col, cex = labels.cex)
    
    if (legend) {
      types <- c(MTP_MODEL, MTP_SPIDER, MTP_INSECT, MTP_NON_MIMIC)
      # Work around a bug - the top border of the legend is sometimes clipped, even though it's actually inside the plot area
      # Unfortunately, this messes up the coordinate system
      oldPar <- par(xpd=TRUE)
      legend("topright", 
             legend=typeToLabel(types), 
             pt.cex = cex,
             cex = labels.cex,
             col = typeToCol(types),
             pt.bg = typeToCol(types),
             pch = typeToPch(tolower(types)),
             y.intersp = legendYIntersp,
             inset = legendInsetFn())
      par(oldPar)
    }
  }
  
  if (!is.null(extraPlottingFn)) {
    extraPlottingFn(colours)
  }
}

HighlightPcaPoints <- function(pca, idx, radius = 1, labels = NULL) {
  symbols(x = pca$x[idx,1], y = pca$x[idx,2], circles = rep(radius, length(idx)), inches = F, add = T)
  if (!is.null(labels)) {
    text(pca$x[idx,1], pca$x[idx,2], labels = pca$fac[idx, labels], pos = 4, offset = 1)
  }
}

# Mean shapes
OverlayMeanShapes <- function(f) {
  ms <- f %>% mshapes("type")
  ms$shp[,MTP_MODEL] %>% coo_plot(border = "purple")
  ms$shp[,MTP_INSECT] %>% coo_draw(border = "green")
  ms$shp[,MTP_SPIDER] %>% coo_draw(border = "red")
  ms$shp[,MTP_NON_MIMIC] %>% coo_draw(border = "grey")
  legend("topright", lwd=1, col=c("purple", "green", "red", "grey"), legend=typeToLabel(MTP_MODEL, MTP_INSECT, MTP_SPIDER, MTP_NON_MIMIC))
}

# Internal fn
plotInverse <- function(s, f, nharmonics, title = paste("Harmonics 0 to", nharmonics)) {
  inverse <- efourier_i(f, nharmonics, 64)
  plot(s, type="l", asp=1, main = title, col="grey", xlab = '', ylab = '')
  polygon(s, col="grey", border=NA)
  lines(inverse[,1], inverse[,2], type="l")
}

# Shows accuracy of shape specified with differing numbers of harmonics
PlotHarmonics <- function(shape, index, nHarmonics = 9) {
  # From Morphometrics with R

  shape <- shape[[index]]

  f <- efourier(shape, norm = T, nb.h = 60)

  nrow <- floor(sqrt(nHarmonics))
  ncol <- ceiling(nHarmonics / nrow)
  JSetParsTemporarily(mfrow=c(nrow, ncol), mar=c(3,4,3,1))

  for (i in 1:nHarmonics) {
    plotInverse(shape, f, i)
  }
}

# Shows accuracy of fourier characterisation with 9 harmonics for multiple shapes
Plot9thHarmonic <- function(shape, titles, indices, harmonics = 9) {
  oldPars <- par(no.readonly = T)
  on.exit(par(oldPars))

  nrow <- floor(sqrt(length(indices)))
  ncol <- ceiling(length(indices) / nrow)
  par(mfrow=c(nrow, ncol))
  par(mar=c(3, 2, 3, 1))
  
  for(idx in indices) {
    # numeric idx gives different result than integer idx!!!
    s <- shape[as.numeric(idx)]
    f <- efourier(s, norm = T, nb.h = 64)
    plotInverse(s, f, harmonics, title = titles[idx])
  }
}

# Plots morphological distance variations to a PNG file
PlotDistanceVariations <- function(coe, filter, title) {
  points <- coe$x[filter,]
  fac <- coe$fac[filter,]
  
  pngWidth <- 800
  pngHeight <- 800
  
  PlotToPng(file.path(MORPHO_DIR, paste0('morpho-distances-', angle, title, '.png')),
            function() { 
              title <- paste0(capwords(angle), ' Morphometric Distances', title)
              PlotMorphoDistance(points, fac, title) 
            },
            width = pngWidth, height = pngHeight)
}

VisualiseVariation <- function(pf, column) {
  plot(pf, column, morphospace = TRUE, labelsgroups = FALSE, ellipsesax=FALSE, pos.shp = 'full_axes')

  # Draw stars (can't use stars param to plot because many groups contain only 1 point)  
  # Can only plot variation of groups with more than 1 row
  multiples <- table(pf$fac[,column])
  multiples <- multiples[multiples > 1]
  for(id in names(multiples)) {
    idxs <- which(id == pf$fac[,column])
    xy <- pf$x[idxs, 1:2]
    centroid <- apply(xy, 2, mean)
    i <- grep(id, names(multiples))
    col <- MyPallete[[i %% length(MyPallete) + 1]]
    lwd <- 1
    segments(xy[,1], xy[,2], centroid[1], centroid[2], col = col, lwd = lwd)
    if (max(dist(xy)) > .04) {
      # print(sprintf("%d points, dist %g (%g), centroid %g, %g, id %s", length(idxs), max(dist(xy)), max(dist(pf$x[idxs,])), centroid[1], centroid[2], id))
      text(centroid[1], centroid[2], labels = id, col = col)
    }
  }
}

# This produces a plot of "mimetic variation" for spider mimics, insect mimics, and models,
# where "mimetic variation" is calculated as distance from the mean shape of all models
plotVariation <- function(speciesShapes, typeShapes, simple = FALSE,
                          includeTypes = NA, lwd = 2, 
                          title = 'Species-level Mimetic Variation', cex.main = 1.2,
                          leg.cex = 1, legBg = rgb(.95, .95, .95)) {
  # Ignore lax type shapes
  mi <- which(typeShapes$Coe$fac$type == MTP_MODEL & typeShapes$Coe$fac$strictOrLax == 'strict')
  
  lineStyle <- c(lax = 2, strict = 1)
  types <- sort(AsPlottableMimicTypeFactor(unique(speciesShapes$Coe$fac$type)))
  cols <- typeToCol(types)
  # Get indices of lax and strict outlines
  laxIdx <- grepl('lax', speciesShapes$Coe$fac$description, ignore.case = TRUE)
  laxTypes <- list(lax = laxIdx, strict = !laxIdx)
  leg <- character(0)
  legCols <- numeric(0)
  legStyle <- numeric(0)
  
  # Calculate density for each type
  densityInfo <- list()
  for (i in 1:length(types)) {
    typeIdx <- speciesShapes$Coe$fac$type == types[i]
    for(strict in c('lax', 'strict')) {
      # Do we want to plot this line?
      if (is.na(includeTypes) || (paste(types[i], strict) %in% includeTypes)) {
        spIdx <- typeIdx & laxTypes[[strict]]
        specimens <- speciesShapes$Coe$coe[spIdx,]
        if (sum(spIdx) > 1) { # ??? Why doesn't it work for 1 specimen???
          density <- density(CalcMimeticVariation(typeShapes$Coe$coe, specimens, rep(mi, nrow(specimens))))
          densityInfo[[length(densityInfo) + 1]] <- list(density = density,
                                                         col = cols[i],
                                                         lty = lineStyle[strict])
          legText <- ifelse(simple, 
                            sprintf("%s (%d species)", capsentence(paste0(types[i], "s")), nrow(specimens)),
                            sprintf("%s %s (n = %d)", capwords(types[i]), strict, nrow(specimens)))
          leg <- c(leg, legText)
          legCols <- c(legCols, cols[i])
          legStyle <- c(legStyle, lineStyle[strict])
        }
      }
    }
  }

  xlim <- range(sapply(densityInfo, function(di) range(di$density$x)))
  ylim <- range(sapply(densityInfo, function(di) range(di$density$y)))

  # Bare plot
  plot(NULL, 
       xlim = xlim, ylim = ylim, 
       main = title, cex.main = cex.main,
       xaxt = 'n', yaxt = 'n',
       xlab = 'Mimic imperfection', ylab = 'Density')
  # Add lines
  invisible(lapply(densityInfo, function(di) lines(di$density, col = di$col, lty = di$lty, lwd = lwd)))
  
  #leg <- sapply(types, function (type) sprintf("%s (n = %d)", capwords(type), sum(speciesShapes$Coe$fac$type == type)))
  legend("topright", legend = leg, inset = .018, col = legCols, lty = legStyle, lwd = lwd, bg = legBg, cex = leg.cex)
  
}

PlotSpeciesPCAWithHighlights <- function(data, hiIdxFn, hiLabels = "scientificName") {
  pca <- PCA(data$species$Coe)
  hiFn <- function(colours) {
    idxs <- hiIdxFn(pca$fac)
    print(idxs)
    HighlightPcaPoints(pca, idxs, labels = hiLabels, radius = .01)
  }
  plotMorphoPca(pca, labels = "scientificName", extraPlottingFn = hiFn, ellipsesax = FALSE)
}

reportMorphoStatistics.obsolete <- function(speciesShapes, typeShapes) {
  
  ##### NOTE NOTE NOTE
  # I think that these are Rayleigh distributions which violates the t.test and var.test assumptions of normality (I think)
  # Need a multi-dimensional test, see package ks?
  mi <- which(typeShapes$Coe$fac$type == MTP_MODEL & typeShapes$Coe$fac$strictOrLax == 'strict')
  
  dists <- CalcMimeticVariation(typeShapes$Coe$coe, speciesShapes$Coe$coe, rep(mi, length(speciesShapes$Coe$fac$type)))
  data <- cbind(speciesShapes$Coe$fac, dists)
  print(" === spider vs insect mimics ===")
  print(t.test(data[data$type == MTP_SPIDER,]$dists, data[data$type == MTP_INSECT,]$dists, alternative = 'g'))
  print(var.test(data[data$type == MTP_SPIDER,]$dists, data[data$type == MTP_INSECT,]$dists, alternative = 'g'))
  print(" === model vs spider mimics ===")
  print(t.test(data[data$type == MTP_MODEL,]$dists, data[data$type == MTP_SPIDER,]$dists))
  print(" === model vs insect mimics ===")
  print(t.test(data[data$type == MTP_MODEL,]$dists, data[data$type == MTP_INSECT,]$dists))
}

reportMorphoStatistics <- function(dorsalData, lateralData, accuracyCol = 'accuracy', reportWorseThanModels = FALSE) {
  
  dd <- dorsalData$species$Coe$fac
  dl <- lateralData$species$Coe$fac
  # dd <- dorsalData$individual$Coe$fac
  # dl <- lateralData$individual$Coe$fac
  alpha <- 0.05
  
  
  .fr <- function(rng) paste(rng, collapse = ' to ')
  all <- rbind(dorsalData$individual$Coe$fac[,c('type', 'bodylength')],
               lateralData$individual$Coe$fac[,c('type', 'bodylength')])
  all <- na.omit(all)
  cat(sprintf("Specimens range from %s mm\n", .fr(range(all$bodylength))))
  cat(sprintf("Spider mimics range from %s mm\n", .fr(range(all[all$type == MTP_SPIDER, 'bodylength']))))
  cat(sprintf("Insect mimics range from %s mm\n", .fr(range(all[all$type == MTP_INSECT, 'bodylength']))))
  
  cat(sprintf("No. of shapes\tModels\tMimetic spiders\tMimetic insects\tNon-mimics\n"))
  t <- table(dd$type)
  cat(sprintf("Dorsal\t\t%d\t%d\t%d\t%d\n", t[MTP_MODEL], t[MTP_SPIDER], t[MTP_INSECT], t[MTP_NON_MIMIC]))
  t <- table(dl$type)
  cat(sprintf("Lateral\t\t%d\t%d\t%d\t%d\n", t[MTP_MODEL], t[MTP_SPIDER], t[MTP_INSECT], t[MTP_NON_MIMIC]))
  cat('Body lengths by type\n')
  print(t(sapply(MTP_NAMES, function(tp) {
    d <- dorsalData$species$Coe$fac
    bl <- d[d$type == tp, 'bodylength']
    c(mean = mean(bl, na.rm = TRUE), sd = sd(bl, na.rm = TRUE))
  })) %>% signif(2))
  
  .meanDiff <- function(d, angle, col = accuracyCol, leftTypes = MTP_SPIDER, rightType = MTP_INSECT) {
    tt <- t.test(d[d$type %in% leftTypes, col], d[d$type == rightType, col], alternative = 'l')
    cat(sprintf("%s mean lower? %s; t = %g, df = %g, p = %g %s\n", angle, ifelse(tt$p.value < alpha, 'Yes', 'No'), round(tt$statistic, 3), tt$parameter, round(tt$p.value, 4), PValToSigStr(tt$p.value)))
  }
  .varDiff <- function(d, angle, col = accuracyCol, leftTypes = MTP_SPIDER, rightType = MTP_INSECT) {
    ft <- var.test(d[d$type %in% leftTypes, col], d[d$type == rightType, col], alternative = 'l')
    cat(sprintf("%s variance lower? %s; F(%d, %d) = %g, p = %g %s\n", angle, ifelse(ft$p.value < alpha, 'Yes', 'No'), ft$parameter[1], ft$parameter[2], round(ft$statistic, 3), round(ft$p.value, 3), PValToSigStr(ft$p.value)))
  }

  # Don't care about this  
  # cat(sprintf("--- Are mimetic spiders smaller than mimetic insects? ---\n"))
  # cat(sprintf('All body lengths range %s\n', paste(range(dorsalData$species$Coe$fac$bodylength), collapse = ' - ')))
  # #.meanDiff(dd, 'Dorsal', col = 'bodylength')
  # print(TukeyHSD(aov(bodylength ~ type, data = dorsalData$species$Coe$fac)))

  
  cat(sprintf("--- Are spiders less accurate mimics than insects? ---\n"))
  .meanDiff(dd, 'Dorsal')
  .varDiff(dd, 'Dorsal')
  .meanDiff(dl, 'Lateral')
  .varDiff(dl, 'Lateral')

  if (reportWorseThanModels) {
    # Are mimics worse than models?
    cat(sprintf("--- Are mimics less ant-like than ants? ---\n"))
    .meanDiff(dd, 'Dorsal spiders', rightType = MTP_MODEL)
    .meanDiff(dl, 'Lateral spiders', rightType = MTP_MODEL)
    .meanDiff(dd, 'Dorsal insects', leftTypes = MTP_INSECT, rightType = MTP_MODEL)
    .meanDiff(dl, 'Lateral insects', leftTypes = MTP_INSECT, rightType = MTP_MODEL)
    .meanDiff(dd, 'All dorsal mimics', leftTypes = MTP_MIMICS, rightType = MTP_MODEL)
    .meanDiff(dl, 'All lateral mimics', leftTypes = MTP_MIMICS, rightType = MTP_MODEL)
  }
  
  
  dd <- dorsalData$photo$fac[dorsalData$photo$fac$type %in% MTP_MIMICS,]$description
  ld <- lateralData$photo$fac[lateralData$photo$fac$type %in% MTP_MIMICS,]$description
  nid <- sum(grepl('interpolat', dd, ignore.case = TRUE))
  nil <- sum(grepl('interpolat', ld, ignore.case = TRUE))
  cat(sprintf("Proportion of interpolated mimics outlines: dorsal %g%%, lateral %g%%\n", round(100 * nid / length(dd), 0), round(100 * nil / length(ld), 0)))
}

# Adds mimetic accuracy based on LDA for photos, individuals and species to the specified dataset 
incorporateLDAAccuracy <- function(data) {
  
  # Photos
  data$photo$fac$accuracy <- MpLDACalcAccuracy(data$photo)$accuracy
  # Individuals
  data$individual$Coe$fac$accuracy <- MpLDACalcAccuracy(data$individual$Coe)$accuracy
  # Species
  data$species$Coe$fac$accuracy <- MpLDACalcAccuracy(data$species$Coe)$accuracy
  
  data
}

reportPhotosPerIndividual <- function(dorsalData, lateralData) {
  .tab <- function(data) {
    angle <- sub('Data', '', substitute(data))
    d <- data$photo$fac[,c('type', 'imageableid')]
    # Count number of photos per individual
    da <- aggregate(d$imageableid, by = list(d$type, d$imageableid), FUN = length)
    # Calc mean number of photos for each type
    cbind(angle, aggregate(da[,3], list(da[,1]), mean))
  }
  ld <- .tab(lateralData)
  names(ld) <- c('angle', 'type', 'photos.per.individual')
  .ppi <- function(tp) signif(ld[ld$type == tp,]$photos.per.individual, 2)
  cat(sprintf("Photos per individual: lateral spider mimics %g, lateral insect mimics %g\n", .ppi('mimetic spider'), .ppi('mimetic insect')))
}

reportSizeDiff <- function(dorsalData) {
  
}


################################################################################


# wantAngle <- 'dorsal'
# force <- TRUE
# 
# if (force || !exists("pca") || angle != wantAngle) {
#   angle <- wantAngle
#   query <- sprintf("ptype=Outline&imageable_type=Specimen&view_angle=%s", angle)
#   x = LoadPhotos(query)
#   cat(sprintf("Loaded %d %s photos\n", length(x), angle))
#   startTime <- proc.time()
#   
#   # xx <- subset(x, 1:6)
#   # p <- fgProcrustes(xx)
#   
#   p <- fgProcrustes(x)
#   cat(sprintf("Procrustes alignment took %g secs per shape\n", (proc.time() - startTime)[3] / length(p$coo)))
#   minCoords <- min(sapply(x$coo, length)) / 2
#   f <- efourier(p, norm = T, nb.h = minCoords %/% 2)
#   #plotPca(PCA(f), labels = "label")
#   PlotToPng(file.path(MORPHO_DIR, paste0('morpho-', angle, '-individual variation.png')),
#             function () VisualiseVariation(PCA(f), 'individualId'), width = 1600, height = 800)
#   # Average multiple outlines to individual outlines
#   individuals <- mshapes(f, 'individualId')
#   cat(sprintf("Loaded outlines of %d individuals\n", length(individuals$Coe)))
#   #plotPca(PCA(individuals$Coe), labels = "label")
#   # Now average to species
#   ms <- mshapes(individuals$Coe, 'speciesId')
#   cat(sprintf("Loaded %d species/outline type/form\n", length(ms$Coe)))
#   # PCA
#   pca <- PCA(ms$Coe)
#   ShowTime("Calculation time", startTime)
#   
#   
#   # NOTES on mshapes. fac must be a factor
#   
#   #X11()
#   hiFn <- function(colours) {
#     #hiIdx <- which(pca$fac$scientificName == "Myrmarachne sp6")
#     hiIdx <- which(pca$fac$imageableid == 508)
#     HighlightPcaPoints(pca, hiIdx, labels = "scientificName", radius = .01)
#   }
#   plotPca(pca, labels = "scientificName", extraPlottingFn = hiFn)
# 
#   # EXP
#   averageTypes <- mshapes(ms$Coe, 'type')
#   plotPca(PCA(averageTypes$Coe), labels = 'type')
# 
#   ExportVariation(individuals, ms, averageTypes, angle, MORPHO_DIR)
# 
#   PlotToPng(file.path(MORPHO_DIR, paste0('mimetic variation ', angle, '.png')), function () {
#     plotVariation(ms, averageTypes)
#   })
#   
#   ###
#   # Export CSV of mimetic imperfection for each photo
#   UpdateImperfectionCSV(MORPHO_DIR, averageTypes, f)
#   
#   ### 
#   # Simplified plots for poster
#   ResultsPlotForPoster(angle, ms, averageTypes, pca)
#   ###
#   
# 
#   w <- 1600
#   #h <- round(w * 460 / 800)
#   h <- round(w * 500 / 800)
#   # Plot bare axes to help explain meaning of morphometrics plot
#   # PlotToPng(file.path(MORPHO_DIR, paste0('morpho-', angle, '-bare.png')),
#   #           function () { plotPca(pca, bareAxes = TRUE) },
#   #           width = w, height = h)
#   PlotToPng(file.path(MORPHO_DIR, paste0('morpho-', angle, '.png')),
#             function () { plotPca(pca, labels = 'scientificName', cex = 2, labels.col = NULL, labels.cex = 1, legend = TRUE) },
#             width = w, height = h)
#   
#   #PauseUntilWindowIsClosed()
# 
#   #OverlayMeanShapes(f)
# 
#   # Show signifance of the first 9 harmonics in one shape 
#   #shapeIndex <- as.numeric(which(p$fac$id == 620)) # Which shape to plot
#   #PlotHarmonics(p, shapeIndex)
#   
#   PlotToPng(file.path(MORPHO_DIR, paste0('9th harmonic-', angle, '.png')),
#             function() {
#               indices <- order(p$fac$imageableid)
#               Plot9thHarmonic(p, sprintf('Specimen %d, Photo %d', p$fac$imageableid, p$fac$id), indices)
#             },
#             width = 1600, height = 1600
#   )
# }
# 
# 
# ##################################################################################
# # Morphological distances to test reproducability
# 
# pcaf <- PCA(f)
# PlotDistanceVariations(pcaf, pcaf$fac$family == 'Formicidae', ', Ants')
# PlotDistanceVariations(pcaf, pcaf$fac$order == 'Araneae', ', All Spiders')
# PlotDistanceVariations(pcaf, 1:nrow(pcaf$x), ', All')
# 
# # Some hypotheses can be tested with measures of phenotypic variation, which is what this represents
# PlotDistanceVariations(pcaf, pcaf$fac$type == MTP_INSECT, ', Insect mimics')
# PlotDistanceVariations(pcaf, pcaf$fac$type == MTP_SPIDER, ', Spider mimics')
# 
# # Random points to test NULL hypothesis
# pngWidth <- 800
# pngHeight <- 800
# PlotToPng(file.path(MORPHO_DIR, 'morpho-distances-random.png'),
#           function () {
#             numObs <- 5000
#             numSpecies <- 30
#             speciesPerGenus <- 2
#             photosPerIndividual <- 2
#             points <- data.frame(x = rnorm(numObs), y = rnorm(numObs))
#             fac <- GenerateRandomFacData(numObs, speciesPerGenus, numSpecies, photosPerIndividual)
#             PlotMorphoDistance(points, fac, 'Morphometric Distances, randomly generated points')
#           },
#           width = pngWidth, height = pngHeight)
# 

PlotsForFinalTalk <- function(dorsalData, lateralData) {
  TALK_DIR <- 'C:/Users/jim_m/Google Drive/Uni/Classes/MRes Thesis/Final talk'
  
  # Plots
  w <- 900
  h <- round(w * 500 / 800)
  
  # Plot densities function
  .pd <- function(densities, types, main = NA) {
    JPlotDensities(densities, col = typeToCol(types), lwd = 4, lty = typeToLty(types), main = main, cex.main = 1.5, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    MLabelAccuracyAxis(cex = 1.5)
    title(ylab = 'Probability density', line = 1, cex.lab = 1.5)
  }
  
  # Plot prediction function
  .prs <- function(sd1 = 1, sd2 = .6, doLegend = TRUE) {
    set.seed(4)
    np = 300
    spider <- density(rnorm(n = np, sd = sd1, mean = 0))
    insect <- density(rnorm(n = np, sd = sd2, mean = 1))
    .pd(list(spider, insect), MTP_MIMICS, main = 'Prediction')
    if (doLegend)
      legend("topright", c('Spiders', 'Insects'), lwd = 4, cex = 1.5, lty = typeToLty(MTP_MIMICS), col = typeToCol(MTP_MIMICS), inset = c(.01, .01))
  }
  
  # Classic prediction for high accuracy
  PlotToPng(file.path(TALK_DIR, 'Prediction classic.png'), function () {
    par(mar = c(3.5, 2.6, 2, 0.1))
    set.seed(4)
    np = 300
    bad <- density(rnorm(n = np, sd = 1, mean = 0))
    good <- density(rnorm(n = np, sd = .6, mean = 1))
    .pd(list(bad, good), MTP_MIMICS, main = 'Expectation')
    # Rectangle showing predation pressure
    x1 <- -3
    y1 <- 0
    x2 <- 0
    y2 <- .45
    polygon(c(x1, x2, x2, x1), c(y1, y1, y2, y2), col = JTransparentColour(MyPallete$lightGrey, 80))
    text((x1 + x2) / 2, y2, labels = 'Predation pressure', pos = 3)
    legend("topleft", c('Before selection', 'After selection'), lwd = 4, cex = 1.5, lty = typeToLty(MTP_MIMICS), col = typeToCol(MTP_MIMICS), inset = c(.01, .01))
  }, width = 450, height = 450)
  
  
  # Prediction of constraints hypothesis
  PlotToPng(file.path(TALK_DIR, 'Prediction constraints.png'), function () {
    par(mar = c(3.5, 2.6, 2, 0.1))
    .prs(.6, 1)
  }, width = 450, height = 450)
  

  # Explanation of morphospace
  .draw_shape <- function(shape, xy, scale, col) {
    shape <- shape * scale
    shape[,1] <- shape[,1] + xy[1]
    shape[,2] <- shape[,2] + xy[2]
    polygon(shape, col = col)
  }
  
  PlotToPng(file.path(TALK_DIR, 'morpho-explain.png'), function () {
    pca <- PCA(dorsalData$species$Coe)
    .ds <- function(colours) {
      .dsp <- function(species, offset = 2, pos = 1) {
        data <- dorsalData$species
        idx <- which(data$Coe$fac$scientificName == species)[1]
        xy <- pca$x[idx, c(1, 2)]
        .draw_shape(data$shp[[idx]], xy, .04, typeToCol(data$Coe$fac$type[idx]))
        label <- data$Coe$fac$genus[idx]
        label <- sub('[0-9]+', '', label)
        text(xy[1], xy[2], labels = label, pos = pos, offset = offset, cex = 1.5)
      }
      .dsp('Sphinctomyrmex sp1', 1.5)
      .dsp('Damoetas sp3', offset = 2.5)
      #.dsp('Riptortus abdominalis')
      .dsp('Vespoidea6 sp1', pos = 3)
      .dsp('Servaea incana', offset = 2.5)
      mtext('Morphospace', side = 3, adj = 0.05, cex = 2)
    }
    # Unfortunately, momocs pca plot doesn't allow symbol fill colours to be specified
    cex <- 1.6
    plotMorphoPca(pca, points = FALSE, cex = cex, labels = FALSE, labels.col = NULL, ellipses = FALSE, zoom = 1.4, legend = FALSE, extraPlottingFn = .ds)
    xy <- pca$x[, c(1, 2)]
    #text(xy, labels = pca$fac[, 'scientificName'], pos = 1, cex = 1.2, col = typeToCol(pca$fac$type))
    #points(xy, pch = typeToPch(pca$fac$type), col = typeToCol(pca$fac$type), bg = typeToCol(pca$fac$type), cex = cex)
    
  }, width = w, height = h)
  
  ### Results
  
  # Plot of morpho-space
  PlotToPng(file.path(TALK_DIR, 'morpho.png'), function () {
    pca <- PCA(dorsalData$species$Coe)
    # Unfortunately, momocs pca plot doesn't allow symbol fill colours to be specified
    cex <- 1.6
    plotMorphoPca(pca, points = FALSE, cex = cex, labels.cex = 1.2, labels = FALSE, ellipsesax = FALSE, zoom = 1.4, legend = TRUE)
    xy <- pca$x[, c(1, 2)]
    #text(xy, labels = pca$fac[, 'scientificName'], pos = 1, cex = 1.2, col = typeToCol(pca$fac$type))
    points(xy, pch = typeToPch(pca$fac$type), col = typeToCol(pca$fac$type), bg = typeToCol(pca$fac$type), cex = cex)
    
  }, width = w, height = h)
  
  # Plot of imperfection by mimic type
  dd <- dorsalData$species$Coe$fac
  
  .pad <- function(data, main = NA) {
    dd <- data$species$Coe$fac
    dd <- dd[dd$type %in% MTP_MIMICS,]
    types <- unique(dd$type)
    densities <- lapply(types, function(tp) density(dd$accuracy[dd$type == tp], na.rm = TRUE))
    names(densities) <- types
    .pd(densities, types, main)
  }
  
  PlotToPng(file.path(TALK_DIR, 'mimetic variation dorsal.png'), function () {
    JSetParsTemporarily(mar = c(3.5, 2.5, 3.5, .1))
    .pad(dorsalData, 'Dorsal')
    legend("topright", c('Spiders', 'Insects'), lwd = 4, cex = 1.5, lty = typeToLty(MTP_MIMICS), col = typeToCol(MTP_MIMICS), inset = c(.01, .01))
  }, width = 450, height = 450)
  PlotToPng(file.path(TALK_DIR, 'mimetic variation lateral.png'), function () {
    JSetParsTemporarily(mar = c(3.5, 2.5, 3.5, .1))
    .pad(lateralData, 'Lateral')
    legend("topright", c('Spiders', 'Insects'), lwd = 4, cex = 1.5, lty = typeToLty(MTP_MIMICS), col = typeToCol(MTP_MIMICS), inset = c(.01, .01))
  }, width = 450, height = 450)
  
  ###
  # Relaxed selection on spiders prediction + lateral result (again)

  PlotToPng(file.path(TALK_DIR, 'shapes relaxed selection.png'), function () {
    JSetParsTemporarily(mfrow = c(1, 2), mar = c(3.5, 2.5, 2.5, .1))
    .prs()
    .pad(lateralData, main = 'Lateral Result')
    legend("topright", c('Spiders', 'Insects'), lwd = 4, cex = 1.5, lty = typeToLty(MTP_MIMICS), col = typeToCol(MTP_MIMICS), inset = c(.01, .01))
  }, width = w, height = w / 2)
  
  ###
  # Relaxed selection on smaller mimics results

  PlotToPng(file.path(TALK_DIR, 'size relaxed selection.png'), function () {
    JSetParsTemporarily(mar = c(4, 3.5, .1, .1))
    PlotAccuracyVsSize(dorsalData$species$Coe$fac, cex = 1.8, lwd = 4)
  }, width = 600, height = 600)
  
  ###
  # Map of collecting sites
  specs <- DbQuerySpecimens()
  PlotToPng(file.path(TALK_DIR, 'map.png'), function () {
    #specs <- dorsalData$individual$Coe$fac
    specs <- specs[grepl('jim|louis', specs$recordedBy, ignore.case = TRUE),]
    GMapSpecimenSites(specs, xlimFrac = 1, ylimFrac = 1)
  })
}

CombinedMorphoPlots <- function(dorsalData, lateralData) {
  # Plot of PCA morpho-space  
  .ppca <- function(data, label) {
    pca <- PCA(data$species$Coe)
    plotMorphoPca(pca, labels.cex = 1.2, labels = FALSE, ellipsesax = FALSE, zoom = 1.4, legend = TRUE, legendInsetFn = JInset, legendYIntersp = .7)
    # I don't know what's going on with the coordinate system here, 
    # line = 2 should be outside the graph
    mtext(label, line = 2, adj = 0.02)
    par(xpd = NA)
    var <- pca$sdev^2
    var <- signif(100 * var/sum(var), 3)
    mtext(sprintf("PC1 (%g%%)", var[1]), side = 1, adj = 0.95, padj = 5)
    mtext(sprintf("PC2 (%g%%)", var[2]), side = 2, adj = 0.1, padj = -4)
    cat(sprintf("%s: 1st 2 PCA components cover %g%% of total variance\n", label, round(100 * scree(pca, 1:2)[2,3])))
  }

  PlotToPng(file.path(MORPHO_DIR, paste0('morpho-pca.png')), function () {
    JSetParsTemporarily(mfrow = c(2, 1))
    .ppca(dorsalData, "a) Dorsal shapes")
    .ppca(lateralData, "b) Lateral shapes")
  }, width = 500, height = 700)

  ###
  # Plot of distribution of accuracy
  
  doIndividuals <- FALSE
  
  .pacc <- function(data, label, accuracyCol = 'accuracy') {
    dd <- data$species$Coe$fac
    if (doIndividuals)
      dd <- data$individual$Coe$fac
    subset <- dd$type == MTP_SPIDER | dd$type == MTP_INSECT
    meanAnt <- mean(dd[dd$type == MTP_MODEL,accuracyCol])

    # Vertical line at mean value for ants
    .annotate <- function(densities) { 
      segments(meanAnt, 0, meanAnt, 0.8 * par()$usr[4], col = MyPallete$lightGrey, lwd = 2)
      oldPar <- par(srt = -90)
      text(meanAnt, 0.8 * par()$usr[4], 'Mean ant shape', adj = c(0, 1.4))
      par(oldPar)
    }
    MpPlotAccuracyDensities(dd$type, dd[,accuracyCol], subset = subset, showSpeciesCount = FALSE, extrasFn = .annotate, includeInX = meanAnt, mar = NULL)
    mtext(label, adj = 0.01, line = -1.5)
  }
  
  nm <- ifelse(doIndividuals, 'mimetic-accuracy-ind.png', 'mimetic-accuracy.png')
  PlotToPng(file.path(MORPHO_DIR, paste0(nm)), function () {
    JSetParsTemporarily(mar = c(3.5, 2.5, .2, .2), mfrow = c(2, 1))
    #ac <- 'gvAccuracy'
    ac <- 'accuracy'
    .pacc(dorsalData, "a) Dorsal shapes", accuracyCol = ac)
    .pacc(lateralData, "b) Lateral shapes", accuracyCol = ac)
  }, width = 450, height = 700)
}

MorphoAnalysisBigPlots <- function(data, accuracyCol = 'accuracy') {
  angle <- tolower(data$photo$fac$angle[1])

  # Plots
  w <- 1600
  h <- round(w * 500 / 800)

  # Plot of morpho-space
  PlotToPng(file.path(MORPHO_DIR, paste0('morpho-', angle, '.png')), function () {
    plotMorphoPca(PCA(data$species$Coe), cex = 1, labels.cex = 1.2, labels = 'scientificName', ellipsesax = FALSE, zoom = 1.4, legend = TRUE)
  }, width = w, height = h)

  # Plot of imperfection by mimic type
  dd <- data$species$Coe$fac

  PlotToPng(file.path(MORPHO_DIR, paste0('mimetic variation ', angle, '.png')), function () {
    JSetParsTemporarily(mar = c(4.2, 4.2, 3, 1))
    MpPlotAccuracyDensities(dd$type, dd[,accuracyCol])
  }, width = 600, height = 450)

  # Test if selection is relaxed based on size
  #print(summary(PlotAccuracyVsSize(data$type$Coe, data$individual$Coe)))
}

doAnalysis <- function() {
  checkUpToDate = ifelse(JCLIOption('f'), 'force', TRUE)
  dorsalData <- MGetData('dorsal', checkUpToDate)
  lateralData <- MGetData('lateral', checkUpToDate)
  
  # Incorporate LDA accuracy
  dorsalData <- incorporateLDAAccuracy(dorsalData)
  lateralData <- incorporateLDAAccuracy(lateralData)
  # Incorporate GV accuracy
  dorsalData <- GVIncorporateAccuracy(dorsalData)
  lateralData <- GVIncorporateAccuracy(lateralData)

  # Plots for MRes final seminar
  #PlotsForFinalTalk(dorsalData, lateralData)
  
  MExportVariation(dorsalData, lateralData, MORPHO_DIR)
  # Appendix - list of specimens in study, where they came from 
  MExportAppendixIndividuals(dorsalData, lateralData, MORPHO_DIR, 'Table S2')
  MExportAppendixSpecies(dorsalData, lateralData, MORPHO_DIR, 'Table S3')
  # Correlation of lateral to dorsal accuracy, morphometrics and machine learning
  CombinedDorsalLateralPlot(MORPHO_DIR)
  
  MorphoAnalysisBigPlots(dorsalData, 'accuracy')
  MorphoAnalysisBigPlots(lateralData, 'accuracy')
  CombinedMorphoPlots(dorsalData, lateralData)
  cat("\n--- Relaxed selection ---\n")
  CombinedRelaxedSelectionPlots(dorsalData, lateralData)
  cat("\n--- Machine learning ---\n")
  CombinedMorphoGvCorrelation(dorsalData, lateralData)
  
  cat("\n--- Constraints ---\n")
  reportMorphoStatistics(dorsalData, lateralData)
  #reportMorphoStatistics(dorsalData, lateralData, 'gvAccuracy', reportWorseThanModels = TRUE)
  
  # Don't think I'll report this
  #reportPhotosPerIndividual(dorsalData, lateralData)
}

if (!JIsRStudio()) doAnalysis()

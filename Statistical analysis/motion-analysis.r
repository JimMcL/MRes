#!Rscript
#
# Arthropod trajectory analysis.
# Analyses CSV files of animal trajectories, incorporating individual and 
# species information from the sampleIt database.

#library(rgl)
library(corrplot)
library(MASS)
suppressWarnings(suppressMessages(library(Momocs)))

LoadFns('motion')
LoadFns('download')
LoadFns('rediscretization')
LoadFns('straightness')
LoadFns('motion-stats')
LoadFns('graphics')
LoadFns('motion-morpho')
LoadFns('morpho-distance')

PROP_EXT <- '.prop'

# List of column names for all statistics
STATS_NAMES <- c("maximum.speed.bodylength", "mean.speed.bodylength", "speed.CV",
                 "moving.duration.mean", "moving.duration.CV", 
                 "stopped.duration.mean", "stopped.duration.CV", 
                 "proportion.time.moving", 
                 "straightness", "sinuosity", "E.max", 
                 "first.minimum.C", "first.minimum.delta.s")



#file <- 'file:///C:/Jim/uni/Classes/Thesis/Behaviour/Motion-videos/CIMG9107Trim.csv'
# Well-behaved myrmarachne
#file <- 'file:///C:/Jim/uni/Classes/Thesis/Behaviour/Motion-videos/CIMG9106Trim.csv'
# Stop/start jumping spider
#file <- 'file:///C:/Jim/uni/Classes/Thesis/Behaviour/Motion-videos/708.csv'
#file <- 'file:///C:/Jim/uni/Classes/Thesis/Behaviour/Motion-videos/947-900x600.csv'
#file <- 'file:///C:/Jim/uni/Classes/Thesis/Behaviour/Motion-videos/661.csv'
#file <- 'C:/Jim/tmp/junk.csv'

# trackVelocities <- function(track, lag = 1, clampNoise = TRUE) {
#   # Differentiate displacement to get speed
#   diffs <- diff(track$polar, lag = lag)
#   displacement <- Mod(diffs)
#   # Calculate change in angle from one segment to the next.
#   # There is 1 fewer angular changes than displacements
#   angularChange <- diff(Arg(diffs))
# 
#   timeIntervals <- diff(track$Time, lag = lag)
#   speed <- displacement / timeIntervals
#   angularVelocity <- angularChange / timeIntervals[2:length(timeIntervals)]
#   
#   if (clampNoise) {
#     qnt <- quantile(angularVelocity, probs=c(.15, .85))
#     angularVelocity[angularVelocity < qnt[1]] <- NA
#     angularVelocity[angularVelocity > qnt[2]] <- NA
#   }
#   list(speed = speed, angularVelocity = angularVelocity)
# }

layoutPlotsInSqr <- function(np, rowToColfactor = 1) {
  nrow <- floor(sqrt(np))
  ncol <- ceiling(np / nrow)
  par(mfrow=c(nrow, ncol))
}

testTrack <- function() {
  # Displacement:   -, 1, 1, 1, 0, sqrt(2), sqrt(2), 0,  1,  1
  # Direction:      -, 0, 0, 0, 0,      45,      45, 0,   90, 90
  # Change in dirn: -, -, 0, 0, 0,      45,       0, -45, 90,  0
  track <- data.frame(
    x = c(0, 1, 2, 3, 3, 4, 5, 5, 5, 5),
    y = c(0, 0, 0, 0, 0, 1, 2, 2, 3, 4)) 
  track$polar <- complex(real = track$x, imaginary = track$y)
  track$Frame <- 1:nrow(track)
  # 2 fps
  track$Time <- 1:nrow(track) / 2
  
  v <- trackVelocities(track, clampNoise = FALSE)
  stopifnot(all.equal(v[[1]], c(2, 2, 2, 0, 2 * sqrt(2), 2 * sqrt(2), 0, 2, 2)))
  stopifnot(all.equal(v[[2]], 2 * c(0, 0, 0, pi / 4, 0, -pi / 4, pi / 2, 0)))
}

plotTrack <- function(track, ...){
  plot(y ~ x, data = track, type = 'l', asp = 1, ...)
  points(track$x[1], track$y[1], pch = 16, cex = .8)
}

# Reads all of the specified tracks
readTracks <- function(trackVideos) {
  mimicTypes <- mimicType(trackVideos)
  
  .rowToTrack <- function(rowIndex) {
    MReadTrack(trackVideos[rowIndex,]$trackFile, mimicTypes[rowIndex], TRUE)
  }
  lapply(1:nrow(trackVideos), .rowToTrack)
}

# Reads all tracks in a directory, assuming they are named <number>.csv
dirTracks <- function(dir, specimens) {
  #files <- list.files(dir, '[[:digit:]]+.csv', full.names = TRUE)
  files <- list.files(dir, '.csv', full.names = TRUE)
  lapply(files, readTrack, specimens)
}

trackToCol <- function(track) {
  typeToCol(attr(track, 'type'))
}

plotTracks <- function(tracks) {
  oldPars <- par(no.readonly = T)
  on.exit(par(oldPars))
  par(mar=c(0.2, 0.2, 3, 0.2))
  
  xlim <- range(sapply(tracks, function(t) range(t$x)))
  ylim <- range(sapply(tracks, function(t) range(t$y)))
  plot(NULL, main = 'Arthropod tracks', 
       xlim = xlim, ylim = ylim, 
       asp = 1, 
       xaxt = 'n', xlab = '', yaxt = 'n', ylab = '')
  box(lwd = 3)
  
  .plotTrack <- function(t) {
    col <- trackToCol(t)
    points(t$x[1], t$y[1], col = col, pch = 16, cex = .6)
    lines(y ~ x, data = t, col = col)
  }
  
  sapply(tracks, .plotTrack)
  leg <- levels(attr(tracks[[1]], 'type'))
  if (length(leg) > 0)
    legend("topleft", leg, lwd = 2, col = c(MyPallete$darkGreen,MyPallete$darkRed, MyPallete$darkBlue, MyPallete$darkPurple), inset = c(0.01, 0.01))
}

plotLineFn <- function(tracks, fn) {
  l <- list()
  for(i in 1:length(tracks)) {
    l[[i]] <- fn(tracks[[i]])
  }
  xlim <- c(0, length(l[[1]]))
  ylim <- range(sapply(l, range))
  plot(NULL, xlim = xlim, ylim = ylim)
  for(i in 1:length(tracks)) {
    track <- tracks[[i]]
    lines(l[[i]], col = trackToCol(track), lwd = 2)
  }
  legend("topleft", levels(attr(tracks[[1]], 'type')), lwd = 2, col = cols, inset = c(0.01, 0.01))
}

plotAllTracksInFile <- function (file) {
  tracks <- read.csv(file, comment.char = '#')
  xlim <- range(tracks$x)
  ylim <- range(tracks$y)
  plot(NULL, xlim = xlim, ylim = ylim, asp = 1)
  
  polar <- complex(real = tracks$x, imaginary = tracks$y)
  
  # Determine the longest track
  tb <- table(tracks$TrackId)
  btid <- names(tb)[which(tb == max(tb))[1]]
  
  .plotTrack <- function(tid, col = 'black', lwd = 1) {
    pts <- tracks[tracks$TrackId == tid,]
    pp <- polar[tracks$TrackId == tid]
    cat(sprintf("Track id %d, length %g, distance %g (%d points)\n", tid, MTrackLength(pp), MTrackLength(pp[c(1,length(pp))]), nrow(pts)))
    # points(pts$x[1], pts$y[1], pch = 16, col = col)
    # points(pts$x[length(pts)], pts$y[length(pts)], pch = 4, col = col)
    # lines(pts$x, pts$y, col = col, lwd = lwd)
    plot3d(pts$x, pts$y, pts$Time, type = 'l', col = col, lwd = lwd)
  }
  
  for (tid in unique(tracks$TrackId)) { 
    pts <- tracks[tracks$TrackId == tid,]
    pp <- polar[tracks$TrackId == tid]
    cat(sprintf("Track id %d, length %g, distance %g (%d points)\n", tid, MTrackLength(pp), MTrackLength(pp[c(1,length(pp))]), nrow(pts)))
    if (nrow(pts) > 1) {
      points(pts$x[1], pts$y[1], pch = 16, col = tid + 1)
      lines(pts$x, pts$y, col = tid + 1)
    } else if (nrow(pts) == 1) {
        points(pts$x, pts$y)
    }
  }
}

# Checks for the existence of a CSV file for every .prop file.
# Lists any that are missing
ListMissingCsvs <- function(dir) {
  props <- list.files(dir, PROP_EXT, full.names = TRUE)
  for(file in props) {
    csv <- paste0(tools::file_path_sans_ext(file), '.csv')
    if (!file.exists(csv)) {
      # Emit the command used to generate the CVS file, with feedback
      cat(sprintf("/cygdrive/c/Jim/uni/apps/VideoAnalysis/run.bat --defaults %s -g -s 10 --display-track -v\n", basename(file)))
    }
  }
}

# Used to sanity check a track file's contents.
# Displays all tracks, then the selected "best" track, 
# then speed and acceleration for the best track
debugFile <- function(file) {
  plotAllTracksInFile(file)

  invisible(readline(prompt="Press [enter] to continue"))
  
  track <- MReadTrack(file, 'mimic', TRUE)
  cat(sprintf("Stopped for %d%% of the time", stoppedPercent(track)))
  plotTrack(track)

  invisible(readline(prompt="Press [enter] to continue"))
  
  # Convert cartesian to polar coordinates
  polar <- track$polar
  # Differentiate displacement to get speed
  displacement <- diff(polar)
  t1 <- track$Time[2:nrow(track)]
  sm <- sm.spline(t1, Mod(displacement))
  smr <- sm.spline(t1, Arg(displacement))
  speed <- sm$y[2:length(sm$y)] / diff(sm$x)
  angularVelocity <- smr$y[2:length(smr$y)] / diff(smr$x)
  t1 <- sm$x[2:length(sm$x)]
  acc <- diff(speed) / diff(t1)
  t2 <- t1[2:length(t1)]
  plot(acc ~ t2, pch = 16, cex = .5, col = 'red', ylim = c(-100, 250))
  points(speed ~ t1, type = 'l', pch = 16, col = 'blue', cex = .4)
  abline(h = 0, col = MyPallete$midGrey, lty = 2)
}

plotDerivSpeed <- function(deriv, bodyLength, title = NULL, xlim = NULL, ylim = NULL) {
  plot(deriv$speed / bodyLength ~ deriv$speedTimes, 
       type = 'l', 
       xlim = xlim, ylim = ylim,
       main = title, 
       xlab = 'Time (sec)', ylab = 'Speed (body length / sec)')
  abline(h = STOPPED_SPEED / bodyLength, col = MyPallete$midRed)
}

plotSpeed <- function(derivs, labels, xlim = c(0, 40), ylim = c(0, 80)) {
  oldPars <- par(no.readonly = T)
  on.exit(par(oldPars))
  
  nplots <- length(derivs)
  # Want more rows than columns
  layoutPlotsInSqr(nplots, 3)

  if(is.null(xlim) || is.na(xlim))
    xlim <- range(sapply(derivs, function(d) range(d$speedTimes)))
  if(is.null(ylim) || is.na(ylim))
    ylim <- range(sapply(derivs, function(d) range(d$speed)))
  
  for (i in 1:length(derivs)) {
    d <- derivs[[i]]
    plot(d$speed ~ d$speedTimes, type = 'l', 
         xlim = xlim, ylim = ylim,
         main = labels[i], 
         xlab = 'Time', ylab = 'Velocity')
    abline(h = STOPPED_SPEED, col = MyPallete$midRed)
  }
}

plotAngles <- function(tracks, labels) {
  oldPars <- par(no.readonly = T)
  on.exit(par(oldPars))
  
  dens <- lapply(tracks, function(t) density(MTrackAngles(MRediscretizeTrack(t, 1))))
  
  layoutPlotsInSqr(length(dens), 3)

  xlim <- range(sapply(dens, function(d) range(d$x)))
  ylim <- range(sapply(dens, function(d) range(d$y)))
  
  for (i in 1:length(dens)) {
    d <- dens[[i]]
    plot(d, 
         xlim = xlim, #ylim = ylim,
         main = labels[i], 
         xlab = 'Angle', ylab = 'Density')
    #abline(v = 0, col = MyPallete$lightGrey)
  }
}

plotSpeedForType <- function(derivs, videos, type) {
  idx <- which(videos$type == type)
  idx <- idx[order(videos[idx,]$scientificName)]
  plotSpeed(derivs[idx], sprintf("%s, video %d", videos[idx,]$scientificName, videos[idx,]$id))
}

plotAnglesForType <- function(videos, type) {
  idx <- which(videos$type == type)
  idx <- idx[order(videos[idx,]$scientificName)]
  plotAngles(videos[idx,]$track, sprintf("%s, video %d", videos[idx,]$scientificName, videos[idx,]$id))
}

exportSpeedPlotsForTypes <- function(derivs, videos, dir) {
  .plotType <- function() { plotSpeedForType(derivs, videos, type) }
  for (type in unique(videos$type)) {
    filename <- file.path(dir, sprintf("speed-%s.png", type))
    PlotToPng(filename, .plotType, width = 1800, height = 1800)
  }
}

exportAnglePlotsForTypes <- function(videos, dir) {
  .plotType <- function() { plotAnglesForType(videos, type) }
  for (type in unique(videos$type)) {
    filename <- file.path(dir, sprintf("angle-%s.png", type))
    PlotToPng(filename, .plotType, width = 1800, height = 1800)
  }
}

plotTypeTracks <- function(videos) {
  # Order types appropriately
  types <- AsPlottableMimicTypeFactor(videos$type)
  ntypes <- length(unique(types))

  oldPars <- par(no.readonly = T)
  on.exit(par(oldPars))
  par(mar=c(0.2, 0.2, 0.2, 0.2))
  par(mfrow = c(1, ntypes))
  
  tracks <- videos$track
  xlim <- range(sapply(tracks, function(t) range(t$x)))
  ylim <- range(sapply(tracks, function(t) range(t$y)))
  
  .plotType <- function(type) {
    idx <- which(types == type)
    col <- typeToCol(type)

    plot(NULL,
         xlim = xlim, ylim = ylim, 
         asp = 1, 
         xaxt = 'n', xlab = '', yaxt = 'n', ylab = '')
    sapply(tracks[idx], function(t) {
      #points(t$x[1], t$y[1], col = col, pch = 16, cex = .6)
      lines(y ~ x, data = t, col = col, cex = 1.1)
    })
    points(0, 0, pch = 16)
    mtext(sprintf("%s (n = %d)", capsentence(paste0(type, 's')), length(idx)), 1, padj = -1.4)
  }
  invisible(sapply(rev(unique(types)), .plotType))
}

plotTurnDistribution <- function(track) {
  rediscretizationFactor <- 0.1
  segLen <- track$bodylength * rediscretizationFactor
  # Tracks are actually lists with a single element which is the track! Weird!
  rt <- MRediscretizeTrack(track$track[[1]], segLen)
  plot(density(cos(diff(Arg(tail(rt$displacement, -1))))), main = 'Turns', xlim = c(0.8, 1))
}

plotEverythingForTrack <- function(video, prefix = c('a', 'b', 'c'), specimenFmt = '%s) Specimen %d, %s') {
  # Plot the track itself
  tt <- video$track[[1]]
  plotTrack(tt, main = sprintf(specimenFmt, prefix[1], video$imageableid, video$scientificName), xlab = 'x (mm)', ylab = 'y (mm)')
  # CalculatePlot autocorrelation
  corr <- MGetDirnAutocorrelation(video)
  MPlotAutocorrelations(list(corr), video$type, title = paste0(prefix[2], ") Direction autocorrelation"))
  # Plot speed over time
  derivs <- MTrackDerivatives(tt)
  plotDerivSpeed(derivs, video$bodylength, paste0(prefix[3], ') Speed / time'), ylim = c(0, max(derivs$speed / video$bodylength)))
}

ExportEverythingForAll <- function(videos, dir) {
  .plotEverythingForType <- function(type) {
    # Subset tracks to those for type
    idx <- which(videos$type == type)
    idx <- idx[order(videos[idx,]$scientificName)]
    tt <- videos[idx,]

    JSetParsTemporarily(mfrow = c(nrow(tt), 3))
    
    for (i in 1:nrow(tt)) {
      plotEverythingForTrack(tt[i,])
    }
  }
  
  # Break into types
  for (type in unique(videos$type)) {
    PlotToPng(file.path(dir, sprintf("%s tracks.png", type)),
              function () .plotEverythingForType(type),
              width = 800, height = sum(videos$type == type) * 300)
  }
  
  # Export a single plot with an example of each type
  PlotToPng(file.path(dir, "EG tracks.png"), function () {
    JSetParsTemporarily(mfrow = c(3, 3), cex = 1, mar = c(4, 4, 2, 1), mgp = c(2.5, 1, 0))
    plotEverythingForTrack(videos[videos$imageableid == 968,], prefix = c('a', 'b', 'c'), specimenFmt = '%s) Ant Trajectory')
    plotEverythingForTrack(videos[videos$imageableid == 863,][1,], prefix = c('d', 'e', 'f'), specimenFmt = '%s) Mimic Trajectory')
    plotEverythingForTrack(videos[videos$imageableid == 582,], prefix = c('g', 'h', 'i'), specimenFmt = '%s) Non-mimic Trajectory')
  }, width = 800, height = 780, res = 80)
}

# Returns list(var, contrib), var is variance explained by first 2 components, contrib is matrix of relative contributions of each variable to the first 2 components
plotMotionPCA <- function(stats, videos, labels = TRUE, cex = 1, pt.cex = 1, fill = FALSE) {
  .cHull <- function(pts, lineColour, fillColour) {
    ch <- chull(pts)
    # Close the polygon
    ch <- c(ch, ch[1])
    
    lines(pts[ch,], col = lineColour)
    polygon(pts[ch,], border = NA, col = fillColour)
  }
  
  .typeContour <- function(pts, type, ...) {
    contour(.typeDensity(pts, type), add = TRUE, col = typeToCol(type), nlevels = 6, drawlabels = FALSE, ...)
  }
  
  .typeHeatmap <- function(pts, types) {
    image(.typeDensity(pts, types), add = TRUE)
  }
  

  # Run the PCA. It is very important to scale the values
  pca <- prcomp(~ ., data = removeStatsNAs(stats[,STATS_NAMES]), scale. = TRUE)
  
  # Plot
  pcaTypes <- AsPlottableMimicTypeFactor(videos$type)
  pts <- pca$x[,1:2]
  typeList <- unique(pcaTypes)
  xlim <- extendrange(pts[,1])
  ylim <- range(pts[,2])

  .typeCHull <- function(pts, type) {
    c <- typeToCol(type)
    .cHull(pts[pcaTypes == type,], lineColour = c, fillColour = JTransparentColour(c, 12))
  }
  .typeDensity <- function(pts, type) {
    kde2d(pts[pcaTypes == type,1], pts[pcaTypes == type,2], n = 100)
  }
  .confidenceEllipse <- function(pts, type, conf = .95) {
    col <- typeToCol(type)
    pts <- pts[pcaTypes == tp, ]
    if (nrow(pts) > 2) {
      ell <- conf_ell(x = pts, conf = conf)
      lines(coo_close(ell$ell), col = typeToCol(type))
      centroid <- coo_centpos(pts)
      points(centroid[1], centroid[2], pch = 3, cex = 0.3, col = col)
    }
  }
  
  plot(NULL, xlim = xlim, ylim = ylim, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
  # for(tp in typeList)
  #   .typeHeatmap(pts, tp)
  # for(tp in typeList)
  #   .typeCHull(pts, tp)
  # for(tp in typeList)
  #   .typeContour(pts, tp, xlim = xlim, ylim = ylim)
  bg <- NULL
  if (fill)
    bg <- typeToCol(pcaTypes)
  points(pts, pch = typeToPch(pcaTypes), col = typeToCol(pcaTypes), bg = bg, cex = pt.cex)
  if (labels)
    text(pts, videos[,]$scientificName, col = typeToCol(pcaTypes), cex = .75, pos = 1)
  var <- pca$sdev^2
  var <- signif(100 * var/sum(var), 3)
  title(xlab = sprintf('PC1 (%g%%)', var[1]), ylab = sprintf('PC2 (%g%%)', var[2]), line = 0.5 * cex, cex.lab = cex)

  for(tp in typeList)
    .confidenceEllipse(pts, tp)
  
  MLegend(typeList, sapply(typeList, function(tp) sum(pcaTypes == tp)), 
          pch = typeToPch(sort(AsPlottableMimicTypeFactor(typeList))),
          legendPos = 'bottomright', fill = fill, pt.cex = pt.cex, cex = cex)
  
  # Relative cotributions of each variable to PC1 and PC2
  absRot <- abs(pca$rotation)
  contrib <- sweep(absRot, 2, colSums(absRot), "/")[,c(1,2)]
  
  # Return variance explained by first 2 components
  pcai <-   summary(pca)$importance
  invisible(list(var = pcai[nrow(pcai), 2], contrib = contrib))
}

plotAllTracksByIndex <- function(videos) {
  # Sort from most to least ant-like
  videos <- videos[order(videos$accuracy),]
  
  oldPars <- par(no.readonly = T)
  on.exit(par(oldPars))
  layoutPlotsInSqr(nrow(videos))
  par(mar = c(4, 2, 2, 1))
  
  JApplyToRows(videos, function(video) {
    plotTrack(MTrackRotate(video$track[[1]], pi / 2), xlab = '', ylab = '')
    t <- sprintf('%s: %s, %g',capsentence(video$type), video$scientificName, video$imperfection)
    title(t)
  })

  invisible(NULL)  
}

# Include NAs as informative data for stats analysis.
# Removes NA from stats which result from a nonexistent first local minimum for direction autocorrelation:
# creates a new column called "noMin" which is true for rows with no local minimum,
# then replaces NA values with replacementValue
removeStatsNAs <- function(stats, replacementValue = 0) {
  # Assume there will only be NAs in first minimum columns
  stats$noMin <- as.numeric(is.na(stats$first.minimum.C))
  stats[is.na(stats)] <- 0
  stats
}

CalcMotionLDANoMimics <- function(stats, videos) {
  
  s <- removeStatsNAs(stats)
  s$type <- AsPlottableMimicTypeFactor(videos$type)

  # Train on models and non-mimics.
  trainingData <- s[which(s$type == MTP_MODEL | s$type == MTP_NON_MIMIC),]
  trainingData$type <- AsPlottableMimicTypeFactor(trainingData$type)  # Prevent warning about unused factor values
  # If we had some idea of the probability of a predator encountering a model vs a non-mimic, we could set the priors.
  # Since I don't know, I will make it 50/50. The default is to calculate priors from the number of samples, which is probably not meaningful
  r <- lda(type ~ ., data = trainingData, prior = c(1, 1) / 2)
  
  p <- predict(r, newdata = s)

  # Table of actual vs predicted
  tab <- table(s$type, p$class)
  names(dimnames(tab)) <- c('Actual', 'Predicted')
  
  # Return imperfection scores for videos
  # accuracy <- p$posterior[,2]
  accuracy <- p$x[,1]
  list(performance = tab, accuracy = accuracy, lda = r)
}

typeToLDAClass <- function(type) {
  as.factor(ifelse(type == MTP_MODEL, 'ant', 'non-ant'))
}

# Calculates a "one against the rest" LDA, whereby 
# the LDA is used to discriminate ants from non-ants.
# 
# priorWeights - prior weighting for encountering an ant or a non-ant. 
#                Can be used to allow for the relative costs of mis-identifying 
#                an ant (by increasing the first value), or a non-ant (increase the 2nd value)
CalcMotionLDAOATR <- function(stats, videos, priorWeights = c(1, 1)) {
  
  s <- removeStatsNAs(stats[, STATS_NAMES])
  
  # # Handle NAs - create a new statistic which is a boolean to record the fact of an NA, then set the NA to 0
  # #NO_LOCAL_MIN_PLACEHOLDER <- mean(s$first.minimum.C, na.rm = TRUE) # Not sure about this, but doesn't make much difference anyway
  # NO_LOCAL_MIN_PLACEHOLDER <- 0
  # s$noMin <- is.na(s$first.minimum.C)
  # s[is.na(s)] <- NO_LOCAL_MIN_PLACEHOLDER

    s$class <- typeToLDAClass(videos$type)

  # Train on ants and non-ants. Convert prior weights to probabilities
  r <- lda(class ~ ., data = s, prior = priorWeights / sum(priorWeights))
  
  p <- predict(r, newdata = s)

  # Table of actual vs predicted
  tab <- table(factor(videos$type), p$class)
  names(dimnames(tab)) <- c('Actual', 'Predicted')
  
  # Return accuracy scores for videos. Negate value since we want ants to score higher than non-ants
  # accuracy <- p$posterior[,1]
  accuracy <- -p$x[,1]
  
  # Calculate sd for ants and non-ants
  antsSd <- apply(s[s$class == 'ant', -15], 2, sd)
  nonAntsSd <- apply(s[s$class == 'non-ant', -15], 2, sd)
  
  # Calculate means and sd for mimics and non-mimics
  mimicsMean <- apply(s[videos$type == 'mimetic spider', -15], 2, mean)
  mimicsSd <- apply(s[videos$type == 'mimetic spider', -15], 2, sd)
  nonMimicsMean <- apply(s[videos$type == 'non-mimic', -15], 2, mean)
  nonMimicsSd <- apply(s[videos$type == 'non-mimic', -15], 2, sd)
  
  list(performance = tab, accuracy = accuracy, lda = r, predicted = p, 
       antsSd = antsSd, nonAntsSd = nonAntsSd, 
       mimicsMean = mimicsMean, mimicsSd = mimicsSd, 
       nonMimicsMean = nonMimicsMean, nonMimicsSd = nonMimicsSd)
}

ExportMotionImperfection <- function(videos) {
  # Where to export to - the sampleIt app so that it is accessible publically
  csvExportDir <- 'C:/Jim/uni/apps/sampleit/public/ac'
  
  # Export per individual
  d <- aggregate(videos$accuracy, by = list(videos$imageableid, videos$type, videos$species, videos$bodylength, videos$disposition), mean)
  names(d) <- c('imageableid', 'type', 'species', 'bodylength', 'disposition', 'accuracy')
  write.csv(d[order(d$accuracy),], file = file.path(csvExportDir, "accuracy-per-individual-motion.csv"), row.names = FALSE)
  
  d <- aggregate(d$accuracy, by = list(d$type, d$species), mean)
  names(d) <- c('type', 'species', 'accuracy')
  write.csv(d[order(d$accuracy),], file = file.path(csvExportDir, "accuracy-per-species-motion.csv"), row.names = FALSE)
}

fixStatsNames <- function(names) { gsub('\\.', ' ', capwords(names)) }

exportLDASummary <- function(behaviourDir, lda) {
  ldaResult <- lda$lda
  scaledMeanDiff <- (ldaResult$means[2,] - ldaResult$means[1,]) * ldaResult$scaling
  ldaOrder <- order(abs(scaledMeanDiff), decreasing = TRUE)
  d <- data.frame(ldaResult$scaling, 
                  ldaResult$means[1,], 
                  lda$antsSd,
                  ldaResult$means[2,],
                  lda$nonAntsSd,
                  #lda$mimicsMean, 
                  #lda$mimicsSd, 
                  #lda$nonMimicsMean, 
                  #lda$nonMimicsSd, 
                  scaledMeanDiff)
  names(d) <- c('Scaling', 'Ants.mean', 'Ants.sd', 'NonAnts.mean', 'NonAnts.sd', 'Means.diff.scaled')
  #names(d) <- c('Scaling', 'Ants.mean', 'Ants.sd', 'NonAnts.mean', 'NonAnts.sd', 'Mimics.mean', 'Mimics.sd', 'NonMimics.mean', 'NonMimics.sd', 'Means.diff.scaled')
  d <- signif(d, 3)
  write.csv(d[ldaOrder,], file = file.path(behaviourDir, "lda-summary.csv"))
}

exportStatsSummary <- function(behaviourDir, lda, videos, stats) {
  ldaResult <- lda$lda
  # Scaled difference of means indicates the relative contribution of each statistic to the analysis
  scaledMeanDiff <- (ldaResult$means[2,] - ldaResult$means[1,]) * ldaResult$scaling
  ldaOrder <- order(abs(scaledMeanDiff), decreasing = TRUE)
  d <- data.frame(ldaResult$means[1,], 
                  lda$antsSd,
                  lda$mimicsMean, 
                  lda$mimicsSd, 
                  lda$nonMimicsMean, 
                  lda$nonMimicsSd)
  names(d) <- c('Ants.mean', 'Ants.sd', 'Mimics.mean', 'Mimics.sd', 'NonMimics.mean', 'NonMimics.sd')
  
  # Add columns showing whether ants differ from mimics and non-mimics
  s <- removeStatsNAs(stats[, STATS_NAMES])
  type <- factor(videos$type)
  ps <- t(sapply(colnames(s), function(col) {
    t <- TukeyHSD(aov(s[,col] ~ type))
    t$type[,'p adj']
  }))
  p.cols <- c('Ants.Mimics.p', 'Mimics.NonMimics.p', 'Ants.NonMimics.p')
  colnames(ps) <- p.cols
  d <- cbind(d, ps)
  #PValToSigStr(ps)
  
  d <- signif(d, 2)
  write.csv(d[ldaOrder,], file = file.path(behaviourDir, "stats-summary.csv"))
  d
}

findBadlyScaledTracks <- function(videos) {
  kansasCity <- sapply(videos$track, function(t) max(t$y))
  tallest <- kansasCity == max(kansasCity)
  print(videos[tallest,]$file)
}

checkCorrelation <- function(formula, stats) {
  plot(formula, data = stats, pch = 16, cex = .8)
  l <- lm(formula, data = stats)
  abline(l, col = 'red')
  pValue <- summary(l)$coefficients[2,4]
  title(bquote(paste(R^2, "= ", .(summary(l)$adj.r.squared), ", p value = ", .(pValue), " ", .(PValToSigStr(pValue)))))
}

plotStatsCorrelations <- function(stats) {
  # Plot correlations between all of the stats
  names(stats) <- fixStatsNames(names(stats))
  c <- cor(na.omit(stats))
  corrplot(c, addCoef.col = "grey", diag = FALSE, main = 'Correlation matrix')
}

# Scatter plot of first local minimum of directional autocorrelation.
# Used for comparison with the equivalent plot in Shamble et al.
firstMinScatter <- function(videos, stats, stars = TRUE) {
  cols <- typeToCol(videos$type)
  plot(first.minimum.C ~ first.minimum.delta.s, data = stats, pch = 16, col = cols, 
       ylab = expression('C('*Delta*s['min']*')'),
       xlab = expression(Delta*s['min']/'body length'))
  if (stars)
    JStars(stats$first.minimum.delta.s, stats$first.minimum.C, videos$type, col = cols, labels = T)
  types <- sort(AsPlottableMimicTypeFactor(unique(videos$type)))
  legend <- paste0(types, 's')
  legend <- sapply(legend, capsentence)
  legend("topright", legend = legend, pch = 16, col = typeToCol(types), inset = c(.01, .01))
}

ExportMotionAppendexIndividuals <- function(videos, stats, dir, basename) {
  cn <- c('imageableid', 'type', 'order', 'scientificName', 'scientificNameAuthorship', 
          'decimalLatitude', 'decimalLongitude', 'coordinateUncertaintyInMeters', 'elevation', 'locationRemarks', 
          'time', 'day', 'month', 'year', 'bodylength', 'accuracy')
  d <- videos[,cn]
  # Combine stats
  d <- cbind(d, stats[,STATS_NAMES])
  d$type <- AsPlottableMimicTypeFactor(d$type)
  # Rename column imageableid to catalogNumber
  names(d)[names(d) == "imageableid"] <- 'catalogNumber'
  write.csv(d[order(d$type, d$scientificName),], file = file.path(dir, paste0(basename, ".csv")), na = '', row.names = FALSE)
  
  dr <- paste(range(as.Date(DBDate(d))), collapse = ' and ')
  cat(sprintf("Summary of videoed specimens collected between %s\n", dr))
  spt <- aggregate(d$scientificName, list(d$type), function(c) length(unique(c)))
  names(spt) <- c('type', 'count')
  st <- table(d$type)
  st <- st[st > 0]
  summaryTable <- rbind(species = spt$count, individuals = st)
  print(cbind(summaryTable, total = rowSums(summaryTable)))
}

ExportMotionAppendexSpecies <- function(species, dir, basename) {
  d <- species[,c("type", "order", "scientificName", "scientificNameAuthorship", 
                 "bodylength", "accuracy", STATS_NAMES)]
  # Merge in justification for mimicry classification
  d <- cbind(d, .MGetMimicryDecisions(species))
  d$type <- AsPlottableMimicTypeFactor(d$type)
  write.csv(d[order(d$type, d$scientificName),], file = file.path(dir, paste0(basename, ".csv")), na = '', row.names = FALSE)
}

plotStatsDensities <- function(stats) {
  densities <- apply(stats, 2, function(s) density(scale(s), na.rm = TRUE))
  GVPlotDensities(densities, 'Statistics values', 1:length(densities))
}

# Given statistics for individual trajectories, averages them to get average values for species
averageStatsToSpecies <- function(videos, stats) {
  ss <- aggregate(cbind(stats, bodylength = videos$bodylength),
                  by = list(videos$type, videos$order, videos$family, videos$genus, 
                            videos$scientificName, videos$scientificNameAuthorship), 
                  mean)
  names(ss)[1:6] <- c('type', 'order', 'family', 'genus', 'scientificName', 'scientificNameAuthorship')
  sc <- aggregate(videos$scientificName, by = list(videos$scientificName), length)
  names(sc) <- c('scientificName', 'trajectory.count')
  merge(sc, ss)
}

reportLDAPerformance <- function(ldar, input) {
  print(ldar$performance)
  cat("Classification was wrong for \n")
  print(input[typeToLDAClass(input$type) != ldar$predicted$class,]$scientificName)
}

plotCombinedPCA <- function(outDir, stats, videos, species) {
  PlotToPng(JUCFile("stats-pca.png", outDir), function () { 
    JSetParsTemporarily(mfrow=c(1, 2), mar = c(1.8, 1.8, 0.1, 0.1))
    varInd <- plotMotionPCA(stats, videos, labels = FALSE)
    mtext("a) Individuals", line = -1.1, adj = .97)
    cat(sprintf("PC1 & PC2 explain %g%% of variance in individual PCA plot\n", signif(varInd$var * 100, 3))) 
    varSp <- plotMotionPCA(species, species, labels = FALSE); 
    mtext("b) Species", line = -1.1, adj = .97)
    cat(sprintf("PC1 & PC2 explain %g%% of variance in species PCA plot\n", signif(varSp$var * 100, 3))) 
    
    # Export contributions of each variable to 1st 2 principal components
    contrib <- cbind(varSp$contrib, varInd$contrib)
    colnames(contrib) <- c('Species PC1', 'Species PC2', 'Individuals PC1', 'Individuals PC2')
    rowIdx <- order(contrib[,1], decreasing = TRUE)
    contrib <- signif(100 * contrib, 1)
    for (r in 1:nrow(contrib)) {
      for (c in 1:ncol(contrib))
        contrib[r,c] <- paste0(contrib[r,c], '%')
    }
    write.csv(contrib[rowIdx,], file = JUCFile("pca-components.csv", outDir))
  }, width = 600, height = 300)
}

FinalTalkMotionPlots <- function(motionDir, morphoDir, stats, videos) {
  TALK_DIR <- 'C:/Users/jim_m/Google Drive/Uni/Classes/MRes Thesis/Final talk'
  
  ## Motion PCA plot
  PlotToPng(file.path(TALK_DIR, "motion-pca-individual.png"), function () {
    par(mar = c(3, 3, .1, .1))
    plotMotionPCA(stats, videos, fill = TRUE, cex = 2, pt.cex = 2, labels = FALSE)
  }, width = 900, height = 700)
  
  ### Motion morpho plot
  both <- .MnReadMotionMorpho(motionDir, c('accuracy', 'scientificName'), 
                              morphoDir, c('type', 'scientificName',
                                           'dorsalAccuracyMorphometric', 'lateralAccuracyMorphometric',
                                           'dorsalAccuracyMachineLearning', 'lateralAccuracyMachineLearning'))
  .plot <- function(x, y, ...) {
    #.pf <- function(mm) {text(mm[,x], mm[,y], mm$scientificName, cex = 0.8, pos = 1, font = 3)}
    .pf <- function(mm) {}
    f <- reformulate(x, response = y)
    xlim <- extendrange(both[both$type == MTP_SPIDER,x], f = .16) + 0.07
    ylim <- extendrange(both[both$type == MTP_SPIDER,y], f = .04) - 0.01
    l <- MPlotScatterRegression(both, formula = f, types = MTP_SPIDER,
                                xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', 
                                legendPos = NULL, extraFn = .pf, 
                                xlim = xlim, ylim = ylim, ...)[[1]]
    MLabelAccuracyAxis('Morphological accuracy', cex = 1.5)
    MLabelAccuracyAxis(axis = 'y', 'Behavioural accuracy', cex = 1.5)
  }

  PlotToPng(file.path(TALK_DIR, "motion-morpho-corr.png"), function () {
    JSetParsTemporarily(mar = c(4, 4, .1, .1))
    .plot('dorsalAccuracyMorphometric', 'accuracy', cex = 2, lwd = 4)
  }, width = 600, height = 500)  
  
  
  ### Speed morpho plot
  both <- .MnReadMotionMorpho(motionDir, c('mean.speed.bodylength', 'maximum.speed.bodylength', 'scientificName'), 
                              morphoDir, c('type', 'scientificName',
                                           'dorsalAccuracyMorphometric', 'lateralAccuracyMorphometric',
                                           'dorsalAccuracyMachineLearning', 'lateralAccuracyMachineLearning'))
  EILICA <- 'Eilica sp1'
  normal <- both[!(both$scientificName %in% c(EILICA)),]
  .plot <- function(xCol, yCol, ...) { 
    x <- both[,c(xCol, yCol, 'type', 'scientificName')]
    names(x) <- c('x', 'y', 'type', 'scientificName')
    nx <- normal[normal$type == c(MTP_SPIDER), c(xCol, yCol, 'type')]
    names(nx) <- c('x', 'y', 'type')
    .normalLn <- function(mm) { 
      eilica <- mm[mm$scientificName == EILICA,]
      text(eilica$x, eilica$y, 'Eilica', cex = 1.5, pos = 3, font = 3)
    }
    xlim = extendrange(both[both$type == MTP_SPIDER, xCol], f = 0.15) - .3
    ylim = extendrange(both[both$type == MTP_SPIDER, yCol], f = 0.02)
    ll <- MPlotScatterRegression(x, y ~ x, types = MTP_SPIDER, 
                                 yaxt = 'n', legendPos = NULL, xlim = xlim, ylim = ylim, ...)
    .normalLn(x)
    MLabelAccuracyAxis(axis = 'y', 'Morphological accuracy', cex = 1.5)
    title(xlab = 'Walking speed (body lengths/sec)', line = 3, cex.lab = 1.5)
  }
  
  PlotToPng(file.path(TALK_DIR, 'speed-morpho-corr.png'), function() {
    JSetParsTemporarily(mar = c(4.5, 4, .1, .1))
    .plot('mean.speed.bodylength', 'dorsalAccuracyMorphometric', cex = 1.5, cex.lab = 1.5, lwd = 4)
  }, width = 600, height = 500)    
}

#########################################################################

debuggingCrap <- function () {
  ## Plot rotational velocity density
  #plotLineFn(videos, function(t) as.numeric(density(abs(trackVelocities(t)[["angularVelocity"]]))$y))
  
  #file <- "C:/Jim/uni/Classes//Thesis/Behaviour/Motion-videos/cache2b28a073691.csv"
  #plotAllTracksInFile(file)
  
  #debugFile("C:/Jim/uni/Classes//Thesis/Behaviour/Motion-videos/cache2b2821477c72.csv")
  
  #track <- readTrack(file, allSpecimens)
  #cat(sprintf("Stopped for %d%% of the time", stoppedPercent(track)))
  #plotTrack(track)
  
  ## Convert cartesian to polar coordinates
  polar <- track$polar
  ## Differentiate displacement to get speed
  displacement <- diff(polar)
  t1 <- track$Time[2:nrow(track)]
  sm <- sm.spline(t1, Mod(displacement))
  smr <- sm.spline(t1, Arg(displacement))
  speed <- sm$y[2:length(sm$y)] / diff(sm$x)
  angularVelocity <- smr$y[2:length(smr$y)] / diff(smr$x)
  t1 <- sm$x[2:length(sm$x)]
  plot(speed ~ t1, type = 'p', pch = 16)
  
  acc <- diff(speed) / diff(t1)
  t2 <- t1[2:length(t1)]
  plot(acc ~ t2, pch = 16, cex = .5, col = 'red', ylim = c(-100, 250))
  points(speed ~ t1, type = 'l', pch = 16, col = 'blue', cex = .4)
  abline(h = 0, col = MyPallete$midGrey, lty = 2)
  
  plot(track$x, track$y, type = 'l', col = 'black', lwd = 2, main = basename(file), asp = 1)
  
  plot3d(track$x, track$y, track$Time, type = 'l', col = 'black', lwd = 2, main = basename(file), asp = 1)
  
  plot3d(speed, angularVelocity, track$Time, type = 'l', col = 'black', lwd = 2, main = basename(file), asp = 1)
  plot(speed, angularVelocity, type = 'l', col = 'blue', lwd = 2, main = basename(file))
  
  ###
  ## FFT
  
  library(GeneCycle)
  
  plot.frequency.spectrum <- function(X.k, xlimits=c(0,length(X.k))) {
    plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))
    
    ## TODO: why this scaling is necessary?
    plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2] 
    
    plot(plot.data, t="h", lwd=2, main="", 
         xlab="Frequency (Hz)", ylab="Strength", 
         xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
  }
  
  ## Remove trends
  par(mfrow=c(2,1))
  lx <- lm(x ~ Time, data = track)
  ly <- lm(y ~ Time, data = track)
  plot(lx$residuals)
  plot(ly$residuals)
  plot.frequency.spectrum(fft(lx$residuals), xlimits = c(1, 20))
  plot.frequency.spectrum(fft(ly$residuals), xlimits = c(1, 20))
  par(mfrow=c(1,1))
}

#####################################################################################################

# Motion analysis
runMotionAnalysis <- function() {
  outDir <- 'Thesis/Behaviour'
  behaviourDir <- JUCFile(subDir = outDir)
  
  # Read them all in
  videos <- MGetAllTracks()
  
  # Calculate imperfection and output stats and index
  stats <- MMotionStatsForVideos(videos)
  #index <- buildMotionIndexFromStats(stats, videos$type)
  #videos$imperfection <- imperfectionFromStats(index, stats)
  #ldar <- CalcMotionLDANoMimics(stats, videos)
  #videos$accuracy <- ldar$accuracy
  ldar2 <- CalcMotionLDAOATR(stats, videos, priorWeights = c(20, 1))
  videos$accuracy <- ldar2$accuracy
   # cat("LDA performance for trajectories:\n")
   # reportLDAPerformance(ldar2, videos)
  
  # Aggregate to species
  species <- averageStatsToSpecies(videos, stats)
  slda <- CalcMotionLDAOATR(species, species)
  species$accuracy <- slda$accuracy
  cat("LDA performance for species:\n")
  reportLDAPerformance(slda, species)
  
  # Dump species stats
  # write.csv(stats, JUCFile("ind-stats.csv", outDir))
  # write.csv(species, JUCFile("species-stats.csv", outDir))
  
  # generate plots for final talk
  
  
  # Is there evidence that mimics are more ant-like than non-mimics? Assume ants mean < non-mimics mean
  cat(sprintf("\nAre mimics significantly more ant-like than non-mimics?\n"))
  #tt <- t.test(videos[videos$type == MTP_SPIDER,]$accuracy, videos[videos$type == MTP_NON_MIMIC,]$accuracy, alternative = 'greater')
  tt <- t.test(species[species$type == MTP_SPIDER,]$accuracy, species[species$type == MTP_NON_MIMIC,]$accuracy, alternative = 'greater')
  print(tt)
  # Mann-Whitney test
  tt <- wilcox.test(species[species$type == MTP_SPIDER,]$accuracy, species[species$type == MTP_NON_MIMIC,]$accuracy, alternative = 'greater')
  print(tt)

  exportLDASummary(behaviourDir, slda)
  summ <- exportStatsSummary(behaviourDir, slda, species, species)
  cat(sprintf("\nMimics are less ant-like than non-mimics for %s\n", paste(rownames(summ)[sign(summ$NonMimics.mean - summ$Ants.mean) != sign(summ$NonMimics.mean - summ$Mimics.mean)], collapse = ', ')))
  cat(sprintf("Mimics are _more_ ant-like than ants for %s\n\n", paste(rownames(summ)[sign(summ$Mimics.mean - summ$Ants.mean) != sign(summ$NonMimics.mean - summ$Ants.mean)], collapse = ', ')))
  
  # Plot them all
  #plotTracks(videos$track)
  
  PlotToPng(JUCFile("tracks-by-type.png", outDir), function() plotTypeTracks(videos), width = 660, height = 350)
  
  #derivs <- lapply(videos$track, MTrackDerivatives)
  #exportSpeedPlotsForTypes(derivs, videos, behaviourDir)
  #exportAnglePlotsForTypes(videos, behaviourDir)
  
  ExportEverythingForAll(videos, behaviourDir)
  PlotToPng(JUCFile("characteristics.png", outDir), function () { MPlotAllStats(videos) }, width = 1200, height = 1200)
  
  PlotToPng(JUCFile("tracks-by-antness.png", outDir), function () { plotAllTracksByIndex(videos) }, width = 2400, height = 2400)
  
  #PlotToPng(JUCFile("motion-LDA-imperfection.png", outDir), function () {  MpPlotAccuracyDensities(videos$type, ldar$accuracy); abline(v = .05, col = MyPallete$lightGrey) }, width = 800, height = 600)
  PlotToPng(JUCFile("motion-imperfection-ind.png", outDir), function () { MpPlotAccuracyDensities(videos$type, videos$accuracy); abline(v = .05, col = MyPallete$lightGrey) }, width = 800, height = 600)
  PlotToPng(JUCFile("motion-imperfection-species.png", outDir), function () { MpPlotAccuracyDensities(species$type, species$accuracy) }, width = 500, height = 450)
  
  PlotToPng(JUCFile("stats-correlations.png", outDir), function () { plotStatsCorrelations(stats) }, width = 800, height = 800)
  #PlotToPng(JUCFile("first-local-minimum.png", outDir), function () { firstMinScatter(videos, stats) })

  # PlotToPng(JUCFile("stats-pca-individual.png", outDir), function () { plotMotionPCA(stats, videos) }, width = 800, height = 600, res = 100)
  # PlotToPng(JUCFile("stats-pca-species.png", outDir), function () { varE <- plotMotionPCA(species, species); cat(sprintf("PC1 & PC2 explain %g%% of variance in species PCA plot\n", signif(varE * 100, 3))) }, width = 560, height = 380)
  plotCombinedPCA(outDir, stats, videos, species)
  
  ExportMotionImperfection(videos)
  
  ExportMotionAppendexIndividuals(videos, stats, behaviourDir, 'Table S4')
  ExportMotionAppendexSpecies(species, behaviourDir, 'Table S5')
  

  # MpPlotAccuracyDensities(videos, CalcMotionLDA2(stats, videos, priorWeights = c(20, 1))$accuracy, ylim = c(0, 5))

  
  MnCmpMotionMorpho(behaviourDir, MORPHO_DIR, behaviourDir)
  MnCmpSpeedMorpho(behaviourDir, MORPHO_DIR, behaviourDir)

  
  #FinalTalkMotionPlots(behaviourDir, MORPHO_DIR, stats, videos)
  
  PlotToPng(JUCFile("map.png", outDir), function () {
    par(mar = c(0, 0, 0, 0))
    GMapSpecimenSites(videos, xlimFrac = 2, ylimFrac = 2)
  }, width = 650, height = 600)
}

# Only call the function on sourcing if not running within RStudio
if (!JIsRStudio())
  runMotionAnalysis()

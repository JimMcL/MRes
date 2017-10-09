#!Rscript

LoadFns('rediscretization')

# Returns the angles between pairs of steps in the track.
# By default, adjacent pairs are compared.
#
# track - track whose turning angles are to be returned
# deltaS - steps separated by this distance are compared (default 1)
MTurningAngles <- function(track, deltaS = 1) {
  # it's confusing, but angles in displacement are tangents to the path segments
  diff(Arg(track$displacement), lag = deltaS)
}


# Directional autocorrelation functions ##########################################################

## Autocorrelation based on (Shamble et al., 2017)
##
## I have tried to use variable names which are suggestive 
## of the names in the article

LOCAL_MINIMA_WINDOW_SIZE <- 50


# Calculates the autocorrelation of the track for the given delta s.
# Track must have been rediscretized so that all segments have the same length.
# deltaS is specified in number of segments.

# Calculates autocorrelations for delta s values ranging from 1 to length(track)/2.
# Track should already be rediscretized.
MDirectionAutocorrelations <- function(track, deltaSMax = round(nrow(track) / 2)) {

  # The guts of the autocorrelation function
  # Calculates autocorrelation for a single delta s
  .deltaCorr <- function(deltaS, track) {
    # Calculate difference in angle for every pair of segments which are deltaS apart, 
    # take cos, then mean
    c <- sapply(seq(1, length.out = deltaS), function(offset) {
      t <- track[offset:nrow(track),]
      cos(MTurningAngles(t, deltaS))
    })
    mean(unlist(c))
  }
  
  deltaSs <- 1:deltaSMax
  data.frame(deltaS = deltaSs,
             C = sapply(deltaSs, .deltaCorr, track))
}

# Runs the autocorrelation function for direction on a track.
# Track is for rediscretized to 1/10th body length.
# value - data frame with columns deltaS and C. DeltaS units are body lengths
MAutocorrelateTrack <- function(video, rediscretizationFactor = .1) {
  segLen <- video$bodylength * rediscretizationFactor
  # Tracks are actually lists with a single element which is the track! Weird!
  rt <- MRediscretizeTrack(video$track[[1]], segLen)
  corr <- MDirectionAutocorrelations(rt)
  # Convert length units from discretized path length to body lengths
  corr$deltaS <- corr$deltaS * rediscretizationFactor
  corr
}

MFindFirstMinimum <- function(corr) {
  # Ignore local minima if it's the end of the track
  windowSize <- min(length(corr$C) - 1, LOCAL_MINIMA_WINDOW_SIZE)
  minima <- JLocalMaxima(-corr$C, windowSize, endIndex = length(corr$C) - windowSize)
  if (length(minima) > 0) {
    c(deltaS = corr$deltaS[minima][1], C = corr$C[minima][1])
  }
}

MFindFirstMaximum <- function(corr) {
  windowSize <- min(length(corr$C) - 1, LOCAL_MINIMA_WINDOW_SIZE)
  # Ignore local maxima if it's the start of the track
  maxima <- JLocalMaxima(corr$C, LOCAL_MINIMA_WINDOW_SIZE, startIndex = windowSize)
  if (length(maxima) > 1) {
    c(deltaS = corr$deltaS[maxima][2], C = corr$C[maxima][2])
  }
}

MPlotAutocorrelations <- function(corrs, type, title = 'Autocorrelation') {
  xlim <- range(sapply(corrs, function(c) c$deltaS))
  ylim <- range(sapply(corrs, function(c) c$C))
  plot(NULL, 
       xlim = xlim, ylim = ylim,
       main = title, 
       xlab = expression(paste(Delta, 's / body length')),
       ylab = expression(paste('C(', Delta, 's)')))
  idx <- 1
  sapply(corrs, function(c) {
    col <- typeToCol(type[idx])
    idx <<- idx + 1
    lines(c$deltaS, c$C, col = col)
    minima <- JLocalMaxima(-c$C, LOCAL_MINIMA_WINDOW_SIZE)
    # Ignore local minima if it's the end of the track
    if (length(minima) > 0 && minima[1] < length(c$C) - LOCAL_MINIMA_WINDOW_SIZE) {
      points(c$deltaS[minima][1], c$C[minima][1], pch = 16, col = col, lwd = 2)
      points(c$deltaS[minima][1], c$C[minima][1], col = 'black', lwd = 2)
    }
  })
  
  if (length(unique(type)) > 1) {
    legend("topright", legend = unique(type), lwd = 1, col = typeToCol(unique(type)))
  }
  NULL
}

# Returns the directional autocorrelation for a video.
# Attempts to read it from a pre-existing CSV file, otherwise calculates it 
# and writes it to the CSV file before returning.
MGetDirnAutocorrelation <- function(video, verbose = TRUE) {
  csvFile <- gsub('.csv$', '-corr.csv', video$trackFile)
  corrInfo <- file.info(csvFile)
  trackInfo <- file.info(video$trackFile)
  # If correlation file exists, is non-empty, and is newer than the track file...
  if (!is.na(corrInfo$size) && corrInfo$size > 0 && corrInfo$mtime > trackInfo$mtime) {
    # Just read it
    corr <- read.csv(csvFile)  
  } else {
    # Otherwise (re)calculate it
    if (verbose) 
      cat(sprintf("Calculating autocorrelation for video %d, specimen %d %s\n", video$id, video$imageableid, video$scientificName))
    corr <- MAutocorrelateTrack(video)
    # Save it for next time
    write.csv(corr, csvFile, row.names = FALSE)
  }
  corr
}
  
.MExtractFirstMinima <- function(corrs, types = NULL) {
  deltaS <- numeric()
  C <- numeric()
  type <- numeric()
  for(i in 1:length(corrs)) {
    min <- MFindFirstMinimum(corrs[[i]])
    if (!is.null(min)) {
      idx <- length(C) + 1
      deltaS[idx] <- min['deltaS']
      C[idx] <- min['C']
      type[idx] <- ifelse(is.null(types), 0, types[i])
    }
  }
  data.frame(deltaS, C, type)
}

# TESTING ########################################################################################

testAutoCorr <- function() {
  # Don't rediscretize, which means that points are evely spaced along the x-axis but at different distance,
  # result is that correlation is aligned with track on x-axis
  track1 <- .MtrackFillIn(data.frame(x = 1:100, y = 50 * cos((1:100) / 5), Time = c(1:100)))
  track2 <- .MtrackFillIn(data.frame(x = 1:100, y = 6 * cos((1:100) / 5), Time = c(1:100)))
  plotTrack(track1, col = 'red')
  lines(y ~ x, data = track2, col = 'blue')
  corr1 <- MDirectionAutocorrelations(track1, deltaSMax = nrow(track1))
  corr2 <- MDirectionAutocorrelations(track2, deltaSMax = nrow(track2))
  lines(corr1, col = 'red', lty = 2)
  lines(corr2, col = 'blue', lty = 2)
}

test <- function() {
    #windows()

    
  corrs <- apply(tracks, 1, MGetDirnAutocorrelation)
  MPlotAutocorrelations(corrs, tracks$type)
  
  antIdx <- which(tracks$type == MTP_MODEL)
  spiderMimicIdx <- which(tracks$type == MTP_SPIDER)
  insectMimicIdx <- which(tracks$type == MTP_INSECT)
  nonMimicIdx <- which(tracks$type == MTP_NON_MIMIC)
  
  m <- .MExtractFirstMinima(corrs, tracks$type)
  plot(C ~ deltaS, data = m, pch = 16, col = typeToCol(minType))
  legend("topright", legend = unique(tracks$type), pch = 16, col = typeToCol(unique(tracks$type)))
  
  .plotCHull <- function(corrs, col) {
    m <- .MExtractFirstMinima(corrs)
    ch <- chull(m[,1:2])
    # Close the polygon
    ch <- c(ch, ch[1])
    lines(m$deltaS[ch], m$C[ch], col = col)
  }
  .plotCHull(corrs[antIdx], "red")
  .plotCHull(corrs[spiderMimicIdx], "purple")
  .plotCHull(corrs[nonMimicIdx], "blue")
  
  .mm <- function(corrs) { as.matrix(.MExtractFirstMinima(corrs))[,1:2] }
  kde.test(.mm(corrs[antIdx]), .mm(corrs[nonMimicIdx]))$pvalue
  kde.test(.mm(corrs[antIdx]), .mm(corrs[spiderMimicIdx]))$pvalue
  kde.test(.mm(corrs[nonMimicIdx]), .mm(corrs[spiderMimicIdx]))$pvalue
  
  
    numerator <- 10
    bodyLen <- 4
    x <- 0:200
    .mt <- function(y) {
        MRediscretizeTrack(.MtrackFillIn(data.frame(x = x, y = y, Time = x)), bodyLen / numerator)
    }

    .plotTracks <- function(tracks) {
        #dev.new()
        xlim <- range(sapply(tracks, function(t) t$x / bodyLen))
        ylim <- range(sapply(tracks, function(t) t$y))

        par(mfrow = c(2, 1))
        
        plot(NULL, xlim = xlim, ylim = ylim, main = 'Tracks')
        col <- c('red', 'blue', MyPallete$darkGreen, MyPallete$darkOrange, MyPallete$midPurple, 'black')
        i <- 1
        sapply(tracks, function(t) {
            lines(t$x / bodyLen, t$y, col = col[i])
            points(t$x[1], t$y[1], pch = 16, cex = .8, col = col[i])
            i <<- i + 1
        })

        corr <- lapply(tracks, MDirectionAutocorrelations)
        xlim <- range(sapply(corr, function(c) c$deltaS / numerator))
        ylim <- range(sapply(corr, function(c) c$C))
        plot(NULL, xlim = xlim, ylim = ylim, main = 'Autocorrelations')
        i <- 1
        sapply(corr, function(c) {
            lines((c$deltaS / numerator), c$C, col = col[i])
            i <<- i + 1
        })
    }

    tracks <- list(.mt(sin(x / 6)),
                   .mt(4 * sin(x / 6)),
                   .mt(sin(x / 9)),
                   .mt(4 * sin(sqrt(x))))
    .plotTracks(tracks)

    PauseUntilWindowIsClosed()
}

#test()

# E-MAX ########################################################################################

## E-max from Cheung et al., (2007)

.Mbeta <- function(points) {
  # Calculate difference in angle for every pair of segments which are deltaS apart, 
  # take cos, then mean
    mean(cos(diff(Arg(points))))
}

M.Emax <- function(track, eMaxB = FALSE) {
  
  # Rediscretize to a constant segment length
  #rt <- MRediscretizeTrack(track, 1)

  # Smooth (loses time values since there may be fewer points after smoothing)
  x <- sm.spline(track$x, track$y)
  st <- .MtrackFillIn(data.frame(x = x$x, y = x$y))

  rt <- MRediscretizeTrack(st, 1)
  
  b <- .Mbeta(rt$polar)
  
  # If it's E max b, multiple by mean path length
  f <- ifelse(eMaxB, 1, 1)
  
  f * b / (1 - b)
}

test <- function() {
  track <- .MtrackFillIn(data.frame(x = c(0, 0, -1), y = c(0, 1, 1), Time = c(1:3)))
  track <- .MtrackFillIn(data.frame(x = c(0, 0, -1), y = c(0, 1, 0), Time = c(1:3)))
  
  plotTrack(track)  
  rt <- MRediscretizeTrack(track, 1)
  lines(y ~ x, data = rt, type = 'l', col = 'red')
  points(y ~ x, data = rt, col = 'blue', pch = 16, cex = .5)

  .Mbeta(rt$displacement)
}

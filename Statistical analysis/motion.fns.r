# Functions for reading and analysing CSV tracks produced from videos
library(pspline)
library(tools)
LoadFns("sample-dbs")

# Tracks are data frames with the following values:
# x, y - cartesian coordinates of points in track
# polar - complex representation of x, y
# displacement - 0 + polar[i] - polar[i - i]. Initial 0 is so that length(displacement) == length(polar)


# Reads a set of points from a file. The points come from multiple tracks 
# due to noise in the video conversion process. 
# The longest track is the one we are interested in
#
# Value - data frame with values x & y, and an attribute "numFrames" which records the number of frames in the source video
.MreadPoints <- function(file) {
  points <- read.csv(file, comment.char = '#')
  
  # Save the number of frames in the file in case the track doesn't extend until the end
  maxFrame <- max(points$Frame)
  
  # Convert left-handed coord-system to right (or is it the other way around?)
  points$y <- max(points$y, na.rm = TRUE) - points$y

  # Pick the track with the largest displacement between start and finish
  .pickBestTrack <- function() {
    polar <- complex(real = points$x, imaginary = points$y)
    tids <- unique(points$TrackId)
    dists <- sapply(tids, function(tid) MTrackDistance(polar[points$TrackId == tid]))
    bestId <- tids[which(dists == max(dists))]
    bestId
  }
  tid <- .pickBestTrack()

  points <- points[points$TrackId == tid,]

  # Unfortunately can't smooth like this because it can change the number of points,
  # so correlation between points and Times is lost. 
  # Presumably no. of x & y points can also differ
  # if (smooth) {
  #   points$x <- sm.spline(points$Time, points$x)$y
  #   points$y <- sm.spline(points$Time, points$y)$y
  # }
  
  # Save number of frames
  attr(points, 'numFrames') <- maxFrame
  
  points  
}

.MtrackFillIn <- function(track) {
  # Get polar coordinates
  track$polar <- complex(real = track$x, imaginary = track$y)
  
  # Calculate displacements from each point to the next
  track$displacement <- c(0, diff(track$polar))
  # Get times associated with displacements, with the first point at time 0
  if ('Time' %in% names(track))
    track$displacementTime <- track$Time[1:nrow(track)] - track$Time[1]

  track  
}

# Reads in a track, given its track file name and mimic type
# 
# translate - if TRUE, the track is translated so that its initial point is at (0, 0)
MReadTrack <- function(file, type, translate) {
  track <- .MreadPoints(file)
  if (translate) {
    track$x <- track$x - track$x[1]
    track$y <- track$y - track$y[1]
  }

  # Fill in polar coords, displacement
  track <- .MtrackFillIn(track)

  attr(track, 'file') <- file
  attr(track, 'type') <- type
  track
}

# Calculates speed and linear acceleration over time.
# If smooth is TRUE, the returned values may not have the same number of data points as the input track.
MTrackDerivatives <- function(track, smooth = TRUE) {
  # Note that displacements are the (polar) displacements from 1 point to the next
  d <- Mod(track$displacement)
  t <- track$displacementTime
  
  if (smooth) {
    # Smooth to reduce noise
    sm <- sm.spline(t, d, norder = 3)
    d <- sm$ysmth
    # Smoothing can produce "displacements" < 0. Is this correction dangerous?
    d[d < 0] <- 0
    t <- sm$x
  }
  
  # Calculate speed
  v <- d[2:length(d)] / diff(t)
  vt <- t[2:length(t)]
  # Calculate linear acceleration
  a <- diff(v) / diff(vt)
  at <- vt[2:length(vt)]
  
  list(speed = v, speedTimes = vt, acceleration = a, accelerationTimes = at)
}

# Calculates the cumulative length of a set of points in polar form,
# or the distance travlled along a path
MTrackLength <- function(polarPoints) {
  sum(Mod(diff(polarPoints)))
}

# Calculates the distance between the start and end of a set of points in polar form.
# Also called the diffusion distance, or bee-line from start to finish.
MTrackDistance <- function(polarPoints) {
  MTrackLength(polarPoints[c(1,length(polarPoints))])
}

# Rotates a track so that angle(finish - start) == angle
MTrackRotate <- function(track, angle = 0) {
  # Calculate current orientation
  orient <- Arg(track$polar[length(track$polar)] - track$polar[1])
  # Calculate required rotation
  alpha <- angle - orient
  # Rotation matrix
  rm <- matrix(c(cos(alpha), sin(alpha), -sin(alpha), cos(alpha)), ncol = 2)

  # New track is old track rotated
  nt <- as.data.frame(t(rm %*% (t(track[,c('x', 'y')]))))
  colnames(nt) <- c('x', 'y')
  nt <- .MtrackFillIn(nt)
  
  nt
}

# Calculate percentage of completely unmoving frames 
stoppedPercent <- function(track) {
  
  stoppedIdx <- which(track$ValueChanged == "false")
  # Sequences of unchanged positions are not written to the CSV file,
  # only the first and last frames are, so count total number of stopped frames
  stoppedFrameCount <- track$Frame[stoppedIdx] - track$Frame[stoppedIdx - 1] + 1
  round(sum(stoppedFrameCount) * 100 / attr(track, 'numFrames'))
}

# Reads all of the specified tracks
MReadVideoTracks <- function(trackVideos) {
  mimicTypes <- mimicType(trackVideos)
  
  .rowToTrack <- function(rowIndex) {
    MReadTrack(trackVideos[rowIndex,]$trackFile, mimicTypes[rowIndex], TRUE)
  }

  lapply(1:nrow(trackVideos), .rowToTrack)
}

# Reads in track files, and combines them with specimen data from the sampleIt database.
# Track files are assumed to be named "<video id>.csv"
MGetAllTracks <- function(trackDir = file.path(JUCFile(subDir = 'Thesis/Behaviour'), 'Motion-videos')) {
  # When the video files are downloaded, they need to be accessible to a person for track extraction, 
  # so put them in a known location and give them meaningful names. To do that, define a function which 
  # takes information about the photos to be downloaded, and returns _a function_ which returns the 
  # files names to download them to
  .tempfileFnFn <- function(photos) {
    # The function has to accept the same arguments as tempfile
    function(pattern, tmpdir, fileext) {
      # Just call the file <id>.<ext>
      file.path(tmpdir, paste0(photos$id, fileext))
    }
  }
  
  # Get info about all motion videos
  allVideos <- DbGetPhotoData('ptype=Motion&ftype=video', dir = trackDir, tempfileFnFn = .tempfileFnFn)
  allVideos$trackFile <- file.path(trackDir, paste0(allVideos$id, '.csv'))

  # Find videos with a track file
  existingNames <- list.files(trackDir, '[0-9].csv', full.names = FALSE)
  withTracks <- which(basename(allVideos$trackFile) %in% existingNames)
  trackVideos <- allVideos[withTracks,]

  # Read them in
  trackVideos$track <- MReadVideoTracks(trackVideos)

  trackVideos
}

MpPlotAccuracyDensities <- function(type, accuracy, subset = NULL, showSpeciesCount = TRUE, legendPos = 'topright', legendInsetFn = JInset, extrasFn = NULL, xlabFn = MLabelAccuracyAxis, mar = c(3, 2, 0.2, 0.2), ...) {

  if (!is.null(subset)) {
    type <- type[subset]
    accuracy <- accuracy[subset]
  }
  
  types <- unique(type)
  densities <- lapply(types, function(tp) density(accuracy[type == tp], na.rm = TRUE))
  names(densities) <- types

  .plotDensityVert <- function(density, x, ...) {
    .yForX <- function(x) { density$y[which.min(abs(density$x - x))] }
    segments(x, 0, x, .yForX(x), ...)
  }
  .plotTypeLines <- function(tp) {
    a <- accuracy[type == tp]
    d <- densities[[tp]]
    col = typeToCol(tp)
    m <- mean(a, na.rm = TRUE)
    sda <- sd(a, na.rm = TRUE)
    .plotDensityVert(d, m, col = col)
    .plotDensityVert(d, m + sda, lty = 3, col = col)
    .plotDensityVert(d, m - sda, lty = 3, col = col)
  }
  
  JSetParsTemporarily(mar = mar)
  JPlotDensities(densities, col = typeToCol(types), lty = typeToLty(types), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', ...)
  title(ylab = 'Probability Density', line = .5)
  xlabFn(axis = 'x')

  sapply(names(densities), .plotTypeLines)
  if (!is.null(extrasFn))
    extrasFn(densities)
  
  ltypes <- sort(AsPlottableMimicTypeFactor(types))
  fmt <- ifelse(showSpeciesCount, '%s (n = %d)', '%s')
  labels <- sapply(ltypes, function(tp) sprintf(fmt, typeToLabel(tp), sum(as.character(type) == as.character(tp))))
  legend(legendPos, legend = labels, inset = legendInsetFn(), col = typeToCol(ltypes), lty = typeToLty(ltypes), lwd = 2)
}




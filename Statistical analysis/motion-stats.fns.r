#!Rscript
# 
# Functions to calculate statistics for animal tracks extracted from videos.

#library('effsize')
LoadFns('motion')
LoadFns('straightness')

# Velocity that is counted as not moving, in mm / sec
STOPPED_SPEED <- 5


MMotionStatsForVideo <- function(video) {
  track <- video$track[[1]]
  deriv <- MTrackDerivatives(track)
  
  # Find transitions between stopped/moving
  # To simplify calculations, assume stopped at the start and end of the track
  s0 <- c(0, deriv$speed, 0)
  t0 <- c(deriv$speedTimes, deriv$speedTimes[length(deriv$speedTimes)])
  # Define stopped as <= a constant
  movingChanges <- diff(s0 > STOPPED_SPEED)
  # Now movingChanges contains 1 if motion started, -1 if it stopped
  stopTimes <- t0[which(movingChanges == -1)]
  startTimes <- t0[which(movingChanges == 1)]
  # Calculate  number of seconds moving for each interval
  movingIntervals <- stopTimes - startTimes
  # and total number of seconds stopped for each interval
  stoppedIntervals <-  0
  if(length(startTimes) > 1)
    stoppedIntervals <- startTimes[2:length(startTimes)] - stopTimes[1:(length(stopTimes) - 1)]
  # Total time
  runTime <- deriv$speedTimes[length(deriv$speedTimes)] - deriv$speedTimes[1]

  # Path segments when moving
  moving <- deriv$speed[deriv$speed > STOPPED_SPEED]
  
  pathLen <- 1  
  rt <- MRediscretizeTrack(track, pathLen)

  sinuosity <- 1.18 * sd(diff(Arg(tail(rt$displacement, -1)))) / sqrt(pathLen)

  # Called straightness index or r by Batschelet, length of mean vector of unit vectors along path segments sampled at equal time intervals
  # Convert track displacement to unit vectors (track is sampled at equal time intervals while the animal is moving)
  uv <- track$displacement[Mod(track$displacement) > 0]
  uv <- uv / Mod(uv)
  r <- Mod(sum(uv)) / length(uv) / pathLen
  # Modified straightness index - rediscretize to unit segment length to remove time/speed component
  # Equivalent to straightness, but rounding errors give different numbers.
  # I don't know which is theoretically better
  r2 <- Mod(sum(rt$displacement)) / nrow(rt) / pathLen
  # Called d by Batschelet
  straightness <- MTrackDistance(rt$polar) / MTrackLength(rt$polar)
  
  corr <- MGetDirnAutocorrelation(video)
  min <- unname(MFindFirstMinimum(corr))
  # Try using distance from (0, 1) to first minimum
  firstMinDist <- dist(rbind(c(0, 1), min))[1]
  min <- if(is.null(min)) {
    c(NA, NA)
  } else {
    min
  }

  # I think I'll ignore first maximum since Shamble use first minimum,
  # even though I suspect it is more meaningful than the first minimum
  max <- MFindFirstMaximum(corr)
  # Try using distance from (0, 1) to first maximum
  firstMaxDist <- dist(rbind(c(0, 1), max))[1]
  max <- if(is.null(max)) {
    c(NA, NA)
  } else {
    max
  }
  
  eMax <- M.Emax(track)
  
  # 'Special' standard deviation which returns 0 for sd of a single value
  .sd <- function(v) ifelse(length(v) == 1, 0, sd(v))
  # 'Special' coefficient of variation which returns 0 for if mean and sd are both 0
  .CV <- function(v) { mv <- mean(v); sv <- .sd(v); ifelse(mv == 0 && sv == 0, 0, sv / mv) }
  
  # Note that not all statistics calculated here are part of the analysis; see STATS_NAMES
  c(
    maximum.speed = max(deriv$speed),
    maximum.speed.bodylength = max(deriv$speed) / video$bodylength, 
    mean.speed = mean(moving), 
    mean.speed.bodylength = mean(moving) / video$bodylength, 
    speed.CV = .CV(moving), 
    moving.duration.mean = mean(movingIntervals),
    moving.duration.CV = .CV(movingIntervals),
    stopped.duration.mean = mean(stoppedIntervals), 
    stopped.duration.CV = .CV(stoppedIntervals),
    proportion.time.moving = sum(movingIntervals) / runTime,
    r2 = r2,
    straightness = straightness,
    sinuosity = sinuosity,
    E.max = eMax,
    first.minimum.C = ifelse(is.null(min), NA, min[2]),
    first.minimum.delta.s = min[1]
  )
}

MMotionStatsForVideos <- function(videos) {
  as.data.frame(t(JApplyToRows(videos, MMotionStatsForVideo)))
}

PValToSigStr <- function(p) {
  ifelse(p < .001, '***',
         ifelse(p < .01, '**',
                ifelse(p < .05, '*', 
                       ifelse(p < .1, '.', ' ')
                )
         )
  )
}

# Plots all candidate motion statistics in 1 big plot
MPlotAllStats <- function(videos) {
  stats <- MMotionStatsForVideos(videos)
  
  oldPars <- par(no.readonly = T)
  on.exit(par(oldPars))
  ncols <- round(sqrt(ncol(stats)))
  nrows <- ceiling(ncol(stats) / ncols)
  par(mfrow = c(nrows, ncols))
  
  .tidyText <- function(t) { capsentence(gsub('\\.', ' ', t)) }
  
  for (col in colnames(stats)) {
    # Weird use of AsPlottableMimicTypeFactor is to remove unused "insect mimic" level and sort appropriately
    plot(stats[,col] ~ AsPlottableMimicTypeFactor(videos$type), main = .tidyText(col), xlab = '', ylab = col)
    
    # Stats - determine whether ants are different from non-mimics
    tt <- t.test(stats[videos$type == MTP_MODEL,col], stats[videos$type == MTP_NON_MIMIC,col])
    
    # For info, also check whether 3 types differ
    an <- anova(lm(stats[,col] ~ factor(videos$type)))
    title(sub = sprintf("ants - non-mimics != 0? %g %s, all 3 differ? %g %s\n", tt$p.value, PValToSigStr(tt$p.value), an$`Pr(>F)`[1], PValToSigStr(an$`Pr(>F)`[1])))
    
    #cat(sprintf("%s,%g,%s\n", col, tt$p.value, PValToSigStr(tt$p.value)))
  }
}

# Obsolete index functions, replaced by use of LDA
# buildMotionIndexFromStats <- function(stats, types) {
#   
#   ants <- which(types == MTP_MODEL)
#   nonMimics <- which(types == MTP_NON_MIMIC)
#   mimics <- which(types == MTP_SPIDER | types == MTP_INSECT)
#   
#   pLevels <- c(.001, .05)
#   esLevels <-  c(.5, .8)
#   
#   .weightfromPandES <- function(p, es) {
#     # Use a statistic if there is a significant difference between ants and non-mimics.
#     pWeight <- length(pLevels) - findInterval(p, pLevels)
#     esWeight <- findInterval(abs(es), esLevels)
#     w <- pWeight + esWeight - 2
#     ifelse(w > 0, w, 0)
#   }
#   
#   index <- data.frame()
#   # For each potential statistic...
#   for (col in names(stats)) {
#     statsAnts <- stats[ants, col]
#     statsNonMimics <- stats[nonMimics, col]
#     statsMimics <- stats[mimics, col]
#     
#     # Determine whether ants are different from non-mimics
#     tt <- t.test(statsAnts, statsNonMimics)
#     #tt <- wilcox.test(statsAnts, statsNonMimics)
#     
#     # Determine effect size - i.e. are ants different than non-mimics
#     es <- cohen.d(na.omit(statsAnts), na.omit(statsNonMimics))
#     
#     antsMean <- mean(statsAnts, na.rm = TRUE)
#     nonMimicsMean <- mean(statsNonMimics, na.rm = TRUE)
#     mimicsMean <- mean(statsMimics, na.rm = TRUE)
#     
#     # Mimics are more ant-like than non-mimics if ants and mimics lie on the same side of non-mimics
#     mimicsMoreAntLike <- sign(mimicsMean - nonMimicsMean) == sign(antsMean - nonMimicsMean)
#     
#     index <- rbind(index, data.frame(
#       col,
#       .weightfromPandES(tt$p.value, es$estimate),
#       tt$p.value,
#       length(pLevels) - findInterval(tt$p.value, pLevels),
#       es$estimate,
#       findInterval(abs(es$estimate), esLevels),
#       antsMean,
#       mimicsMean, 
#       nonMimicsMean,
#       mimicsMoreAntLike
#     ))
#   }
#   colnames(index) <- c('statistic', 'weight', 'p.value', 'p.weight', 'cohens.d', 'd.weight', 'ants.mean', 'mimics.mean', 'non.mimics.mean', 'mimics.more.ant.like')
#   index
# }
# 
# imperfectionFromStats <- function(index, stats) {
#   # Assume NA is average non-ant-like
#   .naTo1 <- function(v) { v[is.na(v)] <- 1; v }
#   
#   # Sum scores for each row and each statistic
#   #
#   # JApplyToRows loops over each row in index, and inner function
#   # calculates scores for that statistic for all rows. Result is a matrix with 
#   # a row for each video, and a column for each statistic. Summing across
#   rowSums(JApplyToRows(index, function(idxRow) {
#     # Calculate imperfection score for this statistic. 
#     # Value is linearly interpolated between 0 for ant-like to 1 for non-mimic-like
#     # Interpolated value is then weighted according to index
#     idxRow$weight * .naTo1((stats[,idxRow$statistic] - idxRow$ants.mean) / (idxRow$non.mimics.mean - idxRow$ants.mean))
#   }))
# }


# Testing ######################################################################

test <- function() {
  .calculateImperfection <- function(videos) {
    stats <- MMotionStatsForVideos(videos)
    index <- buildMotionIndexFromStats(stats, videos$type)
    imperfectionFromStats(index, stats)
  }

  videos <- MGetAllTracks()
  videos$imperfection <- MICalculateImperfection(videos)
  
  i <- videos$imperfection
  videos <- videos[order(i),]

  print(videos[,c('type', 'scientificName', 'imperfection', 'id')])
  
  bySpecies <- aggregate(videos$imperfection, by = list(videos$species, videos$type), mean)
  names(bySpecies) <- c('species', 'type', 'imperfection')
  bySpecies <- bySpecies[order(bySpecies$imperfection),]

  MpPlotAccuracyDensities(videos, i)
}

#test()

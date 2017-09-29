# Some basic utilities
# 
# Path functions assume variables:
# JSecureDir, JScriptDir, JInsecureDir, JDataDir, JClassDir, JCacheDir
# (currently I set these in Windows in Sys.getenv("R_USER")/.Rprofile)

# Returns the full path to a data file (not backed up).
# Optionally creates the containing directory
JUDFile <- function(filename = "", subDir = "", createDir = FALSE) {
  dir <- file.path(JDataDir, subDir)
  if (createDir && !dir.exists(dir))
    dir.create(dir, recursive = TRUE)
  file.path(dir, filename)
}

# Returns the full path to a file which is on my google drive, i.e. backed up
JUSecureFile <- function(filename = "", subDir = "") {
  file.path(JSecureDir, subDir, filename)
}

# Returns the full path to a script file (backed up)
JUSFile <- function(filename, subDir = "") {
  file.path(JSecureDir, "R", subDir, filename)
}

# Returns the full path to a file in my uni Classes dir (not backed up)
JUCFile <- function(filename = "", subDir = "") {
  file.path(JClassDir, subDir, filename)
}

# Returns the full path to a file in my uni dir (not backed up)
JUUniFile <- function(filename = "", subDir = "") {
  file.path(JInsecureDir, subDir, filename)
}

# Returns a value specified on the comand line, or else NULL.
# Values are specified on the command line as eg:
# ... -v <value> ...
# or
# ... --long-name <value> ...
#
# Example usage: 
# host <- JCLIValue('h', 'host')
JCLIValue <- function(shortName = NA, longName = NA, args = commandArgs(TRUE)) {
  hostFlagIdx <- NA
  if (!is.na(shortName))
    hostFlagIdx <- match(paste0('-', shortName), args)
  if (is.na(hostFlagIdx) && !is.na(longName))
    hostFlagIdx <- match(paste0('--', longName), args)
  if (!is.na(hostFlagIdx) && length(args) > hostFlagIdx) {
    args[hostFlagIdx + 1]
  }
}

# Returns TRUE if an option was specified on the comand line.
# Options are specified on the command line as eg:
# ... -o ...
# or
# ... --long-name <value> ...
#
# Example usage: 
# verbose <- JCLIOption('v', 'verbose')
JCLIOption <- function(shortName = NA, longName = NA, args = commandArgs(TRUE)) {
  (!is.na(shortName) && !is.na(match(paste0('-', shortName), args))) ||
    (!is.na(longName) && !is.na(match(paste0('--', longName), args)))
}

# Returns TRUE if running inside RStudio
JIsRStudio <- function() {
  Sys.getenv("RSTUDIO") == "1"
}

# Searches recursively through the specified root directories for a file with the specied name
FindFirstFile <- function(name, rootDirs) {
  list.files(rootDirs, pattern = name, full.names = T, recursive = T)
}

# Searches for a file named <name.fns.r> in the directory tree rooted at JSecureDir
LoadFns <- function(name) {
  file <- FindFirstFile(paste0(name, ".fns.r"), c(JSecureDir))
  if (length(file) == 0) {
    warning(paste("Unable to locate function file:", name))
  } else {
    source(file)
  }
}

JTrim <- function (x) gsub("^\\s+|\\s+$", "", x)

cleanup <- function () rm(list=ls(envir = .GlobalEnv), envir = .GlobalEnv)

# An apply-like function which iterates over rows in a data frame
JApplyToRows <- function(df, fn, ...) {
  sapply(1:nrow(df), function(idx) fn(df[idx,], ...))
}

# This doesn't belong in this file
OccListAsDataFrame <- function(l) {
  f <- as.data.frame(t(matrix(unlist (l), nrow=4)), stringsAsFactors = F)
  f <- setNames(f, c("name", "occ", "inAus", "note"))
  f[order(as.numeric(f$occ), decreasing=T),]
}

# Captialises the first word in a single string
# Copied (and changed) from the doc for toupper
capsentence <- function(s, strict = FALSE) {
  s <- as.character(s)
  paste(toupper(substring(s, 1, 1)), 
        {s <- substring(s, 2); if(strict) tolower(s) else s},
        sep = "", collapse = " " )
}

capwords <- function(s, strict = FALSE) {
  s <- as.character(s)
  sapply(strsplit(s, split = " "), capsentence, USE.NAMES = !is.null(names(s)))
}

# "Private"
.JplotToDevice <- function(filename, plotFn, onlyIfDoesntExist, openDeviceFn) {
  if (!onlyIfDoesntExist || !file.exists(filename)) {
    openDeviceFn()
    tryCatch({
      plotFn()
    }, finally = {
      dev.off()
    })
  }
}

# Default resolution seems to be 72, so to increase image pixel size without decreasing text size, 
# line width etc, increase resolution accordingly.
# Doesn't work with ggplot since it must be evaluated to work.
# Try using ggsave or print(plotFn())
PlotToPng <- function(filename, plotFn, width=800, height=600, res = NA, onlyIfDoesntExist = F) {
  .JplotToDevice(filename, plotFn, onlyIfDoesntExist, function () {
    # type = 'cairo' seems to produce _much_ nicer graphics with antialiasing
    png(filename, width = width, height = height, type = 'cairo', res = res)
  })
}

# Beware: width and height are in inches, not pixels!
PlotToPDF <- function(filename, plotFn, width=8, height=6, onlyIfDoesntExist = F) {
  .JplotToDevice(filename, plotFn, onlyIfDoesntExist, function () {
    pdf(filename, width = width, height = height)
  })
}

# Call this instead of par().
# Sets graphical parameters temporarily by calling saving existing settings, calling par(...), 
# then restoring the original settings when the calling function returns.
# Only call this function once from a function, or the second call will overwrite the original settings.
JSetParsTemporarily <- function(...) {
  oldPars <- par(no.readonly = T)
  # Some ugly R jiggery pokery.
  # Restore original settings when the calling function exits, not when this one exits
  do.call("on.exit", list(substitute(on.exit(par(oldPars))), add=TRUE), envir=parent.frame())
  par(...)
}

# Plots a list of densities as lines
JPlotDensities <- function(densities, cols, lty = 1, lwd = 2, add = FALSE, xlim = NULL, ylim = NULL, includeInX = numeric(0), ...) {
  
  # Create empty plot
  if (!add) {
    if (is.null(xlim))
      xlim <- range(lapply(densities, function(d) d$x), na.rm = TRUE)
    xlim <- range(c(xlim, includeInX))
    if (is.null(ylim))
      ylim <- range(lapply(densities, function(d) d$y), na.rm = TRUE)
    plot(NULL, xlim = xlim, ylim = ylim, ...)
  }
  
  # Recycle lty if it's a single number
  if (length(lty) == 1)
    lty <- rep(lty, length.out = length(densities))

  # Plot densities as lines
  i <- 1
  for(d in densities) {
    lines(d, col = cols[i], lty = lty[i], lwd = lwd)
    i <- i + 1
  }
}

# Returns a legend inset, given input values in inches
JInset <- function(x = 0.1, y = 0.1) {
  c(x, y) / par()$pin
}


ScaleToRange <- function(values, from = 0, to = 1, ...) {
  range <- to - from
  (range * (values - min(values, ...)) / (max(values, ...) - min(values, ...))) + from
}

# Returns indices of local maxima
# To obtain local minima, call localMaxima(-x)
#
# v - vector of values
# window - number of points on each side which defines what counts as a local maxima
# startIndex - index of first point which can qualify as a maximum
# endIndex - index of last point which can qualify as a maximum
JLocalMaxima <- function(v, window = 1, startIndex = 1, endIndex = length(v))
{
  getWindow <- function(i) {
    # Don't try to look past the ends of the data
    si <- max(1, i - window)
    ei <- min(length(v), i + window)
    v[si : ei]
  }
  
  maxima <- numeric(length(v) / 2)
  nm <- 0
  for (i in startIndex:endIndex) {
    # Is this point a maximum?
    if (v[i] == max(getWindow(i))) {
      nm <- nm + 1
      maxima[nm] <- i
    }
  }
  
  head(maxima, nm)
}

JPluralise <- function(singular, value, plural = paste0(s, "s")) {
  ifelse (value == 1, singular, plural)
}

# General function to apply a Perl style regex to a string and return the matches
# 
# E.g. 
# > RegexprPerl("(\\d+)[^\\d]*(\\d+)", "1 2")
# [1] "1 2" "1"   "2"
# 
# 
# @param expr - perl style regular expression
# @param str - string to match 
# @value [[$0, $1, ...], [$0, $1, ...]]
RegexprPerl <- function(expr, str) {
  match <- regexpr(expr, str, perl=T)
  matches <- list()

  for (mi in 1:length(match)) {
    if (attr(match, 'match.length')[mi] >= 0) {
      capture_start <- attr(match, 'capture.start')[mi,]
      capture_length <- attr(match, 'capture.length')[mi,]
      total_matches <- 1 + length(capture_start)
      matches[[mi]] <- character(total_matches)
      matches[[mi]][1] <- substr(str, match, match + attr(match, 'match.length') - 1)[[mi]]
      if (length(capture_start) > 1) {
        for (i in 1:length(capture_start)) {
          matches[[mi]][i + 1] <- substr(str, capture_start[[i]], capture_start[[i]] + capture_length[[i]] - 1)[mi]
        }
      }
    } else {
      matches[[mi]] <- character(1 + length(attr(match, 'capture.start')[mi,]))
    }
  }
  
  matches
}


#
# Splining a polygon.
#
#   The rows of 'xy' give coordinates of the boundary vertices, in order.
#   'vertices' is the number of spline vertices to create.
#              (Not all are used: some are clipped from the ends.)
#   'k' is the number of points to wrap around the ends to obtain
#       a smooth periodic spline.
#
#   Returns an array of points. 
# 
# From stackoverflow: 
# http://gis.stackexchange.com/questions/24827/how-to-smooth-the-polygons-in-a-contour-map/24929#24929
# 
# NEEDS SOME WORK!
spline.poly <- function(xy, vertices, k=3, ...) {
  # Assert: xy is an n by 2 matrix with n >= k.
  
  # Wrap k vertices around each end.
  n <- dim(xy)[1]
  if (k >= 1) {
    data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
  } else {
    data <- xy
  }
  
  # Spline the x and y coordinates.
  data.spline <- spline(1:(n+2*k), data[,1], n=vertices, ...)
  x <- data.spline$x
  x1 <- data.spline$y
  x2 <- spline(1:(n+2*k), data[,2], n=vertices, ...)$y
  
  # Retain only the middle part.
  cbind(x1, x2)[k < x & x <= n+k, ]
}

# Print elapsed time. Silent noop if startTime is null.
#
# startTime <- proc.time()
# ...
# ShowTime("Long process:", startTime)
# 
# value - elapsed seconds
ShowTime <- function(msg, startTime) {
  if (!is.null(startTime)) {
    elapsedSecs <- (proc.time() - startTime)[3]
    elapsed <- elapsedSecs
    if (elapsed >= 3600) {
      elapsed <- elapsed / 3600
      units <- "hours"
    } else if (elapsed >= 60) {
      elapsed <- elapsed / 60
      units <- "mins"
    } else {
      units <- "secs"
    }
    cat(sprintf(paste(msg, "%g %s\n"), round(elapsed, 2), units))
    invisible(elapsedSecs)
  }
}

# Used to cache the result of a long running function.
# Eg. cachedata("slow-function-data", function() {...})
#
# Data is serialised to the JDataDir/serialised-data directory
cachedata <- function(fileBasename, fun) {
  path <- JUDFile(paste0(fileBasename, ".rds"), "serialised-data")
  if (file.exists(path)) {
    data <- readRDS(path)
  } else {
    data <- fun()
    saveRDS(data, path)
  }
  data
}

# Fairly ugly hack to allow assignment to several variables at once
# https://stat.ethz.ch/pipermail/r-help/2004-June/053343.html
# Eg. list[,Green,Blue]  <- col2rgb("aquamarine")
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

# Blocks until the current device no longer exists, i.e. a window is closed
PauseUntilWindowIsClosed <- function() {
    while (!is.null(dev.list()))
        Sys.sleep(1)
}

# Returns the specified colour with the same RGB components but the specified alpha
# colour - any values that can be passed to rgb()
# alpha - new transparency value, [0, 255]
#
# Note that some graphics devices may not support partial alpha.
JTransparentColour <- function(colour, alpha) {
  c <- col2rgb(colour)
  rgb(c[1,], c[2,], c[3,], alpha, maxColorValue = 255)
}

# Some "nice" colours
MyPallete <- list(midGrey = "#737373",
               midRed = "#f15a60",
               midGreen = "#7ac36a",
               midBlue = "#5a9bd4",
               midOrange = "#faa75b",
               midPurple = "#9e67ab",
               midBrown = "#ce7058",
               midPink = "#d77fb4",
               darkGrey = "#010202",
               darkRed = "#ee2e2f",
               darkGreen = "#008c48",
               darkBlue = "#185aa9",
               darkOrange = "#f47d23",
               darkPurple = "#662c91",
               darkBrown = "#a21d21",
               darkPink = "#b43894",
               lightGrey = "#cccccc",
               lightRed = "#f2afad",
               lightGreen = "#d9e4aa",
               lightBlue = "#b8d2ec",
               lightOrange = "#f3d1b0",
               lightPurple = "#d5b2d4",
               lightBrown = "#ddb9a9",
               lightPink = "#ebc0da",
               cyan = "#00f0f0"
)

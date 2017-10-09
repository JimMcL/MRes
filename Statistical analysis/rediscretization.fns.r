#

LoadFns("motion")

# Function to "rediscretize" a path to a constant path length.
# Based on the appendix in Bovet and Benhamou, (1988).
#
# points - list of points, with x & y values
# R - rediscretization step length
#
# value is vector of complex points which are the rediscretized path
MRediscretizePoints <- function(points, R) {
  # Note I have attempted 
  
  # Simplify distance calculations by using polar coordinates as implemented in complex
  # p contains the original path points (x, y)
  p <- complex(real = points$x, imaginary = points$y)

  # result will contain the points in discretized path points (X, Y)
  result <- complex(128)
  result[1] <- p[1]
  I <- 1
  j <- 2
  
  while (j <= length(p)) {
    # Find the first point k for which |p[k] - p_0| >= R
    k <- NA
    for (i in j:length(p)) {
      d <- Mod(p[i] - result[I])
      #cat(sprintf("j = %d, i = %d, I = %d, d = %g\n", j, i, I, d))
      if (d >= R) {
        k <- i
        break;
      }
    }
    
    #cat(sprintf("Got k = %d\n", k))
    if (is.na(k)) {
      # We have reached the end of the path
      break
    }
    
    # The next point may lie on the same segment
    j <- k

    # The next point lies on the segment p[k-1], p[k]
    XI <- Re(result[I])
    xk_1 <- Re(p[k - 1])
    YI <- Im(result[I])
    yk_1 <- Im(p[k - 1])
    lambda <- Arg(diff(p[c(k - 1, k)])) #+ ifelse(Re(p[k]) <= xk_1, pi, 0)
    cos_l <- cos(lambda)
    sin_l <- sin(lambda)
    U <- (XI - xk_1) * cos_l + (YI - yk_1) * sin_l
    V <- (YI - yk_1) * cos_l - (XI - xk_1) * sin_l
    
    H <- U + sqrt(abs(R ^ 2 - V ^ 2))
    XIp1 <- H * cos_l + xk_1
    YIp1 <- H * sin_l + yk_1
    
    # This is purely to make the code run (significantly) faster
    if (length(result) < I + 1)
      length(result) <- 2 * length(result)
    
    # Save the point
    result[I + 1] <- complex(real = XIp1, imaginary = YIp1)
    # Move on to next segment
    I <- I + 1
  }

  # Truncate result to actual length
  result <- head(result, I)
  
  result
}

MRediscretizeTrack <- function(track, R) {
  rt <- MRediscretizePoints(track, R)
  # Convert from complex to cartesian coords
  rt <- data.frame(x = Re(rt), y = Im(rt))
  # Fill in other track stuff
  rt <- .MtrackFillIn(rt)
  # Copy attributes across
  for (att in c('numFrames', 'type', 'file')) {
    attr(rt, att) <- attr(track, att)
  }
  rt
}

# Returns angles of each segment relative to the previous segment
MTrackAngles <- function(track) {
  angles <- diff(Arg(track$displacement))
  # Normalise so that -pi < angle <= pi
  ii <- angles <= -pi
  angles[ii] <- angles[ii] + 2 * pi
  ii <- angles > pi
  angles[ii] <- angles[ii] - 2 * pi
  angles
}

test <- function() {
  #track <- .MtrackFillIn(data.frame(x = c(0, 3), y = c(0, 1), Time = c(1:2)))
  # track <- .MtrackFillIn(data.frame(x = c(2, 0), y = c(0, 0), Time = c(1:2)))
  #track <- .MtrackFillIn(data.frame(x = c(0, 3, 6, 9, 9, 5), y = c(0, 0, 0, 1, 2, 1), Time = c(1:6)))
  #s22 <- sqrt(2) / 2
  #track <- .MtrackFillIn(data.frame(x = c(0, 1, 1 + s22, 1 + s22, s22, 1 + s22), y = c(0, 0, s22, 1 + s22, 1 + s22, 2), Time = c(1:6)))
  #np <- 20
   #track <- .MtrackFillIn(data.frame(x = runif(np, 0, 100), y = runif(np, 0, 100), Time = 1:np))
  track <- tracks[[1]]
  
  oldPars <- par(no.readonly = T)
  on.exit(par(oldPars))
  par(mfrow = c(2, 1))
  
  plotTrack(track)
  #points(y ~ x, data = track, col = 'green', pch = 3, cex = .8)
  stepLength <- .5
  rt <- MRediscretizeTrack(track, stepLength)
  lines(y ~ x, data = rt, type = 'l', col = 'red')
  points(y ~ x, data = rt, col = 'blue', pch = 16, cex = .5)
  symbols(rt[nrow(rt),c('x', 'y')], circles = stepLength, add = TRUE, fg = 'red', inches = FALSE)
  plot(density(MTrackAngles(rt)))
}

#test()


# Function to produce plots for my MRes poster.

ResultsPlotForPoster <- function(angle, ms, averageTypes, pca) {
  posterDir <- JUSecureFile("Classes\\MRes poster")
  
  # Results plot
  PlotToPng(file.path(posterDir, sprintf('Constraints result %s.png', angle)),
            function () {
              par(mar = c(3, 3, 4, 2))
              plotVariation(ms, averageTypes, simple = TRUE, 
                            includeTypes = c('spider mimic strict', 'insect mimic strict'), 
                            lwd = 8, title = 'Constraints Result', cex.main = 2, 
                            leg.cex = 1.4, legBg = 'white')
              title(xlab = 'Mimetic imperfection', line = 0.8, cex.lab = 1.6)
              title(ylab = 'Actual Density', line = 0.8, cex.lab = 1.6)
            },
            width = 2800, height = 2800, res = 400
  )
  
  # Morphospace plot. Hardwire colours to get what we want
  palette <- function (nColours) {c(MyPallete$darkRed, MyPallete$midPurple, MyPallete$midGreen, MyPallete$darkBlue)}
  legIt <- function(colours) {
    classes <- c("Models", "Spider mimics", "Insect mimics", "Non-mimics")
    # I don't know how colours are assigned to levels in pca$fac$type, so hard-wire them
    legColours <- colours[c(2, 4, 1, 3)]
    # Work around a bug - the top border of the legend is sometimes clipped, even though it's actually inside the plot area
    # Unfortunately, this messes up the coordinate system so it's not possible to add to the plot after this function exits!
    oldPars <- par(no.readonly = T)
    on.exit(par(oldPars))
    par(xpd=TRUE)
    legend("topleft", 
           legend=c(classes, "Actual outline", "Apparent outline"), 
           cex = 2,
           col = c(legColours, 'black', 'black'),
           pch = c(17, 17, 17, 17, 17, 16),
           inset = c(.03, .04))
  }
  PlotToPng(file.path(posterDir, paste0('morpho-', angle, '.png')),
            function () {
              plotPca(pca, labels = FALSE, cex = 3, 
                      labels.col = 'black', labels.cex = 3, ellipses = FALSE,
                      legend = FALSE, extraPlottingFn = legIt, palette = palette) 
            },
            width = 1200, height = 800)
  
  ###
  # Medium-complexity plots plot for handout
  handoutDir <- file.path(posterDir, 'handout')
  PlotToPng(file.path(handoutDir, sprintf('Constraints result %s.png', angle)),
            function () {
              par(mar = c(3, 3, 2, 2))
              plotVariation(ms, averageTypes, simple = TRUE, 
                            includeTypes = c('model strict', 'spider mimic strict', 'insect mimic strict', 'non-mimic strict'), 
                            lwd = 5, title = '', 
                            leg.cex = 1.5, legBg = 'white')
              title(xlab = 'Mimetic imperfection', line = 0.8, cex.lab = 1.6)
              title(ylab = 'Density', line = 0.8, cex.lab = 1.6)
            }
            
  )
  palette <- function (nColours) { c(MyPallete$darkBlue, MyPallete$midOrange, MyPallete$midGreen, MyPallete$darkRed) }
  legIt <- function(colours) {
    classes <- c("Models", "Spider mimics", "Insect mimics", "Non-mimics")
    # I don't know how colours are assigned to levels in pca$fac$type, so hard-wire them
    legColours <- colours[c(2, 4, 1, 3)]
    # Work around a bug - the top border of the legend is sometimes clipped, even though it's actually inside the plot area
    # Unfortunately, this messes up the coordinate system so it's not possible to add to the plot after this function exits!
    oldPars <- par(no.readonly = T)
    on.exit(par(oldPars))
    par(xpd=TRUE)
    legend("topleft", 
           legend=c(classes, "Actual outline", "Apparent outline"), 
           cex = 1.5,
           col = c(legColours, 'black', 'black'),
           pch = c(17, 17, 17, 17, 17, 16),
           inset = c(.03, .04))
  }
  PlotToPng(file.path(handoutDir, paste0('morpho-', angle, '.png')),
            function () {
              plotPca(pca, labels = 'scientificName', 
                      labels.col = NULL, labels.cex = 1,
                      legend = FALSE, extraPlottingFn = legIt, palette = palette) 
            },
            width = 1200, height = 600)
  
}


LoadFns("spectral_analysis")

# Location of files generated from JAZ spectrometer
originalsDir <- JUUniFile(subDir = "Classes/Thesis/Spectra")

by <- 5 # No. of copies of each measurement

# Explore spectra in the usual directory
# Imports all the spectral files in originalsDir, and (optionally) plots raw spectra in groups of "by"
s <- ExploreDir(originalsDir, by = by)

# Combine spectra for repeated readings
spec <- aggspec(s, by = by, FUN = "mean")
# Smooth
spec.sm <- procspec(spec, opt = c("smooth"), span = .12)
plot(spec.sm, col = adjustColour(col2rgb(spec2rgb(spec.sm)), 1, 2, 8))

# Scale
spec.max <- procspec(spec.sm, opt = c('max'))
plot(spec.max, col = adjustColour(col2rgb(spec2rgb(spec.max)), 1, 2, 8))

par(mfrow=c(1,2))
plot(spec.sm, col = adjustColour(col2rgb(spec2rgb(spec.sm)), 1, 2, 8))
plot(spec.max, col = adjustColour(col2rgb(spec2rgb(spec.sm)), 1, 2, 8))

########################################
par(mfrow=c(1,1))

# Bee
vm <- vismodel(spec.sm, visual = 'apis', achromatic = 'none', relative = T)
cs <- colspace(vm, space = 'tri')
plot(cs)

# Bird
vm <- vismodel(spec.sm, visual = 'avg.uv', achromatic = 'none', relative = T)
cs <- colspace(vm, space = 'tcs')
plot(cs)
tcsplot(cs)


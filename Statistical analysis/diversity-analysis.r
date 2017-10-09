
# Plot species diversity
LoadFns("diversity")

THESIS_DIR <- JUCFile(subDir = "Thesis")
ALA_DATA_PATH <- file.path(THESIS_DIR, "Downloaded data/ALA-Araneaea-records-2017-01-09")


################################################################################

# ALA records
spiderFile <- file.path(ALA_DATA_PATH, "records-2017-01-09.csv")
spiders <- ALAReadOccData(spiderFile)
source <- 'ALA'
# Only keep spiders identified to genus or better
spiders <- spiders[spiders$taxonRank %in% c("genus", "species", "subspecies"),]

# GBIF records
# spiders <- ReadGBIFCsvData(file.path(THESIS_DIR, "Downloaded data/GBIF/0056059-160910150852091.csv"))
# spiders <- spiders[spiders$taxonrank %in% c("GENUS", "SPECIES", "SUBSPECIES"),]
# source <- 'GBIF'

# Only keep records on continental Aus
spiders.aus <- CropToAustralianMainland(spiders, T)
# Project to albers
spiders.proj <- OCCProject(spiders.aus, CRS(albers.proj4))

res <- 1
colname = "genus"
humanName = "genera"
PlotDiversityToPng(spiders.proj, res, "Spider", paste(source, "spider", sep = '-'), colname, humanName)

mimicGenera <- c("Micaria", "Myrmarachne", "Amyciaea", "Storena", "Tinytrema", "Judalana", "Ligonipes", "Rhombonotus", "Damoetas", "Myrmarachne", "Castianeira", "Apochinomma")
mimicOccs <- spiders.proj[spiders.proj$genus %in% mimicGenera,]
PlotDiversityToPng(mimicOccs, res, "Mimic", paste(source, "mimic", sep = '-'), colname, humanName = humanName)


######
# Maxent analysis
latLonMimics <- spiders.aus[spiders.aus$genus %in% mimicGenera,]
mdf <- data.frame(species=latLonMimics$genus, longitude = coordinates(latLonMimics)[,1], latitude = coordinates(latLonMimics)[,2])
write.csv(mdf[order(mdf$species),], file.path(ALA_DATA_PATH, "mimics-for-maxent.csv"), row.names = F)
sort(table(mdf$species), decreasing = T)

## Sum up all of the results for each species/genus (whatever is in mdf)
genera <- unique(mdf$species)
generaToPlot <- genera[genera != 'Storena']
rasters <- sapply(generaToPlot, function(sp) {
  print(sp)
  file <- file.path(THESIS_DIR, "Maxent output", "Mimics by genus", paste0(sp, ".asc"))
  print(file)
  raster(readAsciiGrid(file))
})

PlotToPng(file.path(THESIS_DIR, "Diversity", "mimics by maxent.png"),
          function() {
            plot(Reduce("+", rasters))
            PlotAusCities(text = T)
            text(110, -42, paste(generaToPlot, collapse = ", "), adj = 0, cex = .7)
          })

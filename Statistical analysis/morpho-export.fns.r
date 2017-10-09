# Functions to export mimetic imperfection values to CSV files

suppressMessages(library(plyr))
LoadFns("sample-dbs")

# Calculates an imperfection value for every photo
# typeShapes - averaged shapes for all mimic types - used to obtain the average shape of an ant
# shapes - shapes for each photo
MpPCACalcImperfectionPerShape <- function(shapesCoe, typeShapes) {
  mi <- which(typeShapes$Coe$fac$type == MTP_MODEL & typeShapes$Coe$fac$strictOrLax == 'strict')
  dists <- CalcMimeticVariation(typeShapes$Coe$coe, shapesCoe$coe, rep(mi, length(shapesCoe$fac$type)))
  # Record the id of the source photo, not the outline photo
  ids <- SourcePhotoId(shapesCoe$fac$source)
  data.frame(id = ids, 
             mimicType = shapesCoe$fac$type, 
             species = shapesCoe$fac$scientificName, 
             bodyLength = shapesCoe$fac$bodylength,
             imperfection = ScaleToRange(dists),
             angle = shapesCoe$fac$angle,
             disposition = shapesCoe$fac$disposition) 
}


# Calculates "mimetic imperfection" for individuals with a common angle
MpPCACalcIndividualImperfection <- function(individualsShapes, typeShapes) {
  # Calculate average model shape
  mi <- which(typeShapes$Coe$fac$type == MTP_MODEL)
  
  dists <- CalcMimeticVariation(typeShapes$Coe$coe, individualsShapes$Coe$coe, rep(mi, length(individualsShapes$Coe$fac$type)))
  data.frame(id = individualsShapes$Coe$fac$imageableid,
             mimicType = individualsShapes$Coe$fac$type, 
             species = individualsShapes$Coe$fac$scientificName, 
             bodyLength = individualsShapes$Coe$fac$bodylength,
             imperfection = ScaleToRange(dists),
             angle = individualsShapes$Coe$fac$angle,
             disposition = individualsShapes$Coe$fac$disposition) 
}

# Calculates "mimetic imperfection" for species with a common angle
MpPCACalcSpeciesImperfection <- function(speciesShapes, typeShapes) {
  # Calculate average model shape
  mi <- which(typeShapes$Coe$fac$type == MTP_MODEL)
  
  dists <- CalcMimeticVariation(typeShapes$Coe$coe, speciesShapes$Coe$coe, rep(mi, length(speciesShapes$Coe$fac$type)))
  data.frame(mimicType = speciesShapes$Coe$fac$type, 
             species = speciesShapes$Coe$fac$scientificName, 
             angle = speciesShapes$Coe$fac$angle,
             imperfection = dists)
}

# Export CSV files with "mimetic imperfection" for individuals and species
MExportVariation <- function(dorsalData, lateralData, dir) {
  
  # For photos, combine dorsal and lateral
  cols <- c('type', 'species', 'bodylength', 'disposition', 'accuracy')
  d <- rbind(
      cbind(data.frame(id = SourcePhotoId(dorsalData$photo$fac$source),  angle = 'dorsal'),   dorsalData$photo$fac[,cols]),
      cbind(data.frame(id = SourcePhotoId(lateralData$photo$fac$source), angle = 'lateral'), lateralData$photo$fac[,cols]))
  write.csv(d[order(d$accuracy, decreasing = TRUE),], file = file.path(dir, "accuracy-per-photo-morpho.csv"), row.names = FALSE)

  # For individuals, combine dorsal and lateral
  cols <- c('imageableid', 'type', 'species', 'bodylength', 'disposition', 'accuracy')
  d <- rbind(
    cbind(angle = 'dorsal',    dorsalData$individual$Coe$fac[,cols]),
    cbind(angle = 'lateral', lateralData$individual$Coe$fac[,cols]))
  write.csv(d[order(d$accuracy, decreasing = TRUE),], file = file.path(dir, "accuracy-per-individual-angle-morpho.csv"), row.names = FALSE)
  
  # For species, combine dorsal and lateral
  cols <- c('type', 'species', 'bodylength', 'accuracy')
  d <- rbind(
    cbind(angle = 'dorsal',    dorsalData$species$Coe$fac[,cols]),
    cbind(angle = 'lateral', lateralData$species$Coe$fac[,cols]))
  write.csv(d[order(d$accuracy, decreasing = TRUE),], file = file.path(dir, "accuracy-per-species-angle-morpho.csv"), row.names = FALSE)
}

# UpdateImperfectionCSV <- function(dir, averageTypes, photoShapes) {
#   # Calculate imperfection
#   new <- CalcImperfectionPerShape(photoShapes, averageTypes)
#   
#   # Incorporate any pre-calculated photos not in this batch  
#   photoImperfectionCsv <- file.path(dir, 'imperfection-per-photo-morpho.csv')
#   if (file.exists(photoImperfectionCsv)) {
#     old <- read.csv(photoImperfectionCsv)
#     
#     # Append any old values which are missing from new (i.e. probably a different angle)
#     missingIds <- old$Id[!(old$Id %in% new$Id)]
#     if (length(missingIds) > 0) {
#       missingIdxs <- which(old$Id %in% missingIds)
#       new[(nrow(new) + 1):(nrow(new) + length(missingIdxs)),] <- old[missingIdxs,]
#     }
#   }
#   
#   # Write out the result
#   write.csv(new, photoImperfectionCsv, row.names = FALSE)
# }


MExportAppendixIndividuals <- function(dorsalData, lateralData, dir, basename) {
  # Take a subset of the columns because we only want 1 row per specimen, not per photo
  mc <- c("imageableid", "type", "order", "scientificName", 
          "scientificNameAuthorship", "decimalLatitude", 
          "decimalLongitude", "coordinateUncertaintyInMeters", "elevation", 
          "locationRemarks", "time", "day", "month", "year", "bodylength",
          "disposition")
  # Merge the 2 datasets, just keeping info about individual specimens
  d <- merge(dorsalData$individual$Coe$fac[,mc], lateralData$individual$Coe$fac[,mc], all = TRUE, suffixes = c('Dorsal', 'Lateral'))
  
  # Incorporate mimetic accuracy
  .setAccuracyCol <- function(d, data, col, label) {
    d[, label] <- NA
    fac <- data$individual$Coe$fac
    d[match(fac$imageableid, d$imageableid), label] <- fac[,col]
    d
  }
  
  # LDA accuracy
  d <- .setAccuracyCol(d, dorsalData, 'accuracy', 'dorsalAccuracyMorphometric')
  d <- .setAccuracyCol(d, lateralData, 'accuracy', 'lateralAccuracyMorphometric')
  # Google vision accuracy
  d <- .setAccuracyCol(d, dorsalData, 'gvAccuracy', 'dorsalAccuracyMachineLearning')
  d <- .setAccuracyCol(d, lateralData, 'gvAccuracy', 'lateralAccuracyMachineLearning')
  
  # Change name of column imageableid to catalogNumber a la Darwin core,
  # and change disposition to source (also change its value)
  d <- plyr::rename(d, c("imageableid" = 'catalogNumber', 'disposition' = 'source'))
  # Identify where I got the specimen from, i.e. credit the museum for borrowed specimens
  d$source <- ifelse(grepl('on loan', d$source, ignore.case = TRUE), d$source, 'Collected')
  d$source <- sub('on loan', 'Borrowed', d$source, ignore.case = TRUE)
  
  # Change some of the location remarks which are only meaningful to me
  d$locationRemarks <- ifelse(grepl('^PF', d$locationRemarks), 'MQ campus', d$locationRemarks)
  
  write.csv(d[order(d$scientificName),], file = file.path(dir, paste0(basename, ".csv")), row.names = FALSE, na = '')
  cat(sprintf("%d outlines, %d dorsal, %d lateral\n", nrow(dorsalData$photo$fac) + nrow(lateralData$photo$fac), nrow(dorsalData$photo$fac), nrow(lateralData$photo$fac)))
  cat(sprintf("%d individual shapes, %d dorsal, %d lateral\n", nrow(d), sum(!is.na(d$dorsalAccuracyMorphometric)), sum(!is.na(d$lateralAccuracyMorphometric))))
}

MExportAppendixSpecies <- function(dorsalData, lateralData, dir, basename) {
  dd <- dorsalData$species$Coe$fac
  ld <- lateralData$species$Coe$fac
  
  # Incorporate references for decision to count as a mimic
  dd <- cbind(dd, .MGetMimicryDecisions(dd))
  ld <- cbind(ld, .MGetMimicryDecisions(ld))
  
  # Subset of the columns which are meaningful for species
  mc <- c("type", "order", "scientificName", "scientificNameAuthorship", 'mimicryReference', 'mimicryNotes')
  # Merge the 2 datasets, just keeping info about species
  d <- merge(dd[,mc], ld[,mc], all = TRUE, suffixes = c('Dorsal', 'Lateral'))

  # Incorporate mimetic accuracy
  .setAccuracyCol <- function(d, data, col, label) {
    d[, label] <- NA
    fac <- data$species$Coe$fac
    for(i in 1:nrow(d)) {
      row <- d[i,]
      facRow <- which(fac$type == row$type & fac$scientificName == row$scientificName)[1]
      d[i, label] <- fac[facRow, col]
    }
    d
  }
  # LDA accuracy
  d <- .setAccuracyCol(d, dorsalData, 'accuracy', 'dorsalAccuracyMorphometric')
  d <- .setAccuracyCol(d, lateralData, 'accuracy', 'lateralAccuracyMorphometric')
  # Google vision accuracy
  d <- .setAccuracyCol(d, dorsalData, 'gvAccuracy', 'dorsalAccuracyMachineLearning')
  d <- .setAccuracyCol(d, lateralData, 'gvAccuracy', 'lateralAccuracyMachineLearning')
  
  # Body length in dorsalData and lateralData is only calculated from specimens with dorsal or lateral outlines.
  # Instead, calculate from all specimens with known length
  sp <- DbQuerySpecimens(sprintf("sql=taxon_id in (select id from taxa where scientific_name in ('%s'))", paste(unique(d$scientificName), collapse="','")))
  sp <- sp[!is.na(sp$bodylength),]
  sm <- aggregate(sp$bodylength, list(mimicType(sp), sp$scientificName, sp$scientificNameAuthorship), mean)
  names(sm) <- c('type', 'scientificName', 'scientificNameAuthorship', 'bodylength')
  d <- merge(d, sm, all.x = TRUE)
  
  write.csv(d[order(d$scientificName),], file = file.path(dir, paste0(basename, ".csv")), row.names = FALSE, na = '')
  cat(sprintf("%d species shapes, %d dorsal, %d lateral\n", nrow(d), sum(!is.na(d$dorsalAccuracyMorphometric)), sum(!is.na(d$lateralAccuracyMorphometric))))
}

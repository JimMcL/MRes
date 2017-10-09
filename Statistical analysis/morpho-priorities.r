#!Rscript
# 
# Determines priorities specimens for outline conversion
LoadFns("sample-dbs")
LoadFns("morpho-calculation")

# Would like this many individuals per species
GOAL <- 5

SummariseData <- function(con = stdout()) {
  
  occCnt <- function(species) {
    t <- table(species)
    c(goal = sum(t >= 5),
      four = sum(t == 4),
      three = sum(t == 3),
      two = sum(t == 2),
      one = sum(t == 1),
      total = sum(t >= 1))
  }
  
  qrySpecimensWithPhotos <- function(photoQuery) {
    photos <- DbQueryPhotos(photoQuery)
    specimens <- DbQuerySpecimens('id=[', paste(unique(photos$imageableid), collapse=','), ']')
    specimens <- unique(specimens)
    specimens$mimicType <- mimicType(specimens)
    specimens
  }

  # With dorsal photos  
  specimensWithPhotos <- qrySpecimensWithPhotos('ptype=Photo&imageable_type=Specimen&view_angle=dorsal')
  
  # With dorsal outlines  
  specimensWithOutlines <- qrySpecimensWithPhotos('ptype=Outline&imageable_type=Specimen&view_angle=dorsal')
  
  # With motion videos
  specimensWithVideos <- qrySpecimensWithPhotos('ptype=Motion&imageable_type=Specimen')
  
  # Alive
  alive <- DbQuerySpecimens('q=alive')
  # Just in case the word "alive" appears in another field
  alive <- alive[grepl("alive", alive$disposition, ignore.case = TRUE),]
  alive$mimicType <- mimicType(alive)
  
  
  write.csv(rbind(
    spider.mimics.with.photos = occCnt(specimensWithPhotos[specimensWithPhotos$mimicType == MTP_SPIDER,]$species), 
    insect.mimics.with.photos = occCnt(specimensWithPhotos[specimensWithPhotos$mimicType == MTP_INSECT,]$species), 
    all.with.photos = occCnt(specimensWithPhotos$species), 
    spider.mimics.with.outlines = occCnt(specimensWithOutlines[specimensWithOutlines$mimicType == MTP_SPIDER,]$species), 
    insect.mimics.with.outlines = occCnt(specimensWithOutlines[specimensWithOutlines$mimicType == MTP_INSECT,]$species), 
    all.with.outlines = occCnt(specimensWithOutlines$species), 
    spider.mimics.with.videos = occCnt(specimensWithVideos[specimensWithVideos$mimicType == MTP_SPIDER,]$species),
    insect.mimics.with.videos = occCnt(specimensWithVideos[specimensWithVideos$mimicType == MTP_INSECT,]$species),
    all.with.videos = occCnt(specimensWithVideos$species),
    spider.mimics.alive = occCnt(alive[alive$mimicType == MTP_SPIDER,]$species),
    insect.mimics.alive = occCnt(alive[alive$mimicType == MTP_INSECT,]$species),
    ants.alive = occCnt(alive[alive$mimicType == MTP_MODEL,]$species),
    all.alive = occCnt(alive$species)),
    file = con)
}

ListOutlinePriorities <- function(angle, con = stdout()) {
  photos <- DbQueryPhotos('ptype=Photo&imageable_type=Specimen&view_angle=', angle)
  specimensWithPhotos <- DbQuerySpecimens('id=[', paste(unique(photos$imageableid), collapse=','), ']')
  specimensWithPhotos <- unique(specimensWithPhotos)
  specimensWithPhotos$maxRating <- sapply(specimensWithPhotos$id, function(sid) max(c(0,photos[photos$imageableid == sid,]$rating), na.rm = TRUE))
  specimenTypes <- mimicType(specimensWithPhotos)
  specimensWithPhotos$mimicType <- specimenTypes
  
  outlines <- DbQueryPhotos('ptype=Outline&imageable_type=Specimen&view_angle=', angle)
  # Only want strict outlines  
  outlines <- outlines[strictOrLax(outlines) == 'strict',]
  
  outlineIdxs <- na.omit(match(unique(outlines$imageableid), specimensWithPhotos$id))
  specimensWithOutlines <- specimensWithPhotos[outlineIdxs,]
  specimensWithoutOutlines <- specimensWithPhotos[-outlineIdxs,]
  
  # Want specimens with photos not yet outlined, < 5 outlined for species, ordered by max photo rating
  # Priority species
  priSpecimens <- specimensWithPhotos[specimenTypes %in% c(MTP_SPIDER, MTP_INSECT),]
  t <- as.data.frame(table(priSpecimens$species))
  colnames(t) <- c('species', 'individuals')
  t$processed <- sapply(t$species, function(sp) sum(specimensWithOutlines$species == sp))
  t$shortfall <- GOAL - t$processed
  t$maxRating <- sapply(t$species, function(sp) max(c(0,specimensWithPhotos[specimensWithoutOutlines$species == sp,]$maxRating)))
  # Eliminate species which already have >= 5 processed specimens
  t <- t[t$processed < GOAL,]
  # Sort on no. of specimens and max rating
  t <- t[order(t$individuals, t$maxRating, decreasing = TRUE),]
  
  invisible(apply(t, 1, function(row) {
    sp <- row['species']
    sf <- as.numeric(row['shortfall'])
    specs <- specimensWithoutOutlines[specimensWithoutOutlines$species == sp,]
    specs <- specs[order(specs$maxRating, decreasing = TRUE),]
    specs <- head(specs, sf)
    writeLines(sprintf('Specimen %d %s', specs$id, specs$species), con = con)
  }))
}

ListMotionPriorities <- function(con = stdout()) {
  photos <- DbQueryPhotos('ptype=Motion&imageable_type=Specimen')
  only900mm <- TRUE
  if (only900mm)
    photos <- photos[grepl('900mm', photos$description),]
  specimensWithVideos <- DbQuerySpecimens('id=[', paste(unique(photos$imageableid), collapse=','), ']')
  specimensWithVideos <- unique(specimensWithVideos)

  # Get currently alive specimens
  alive <- DbQuerySpecimens('q=alive')
  # Just in case the word "alive" appears in another field
  alive <- alive[grepl("alive", alive$disposition, ignore.case = TRUE),]
  alive$mimicType <- mimicType(alive)
  
  # Eliminate all species which already have 5 or more motion videos
  t <- table(alive$species)
  goalSpecies <- names(t)[t < GOAL & t >= 1]
  alive <- alive[alive$species %in% goalSpecies,]
  # For now, don't worry about non-mimics
  alive <- alive[alive$mimicType != MTP_NON_MIMIC,]

  # Now list remaining specimens without a motion video
  missingIds <- setdiff(alive$id, specimensWithVideos$id)
  
  for (i in 1:length(missingIds)) {
    id <- missingIds[i]
    spec <- alive[alive$id == id,]
    writeLines(sprintf('%d: Specimen %s, %s - %s', i, spec$id, spec$species, spec$mimicType), con = con)
  }
}


.usage <- function() {
  writeLines("Usage : {motion|outlines|summary}", con = stderr())
}

args <- commandArgs(TRUE)
if (length(args) != 1) {
  .usage()
} else {
  cmd <- args[1]
  if (cmd == 'motion')
    ListMotionPriorities()
  else if (cmd == 'outlines')
    ListOutlinePriorities('dorsal')
  else if (cmd == 'summary')
    SummariseData()
  else
    .usage()
}


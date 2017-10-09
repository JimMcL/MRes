# Functions to query the samples database.
# The database contains information about specimens and their photos.
# 
# The database is queried using an HTTP interface.
# The host defaults to localhost, but can be specified on the R command line as 
# Rscript ... -h <host>
# or by setting the global variable MORPHO_HOST.
#
LoadFns("download")
LoadFns('maps')
library(ggmap)


# Returns the name of the host for database queries.
# Check the R command line for an option like '-host <host>'
# otherwise if MORPHO_HOST is set, returns it, otherwise returns 'localhost'.
DbChooseHost <- function(args = commandArgs(TRUE)) {
  host <- JCLIValue('h', 'host')
  if (!is.null(host))
    host
  else if (exists('MORPHO_HOST'))
    MORPHO_HOST
  else
    'localhost'
}

# Returns a URL by pasting arguments on the end of 'http://MORPHO_HOST/'
DbUrl <- function(...) {
  u <- sprintf("http://%s/%s", DbChooseHost(), paste0(...))
  #cat(paste0(u, '\n'))
  u
}

# Returns data from the database
DbQuery <- function(...) {
  read.csv(DbUrl(...), stringsAsFactors = F, strip.white=TRUE)
}

# Returns photo metadata from the database. Any arguments are used to construct the URL
DbQueryPhotos <- function(...) {
  DbQuery('photos.csv?', ...)
}

# Returns specimen metadata from the database. Any arguments are used to construct the URL
DbQuerySpecimens <- function(...) {
  DbQuery('specimens.csv?', ...)
}

# Returns specimen metadata corresponding to a set of photos
DbQuerySpecimensForPhotos <- function(photos) {
  specimenIds <- unique(photos[photos$imageabletype == 'Specimen',]$imageableid)
  DbQuerySpecimens(sprintf("id=[%s]", paste(specimenIds, collapse=',')))
}

# Returns a data frame containing details of all photos in the database which satisfy the passed in query.
# The photo image files are downloaded to a local cache, and the file name is stored in the column "file".
#
# ... - arguments passed to DbQueryPhotos
# tempfileFnFn - optional function which returns a function passed to JDownload as the tempfileFnFn argument. Called with 1 argument, the photos data frame.
# strictOnly - if TRUE (the default), only strict outlines are returned (only meaningful for outlines)
# NOTE: only returns strict outlines
DbGetPhotoData <- function(..., dir = JCacheDir, tempfileFnFn = NULL, strictOnly = FALSE, debug = FALSE) {
  
  # Get the list of matching photos
  photos <- DbQueryPhotos(...)

  # For now, don't analyse lax outlines - need to put more thought into how to create an apparent/perceived outline
  if (strictOnly) {
    sOrL <- strictOrLax(photos)
    photos <- photos[sOrL == 'strict',]
  }
  
  # Allow caller to specify how to name downloaded files
  tempfileFn <- NULL
  if (!is.null(tempfileFnFn))
    tempfileFn <- tempfileFnFn(photos)
  
  # Download the photos
  file <- JDownload(photos$url, verbose = F, cacheDir = dir, tempfileFn = tempfileFn, debug = debug)
  
  ## Read specimen info for each photo, and add to each photo
  specimens <- DbQuerySpecimensForPhotos(photos)
  # Some entries in the dbs have multiple specimens, which we don't care about, so eliminate them
  specimens <- unique(specimens)
  # Match specimens to photos
  fac <- specimens[match(photos$imageableid, specimens$id),]
  # Determine mimc type before removing description column
  mt <- mimicType(fac)
  # Remove columns with duplicate names in both photos and specimens
  fac <- fac[,!colnames(fac) %in% c('id', 'description')]

  sOrL <- strictOrLax(photos)
  
  ## Build photos data frame containing info about photos and specimens
  cbind(photos, 
        file, 
        type = mt, 
        label = sprintf("%s, %s", photos$imageableid, photos$id), 
        individualId = individualIds(photos, sOrL, mt),
        speciesId = speciesIds(fac, sOrL, mt),
        typeId = typeIds(sOrL, mt),
        strictOrLax = sOrL,
        fac,
        stringsAsFactors = F)
}

# Returns a data frame containing details of all photos in the database with the specified viewing angle.
# The photo image files are downloaded to a local cache, and the file name is stored in the column "file".
# angle - 'dorsal' or 'lateral'
DbGetOutlines <- function(angle) {
  DbGetPhotoData("ptype=Outline&imageable_type=Specimen&view_angle=", angle, strictOnly = TRUE)
}

# Plots a map of the specimens' collection sites
DbMapSpecimenSites <- function(specimens, xlimFrac = 0.05, ylimFrac = 0.05, showLegend = FALSE) {
  sites <- specimens[specimens$type == MTP_SPIDER,c('scientificName', 'decimalLatitude', 'decimalLongitude')]
  xlim <- extendrange(sites$decimalLongitude, f = xlimFrac)
  ylim <- extendrange(sites$decimalLatitude, f = ylimFrac)
  MapSpeciesOcc(sites, columnNames = c("scientificName", "decimalLongitude", "decimalLatitude"), xlim = xlim, ylim = ylim, scaleFactor = 5, showLegend = showLegend)
  #PlotAusCities(text = TRUE)
  map.cities(minpop = 50000)
}

# Plots a google map of the specimens' collection sites
GMapSpecimenSites <- function(specimens, xlimFrac = 0.5, ylimFrac = 0.5, showLegend = FALSE) {

  cols <- c('decimalLatitude', 'decimalLongitude')
  sites <- ddply(specimens, cols, nrow)
  names(sites) <- c(cols, 'Count')
  
  xlim <- extendrange(sites$decimalLongitude, f = xlimFrac)
  ylim <- extendrange(sites$decimalLatitude, f = ylimFrac)
  m <- get_map(c(xlim[1], ylim[1], xlim[2], ylim[2]))

  # Explicitly print so that it works when plotting to a device  
  print(ggmap(m) + 
          geom_point(data = sites, aes(x = decimalLongitude, y = decimalLatitude, size = Count), 
                     color = 'black', fill= 'red', pch = 21, alpha = 0.5) + 
          scale_size_continuous(range = c(5, 20)) +
          theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))
}


# Mimicry specific functions #############################################################
# Following are functions which are specific to my usage of the sampleIt database

MTP_MODEL <- 'model'
MTP_SPIDER <- 'mimetic spider'
MTP_INSECT <- 'mimetic insect'
MTP_NON_MIMIC <- 'non-mimic'
MTP_NAMES <- c(MTP_MODEL, MTP_SPIDER, MTP_INSECT, MTP_NON_MIMIC)
MTP_MIMICS <- c(MTP_SPIDER, MTP_INSECT)

# Extracts year, month and day, and optional hour and minute, columns 
# from the specified data frame, and converts them to a date
DBDate <- function(x) {
  if('hour' %in% names(x))
    ISOdate(x$year, x$month, x$day, x$hour, x$minute, 0)
  else
    ISOdate(x$year, x$month, x$day)
}

# Determines whether a specimen is a mimic, model or non-mimic,
# based on columns family, order, class, description and other.
# Mimics are divided into spider mimics and insect mimics.
# Ants are models, mimics are non-ants with description which includes 'mimic',
# everything else is a non-mimic
mimicType <- function(fac) {
  ants <- fac$family == 'Formicidae'
  spiders <- fac$order == 'Araneae'
  insects <- fac$class == 'Insecta'
  mimics <- (grepl("mimic", fac$description, ignore.case = T) | 
               grepl("mimic", fac$other, ignore.case = T)) & 
    !(grepl("non", fac$description, ignore.case = T) | 
        grepl("non", fac$other, ignore.case = T))
  
  # EXPERIMENT TODO SCRAP THIS CODE!!! try excluding Cosmophasis (except bitaeniata) from mimics
  mimics <- mimics & !(fac$genus == 'Cosmophasis' & fac$species != 'Cosmophasis bitaeniata')
  
  as.factor(ifelse (ants, MTP_MODEL, 
                    ifelse(mimics & spiders, MTP_SPIDER, 
                           ifelse(mimics & insects, MTP_INSECT, 
                                  MTP_NON_MIMIC))))
}

# Reads in a CSV file which contains references and notes for why specimens are 
# classed as mimics. 
# Returns a data frame with 2 columns: 'mimicryReference', 'mimicryNotes'
.MGetMimicryDecisions <- function(d) {
  refs <- read.csv(JUSecureFile("mimicry refs.csv", "Classes/MRes Thesis/Traits"), stringsAsFactors = FALSE, na.strings = '')
  # Cleanup whitespace
  for(col in c('species', 'genus', 'family'))
    refs[,col] <- JTrim(refs[,col])
  
  ri <- unlist(JApplyToRows(d, function(r) {
    m <- c()
    if (grepl('mimetic.*', r$type)) {
      # Choose the most specific match
      m <- which(refs$species == r$species)
      if (length(m) == 0 && nchar(r$genus) > 0)
        m <- which(refs$genus == r$genus)
      if (length(m) == 0 && nchar(r$family) > 0)
        m <- which(refs$family == r$family)
    }
    ifelse(length(m) == 0, NA, m)
  }))
  
  refs[ri,c('mimicryReference', 'mimicryNotes')]
}



# Construct a new mimic type factor for plotting. This does 2 things:
# 1. removes any unused level (e.g. "insect mimic")
# 2. orders the factor in the way I want it
# Normally these things don't matter, but they do matter when plotting 
AsPlottableMimicTypeFactor <- function(types) {
  # Order of MTP_NAMES is important here
  levels <- MTP_NAMES[MTP_NAMES %in% types]
  factor(types, levels = levels)
}

.typeToVal <- function(type, val, unknownTypeVal) {
  # I was relying on numeric factor values to do this, but ran into mysterious cases 
  # when it didn't produce the correct results! This mechanism is more general
  idx <- match(type, MTP_NAMES)
  ifelse(is.na(idx), unknownTypeVal, val[idx])
}

# Given a type (e.g. 'model', 'spider mimic', ...), returns a colour to use when plotting the type.
# Works using character comparisons, so factor values numeric values are not assumed.
typeToCol <- function(type) {
  .typeToVal(type, c( MyPallete$darkPurple, MyPallete$darkRed, MyPallete$darkBlue, MyPallete$darkGreen), 'black')
}

# Given a type (e.g. 'model', 'spider mimic', ...), returns a line style to use when plotting the type
# Works using character comparisons, so factor values numeric values are not assumed.
typeToLty <- function(type) {
  .typeToVal(type, c(4, 1, 2, 3), 6)
}

# Given a type (e.g. 'model', 'spider mimic', ...), returns a line style to use when plotting the type
# Works using character comparisons, so factor values numeric values are not assumed.
typeToPch <- function(type) {
  #.typeToVal(type, c(2, 15, 16, 6), 43)
  .typeToVal(type, c(24, 15, 16, 25), 43)
}

typeToLabel <- function(type, plural = TRUE) {
  pl <- c('Models', 'Mimetic spiders', 'Mimetic insects', 'Non-mimics')
  sing <- c('Model', 'Mimetic spider', 'Mimetic insect', 'Non-mimic')
  .typeToVal(type, if(plural) pl else sing, type)
}

# Returns 'strict' or 'lax' as a factor
strictOrLax <- function(photos) {
  as.factor(ifelse(grepl('strict', photos$description), 'strict', 'lax'))
}

# Returns individual ids as a factor.
# For the purpose of plotting, an individual has a unique imageableid 
# and either strict or lax photos, + mimic status in order to differentiate between 
# e.g. nymphs which are mimics and adults which are not
individualIds <- function(photos, strictOrLax, mimicType) {
  as.factor(paste(photos$imageableid, strictOrLax, mimicType, sep = '-'))
}

# Returns "sort-of" species ID as a factor. 
# More accurately, species with further differentiation on strict or lax photos and mimic type.
# Hence, 1 individual with both strict and lax outlines will end up in two "sort-of species".
# Similarly, a nymph and an adult of the same species will end up in 2 different "species" 
# if one is a mimic and the other isn't.
speciesIds <- function(fac, strictOrLax, mimicType) {
  as.factor(paste(fac$scientificName, strictOrLax, mimicType, sep = '-'))
}

# Returns a type ID as a factor. 
# More accurately, mimic type (model, x-mimic, non-mimic) with further differentiation on strict or lax photos.
# Hence, there will be potentially be up to 8 types: strict models, lax models, strict spider mimics etc.
# In practice, there will be no lax models or non-mimics.
typeIds <- function(strictOrLax, mimicType) {
  as.factor(paste(strictOrLax, mimicType, sep = '-'))
}

# Given some photo source field values, returns the ids of the source photos.
# The field is assumed to have the format "Photo <id>"
SourcePhotoId <- function(source) {
  as.numeric(gsub('.*Photo +', '', source))
}


#!Rscript

# Attempts to calculate percentage of time that specimen is stopped in prey videos

LoadFns('motion')
LoadFns('download')

#.host <- '130.56.244.191'
.host <- 'localhost'

.url <- function(path) sprintf('http://%s/%s', .host, path)


#########################################################################


# Load up all prey videos
prey.videos <- read.csv(.url('photos.csv?ptype=Prey&ftype=video'), stringsAsFactors = FALSE)
specimen.ids.with.videos <- unique(prey.videos$imageableid)
specimens <- read.csv(.url(paste0('specimens.csv?id=[', paste(specimen.ids.with.videos, collapse = ","), "]")))

# get videos which have not yet had stopped time calculated
todo <- prey.videos[grep('stopped', prey.videos$description, ignore.case = T, invert = T),]

# Download videos
files <- JDownload(todo$url, verbose = F)

# Now for each video...
for(i in 1:nrow(todo)) {
  # Convert track to CSV
  vidFile <- files[i]
  system2('c:/Jim/uni/apps/VideoAnalysis/run.bat', c('-defaults', 'stopped', vidFile))
  csvFile <- paste0(file_path_sans_ext(vidFile), '.csv')
  track <- readTrack(csvFile, todo, rotate = FALSE)
  # Calculate % time that track is stopped
  sp <- todo[i,]
  cat(sprintf("Specimen %d, video %d, stopped %d%%\n", sp$imageableid, sp$id, stoppedPercent(track)))
}
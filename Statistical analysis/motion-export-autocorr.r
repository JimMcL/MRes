#!Rscript
LoadFns('motion')
LoadFns("sample-dbs")
LoadFns('rediscretization')
LoadFns('straightness')

tracks <- MGetAllTracks()
fn <- gsub('.csv$', '-corr.csv', tracks$trackFile)
x <- sapply(1:nrow(tracks), function(i) {
  cat(sprintf('Correlating track %d of %d (video %d, track length %d)\n', i, nrow(tracks), tracks[i,]$id, length(tracks[i,]$track[[1]]$polar)))
  corr <- MAutocorrelateTrack(tracks[i,])
  write.csv(corr, fn[i], row.names = FALSE)
})

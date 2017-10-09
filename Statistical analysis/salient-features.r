

motion <- read.csv('http://localhost/ac/imperfection-per-species-motion.csv')
morpho <- read.csv('http://localhost/ac/imperfection-per-species-angle-morpho.csv')
morpho <- morpho[morpho$angle == 'Dorsal' & morpho$outline == 'actual',]

speciesInCommon <- intersect(motion$species, morpho$species)
common <- merge(motion, morpho, by = 'species')
common <- common[common$type == MTP_SPIDER,]
l <- lm(imperfection.y ~ imperfection.x, data = common)
plot(imperfection.y ~ imperfection.x, data = common)
abline(l, col = 'red')
summary(l)

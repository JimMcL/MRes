# Code to compare the performance of LDA and PCA for discriminating mimics from models.


library(MASS)

# Comparison using the simple data model from Sherratt, T.N. & Peet-Paré (2017).
# Result: PCA and LDA both perform optimally, just using colour to discriminate between classes.

nMimics <- 500
nModels <- 500

nMimicsWithC <- round(nMimics * 1 / 10)
nModelsWithC <- round(nModels * 9 / 10)
nMimicsWithP <- round(nMimics * 1 / 25)
nModelsWithP <- round(nModels * 3 / 50)

.assignTrait <- function(n, nWith, code) sample(c(rep(paste0(code, '+'), nWith), rep(paste0(code, '-'), n - nWith)))
.build <- function(type, n, withC, withP) data.frame(type = type, colour = .assignTrait(n, withC, 'C'), pattern = .assignTrait(n, withP, 'P'))
.payOff <- function(alg, mimicsAttacked, modelsAttacked) cat(sprintf("Payoff for %s: %d (mimics attacked) - %d (models attacked) = %d\n", alg, mimicsAttacked, modelsAttacked, mimicsAttacked - modelsAttacked))

prey <- rbind(
  .build('mimic', nMimics, nMimicsWithC, nMimicsWithP),
  .build('model', nModels, nModelsWithC, nModelsWithP)
)

### LDA
l <- lda(type ~ ., prey)
lp <- predict(l)
.payOff('LDA', sum(lp$class == 'mimic' & prey$type == 'mimic'), sum(lp$class == 'mimic' & prey$type == 'model'))

### PCA
pc <- prcomp(~ as.numeric(colour) + as.numeric(pattern), data = prey)
.payOff('PCA (PC1)', sum(pc$x[,1] < 0 & prey$type == 'mimic'), sum(pc$x[,1] < 0 & prey$type == 'model'))

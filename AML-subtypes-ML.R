# Ctrl + L to clear console in RStudio
setwd("~/Programming/R&RStudio") # Set working directory

# Load required libraries
library(GEOquery)     
library(limma)        
library(randomForest) 
library(caret)        
library(glmnet)       
library(ggplot2)      
library(pheatmap)     
library(e1071)        
library(smotefamily)  
library(reshape2)     
library(AnnotationDbi)
library(hgu133plus2.db)

# Enable new plots in RStudio
options(device.ask.default = TRUE)

# Load GSE13159 dataset from local directory
gset <- readRDS("GSE13159.rds")

# Extract and clean expression data
exprsData <- exprs(gset)
exprsData <- as.data.frame(exprsData)
exprsData <- na.omit(exprsData)

# Pre-filter data distribution
cat("Summary of rowMeans(exprsData) before filtering:\n")
print(summary(rowMeans(exprsData)))
hist(rowMeans(exprsData), breaks = 50, main = "Distribution of Row Means (Unfiltered)")

# Filter low-expression probes using median threshold
threshold <- median(rowMeans(exprsData))
exprsFiltered <- exprsData[rowMeans(exprsData) > threshold, ]

# Post-filter data distribution
hist(rowMeans(exprsFiltered), breaks = 50, main = "Filtered Distribution of Row Means")

# Load and Filter Phenotype Data
phenoData <- pData(gset)
phenoData <- phenoData[grepl("AML", phenoData$`leukemia class:ch1`), ]

# Match samples between expression data and phenotype data
commonSamples <- intersect(colnames(exprsFiltered), rownames(phenoData))
exprsFiltered <- exprsFiltered[, commonSamples]
phenoData <- phenoData[commonSamples, ]

# Output unique AML subtypes
uniqueClasses <- unique(phenoData$`leukemia class:ch1`)
leukemiaClassCounts <- table(phenoData$`leukemia class:ch1`)
cat("Unique AML subtypes and counts:\n")
print(leukemiaClassCounts)

# Extract AML subtype labels
labels <- as.factor(phenoData$`leukemia class:ch1`)

# LASSO Feature Selection
lassoModel <- cv.glmnet(as.matrix(t(exprsFiltered)), labels, alpha = 1, family = "multinomial")
selectedProbes <- coef(lassoModel, s = "lambda.min")
selectedProbes <- rownames(selectedProbes[[1]])[which(selectedProbes[[1]] != 0)]

exprsFilteredLasso <- exprsFiltered[selectedProbes, ]
exprsFilteredLasso <- exprsFilteredLasso[!rowSums(is.na(exprsFilteredLasso)), ]
cat("Number of selected probes after LASSO:", length(selectedProbes), "\n")

# Probe-to-gene mapping
probeToGene <- mapIds(
  hgu133plus2.db, keys = rownames(exprsFilteredLasso),
  column = "SYMBOL", keytype = "PROBEID", multiVals = "first"
)

rownames(exprsFilteredLasso) <- ifelse(
  is.na(probeToGene), rownames(exprsFilteredLasso),
  paste0(rownames(exprsFilteredLasso), "_", probeToGene)
)

# Train-test split (70-30)
set.seed(100)

trainIndex <- createDataPartition(labels, p = 0.7, list = FALSE)
trainData <- exprsFilteredLasso[, trainIndex]
testData  <- exprsFilteredLasso[, -trainIndex]
trainLabels <- labels[trainIndex]
testLabels  <- labels[-trainIndex]

commonGenes <- intersect(rownames(trainData), rownames(testData))
trainData <- trainData[commonGenes, ]
testData  <- testData[commonGenes, ]

# Train baseline Random Forest (pre-SMOTE)
rfModelPreSmote <- randomForest(x = t(trainData), y = trainLabels, importance = TRUE)
testPredictionsPreSmote <- predict(rfModelPreSmote, newdata = t(testData))
confMatrixPreSmote <- confusionMatrix(testPredictionsPreSmote, testLabels)
print(confMatrixPreSmote)

accuracyPreSmote <- confMatrixPreSmote$overall["Accuracy"]
cat("Pre-SMOTE RF Classification Accuracy:", round(accuracyPreSmote, 3), "\n")

# Apply SMOTE to balance training data
trainDf <- data.frame(t(trainData))
trainDf$Label <- trainLabels

X <- trainDf[, -ncol(trainDf)]
y <- trainDf$Label

smoteResult <- SMOTE(X, y, K = 5, dup_size = 1)
trainDataSmote <- smoteResult$data

# Remove "X" prefixes from column names (if any)
colnames(trainDataSmote) <- sub("^X", "", colnames(trainDataSmote))

colnames(trainDataSmote)[ncol(trainDataSmote)] <- "Label"
trainDataSmote$Label <- as.factor(trainDataSmote$Label)

# Reformat SMOTE-adjusted data
trainDataPostSmote <- t(trainDataSmote[, -ncol(trainDataSmote)])
trainLabelsPostSmote <- trainDataSmote$Label

# Train post-SMOTE RF model
commonGenesSMOTE <- intersect(rownames(trainDataPostSmote), rownames(testData))
trainDataPostSmote <- trainDataPostSmote[commonGenesSMOTE, , drop = FALSE]
testData           <- testData[commonGenesSMOTE, , drop = FALSE]

# Remove more "X" prefixes
rownames(trainDataPostSmote) <- sub("^X", "", rownames(trainDataPostSmote))
rownames(testData)           <- sub("^X", "", rownames(testData))

rfModelPostSmote <- randomForest(x = t(trainDataPostSmote), y = trainLabelsPostSmote, importance = TRUE)

testPredictionsPostSmote <- predict(rfModelPostSmote, newdata = t(testData))
confMatrixPostSmote <- confusionMatrix(testPredictionsPostSmote, testLabels)
print(confMatrixPostSmote)

accuracyPostSmote <- confMatrixPostSmote$overall["Accuracy"]
cat("Post-SMOTE RF Classification Accuracy (Untuned):", round(accuracyPostSmote, 3), "\n")

# Convert data to samples x features for caret
rfTrainX <- t(trainDataPostSmote)
rfTrainY <- trainLabelsPostSmote

# Hyperparameter tuning for RF (10-fold cross-validation)
rfCtrl <- trainControl(method = "cv", number = 10)
rfGrid <- expand.grid(mtry = c(2, 3, 4, 5, 6))

set.seed(100)
rfCaretModel <- train(
  x = rfTrainX, y = rfTrainY,
  method = "rf", metric = "Accuracy", tuneGrid = rfGrid,
  trControl = rfCtrl, ntree = 500
)

# Evaluate final tuned RF on testData
rfPredictionsCaret <- predict(rfCaretModel, t(testData))
confMatrixRfCaret <- confusionMatrix(rfPredictionsCaret, testLabels)
print(confMatrixRfCaret)

accuracyRfCaret <- confMatrixRfCaret$overall["Accuracy"]
cat("Tuned RF Classification Accuracy:", round(accuracyRfCaret, 3), "\n")

# Hyperparameter tuning for SVM (10-fold cross-validation)
svmCtrl <- trainControl(method = "cv", number = 10)
svmGrid <- expand.grid(sigma = 2^(-5:-1), C = 2^(0:4))

set.seed(100)
svmCaretModel <- train(
  x = rfTrainX, y = rfTrainY,
  method = "svmRadial", metric = "Accuracy",
  tuneGrid = svmGrid, trControl = svmCtrl
)

# Evaluate final tuned SVM on testData
svmPredictions <- predict(svmCaretModel, t(testData))
confMatrixSvm <- confusionMatrix(svmPredictions, testLabels)
print(confMatrixSvm)

# Feature importance plot
varImpPlot(rfModelPostSmote, n.var = 20, main = "Top 20 Important Probes (RF)")

# PCA visualization
pcaRf <- prcomp(t(testData))
pcaDataRf <- data.frame(
  Sample = rownames(pcaRf$x),
  PC1 = pcaRf$x[, 1], PC2 = pcaRf$x[, 2],
  Label = testLabels
)

ggplot(pcaDataRf, aes(x = PC1, y = PC2, color = Label)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot - Test Data (Color by True Label)", x = "PC1", y = "PC2")

# Function to normalize confusion matrices row-wise
normalizeConfMatrix <- function(confMatrix) {
  mat <- as.matrix(confMatrix$table)  # Convert to matrix
  rowTotals <- rowSums(mat)           # Sum across rows
  rowTotals[rowTotals == 0] <- 1      # Prevent division by zero
  sweep(mat, 1, rowTotals, FUN = "/") # Row-normalized matrix
}

# Function to plot heatmaps
plotHeatmap <- function(matrix, title, colorPalette) {
  pheatmap(
    matrix, 
    cluster_rows = FALSE, 
    cluster_cols = FALSE,
    color = colorRampPalette(colorPalette)(100),
    display_numbers = TRUE,
    main = title
  )
}

# Normalize confusion matrices
rfPreNorm <- normalizeConfMatrix(confMatrixPreSmote)
rfPostNorm <- normalizeConfMatrix(confMatrixPostSmote)
rfCaretNorm <- normalizeConfMatrix(confMatrixRfCaret)
svmNorm <- normalizeConfMatrix(confMatrixSvm)

# Plot all heatmaps
plotHeatmap(rfPreNorm, "Confusion Matrix - Pre-SMOTE RF (Row-Normalized)", c("white", "red"))
plotHeatmap(rfPostNorm, "Confusion Matrix - Post-SMOTE RF (Row-Normalized)", c("white", "red"))
plotHeatmap(rfCaretNorm, "Confusion Matrix - Tuned RF (Row-Normalized)", c("white", "red"))
plotHeatmap(svmNorm, "Confusion Matrix - Tuned SVM (Row-Normalized)", c("white", "blue"))

# Classification accuracy bar char
allAccuracies <- data.frame(
  Model = c("Pre-SMOTE RF", "Post-SMOTE (Untuned) RF", "Tuned RF", "Tuned SVM"),
  Accuracy = c(
    accuracyPreSmote,
    accuracyPostSmote,
    accuracyRfCaret,
    confMatrixSvm$overall["Accuracy"]
  )
)

ggplot(allAccuracies, aes(x = Model, y = Accuracy)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  ylim(c(0, 1)) +
  labs(title = "Classification Accuracies", x = "Models", y = "Accuracy")

# Print completion message
cat("\nDONE: AML classification pipeline executed successfully.\n")
cat("Pre-SMOTE RF:", round(confMatrixPreSmote$overall["Accuracy"], 4),
    "| Post-SMOTE RF:", round(confMatrixPostSmote$overall["Accuracy"], 4),
    "| Tuned RF:", round(confMatrixRfCaret$overall["Accuracy"], 4),
    "| Tuned SVM:", round(confMatrixSvm$overall["Accuracy"], 4), "\n")

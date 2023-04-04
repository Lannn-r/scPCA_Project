library(zellkonverter)
library(SingleCellExperiment)
library(scuttle)
#Load in raw data
sce = readH5AD("data/tm_brain/mice_brain.h5ad")
sce
# Rewrite column names to match the common schema, etc
sce$ori_cell_id<-sce$cell
sce$cell<-NULL
sce$bio_replicates<-sce$mouse.id
sce$mouse.id<-NULL 
colnames(sce)=sce$ori_cell_id

# Unlog original normalization and us log-transformed Normalization
unlog = exp(as.matrix(assay(sce, 'X')))-1 #coerce to be dense matrix
counts = t(t(unlog/1000)*colData(sce)$n_counts)
counts(sce) = Matrix::Matrix(counts, sparse = TRUE)
sce=logNormCounts(sce)

# Save processed data in the data directory, make sure it gets synced onto Box
saveRDS(sce, "data/processed/tm_brain.rds")

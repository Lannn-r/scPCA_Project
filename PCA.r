library(SingleCellExperiment)
library(bluster)
library(BiocNeighbors)


trainPCA <- function(sce.train, exprs_values = "logcounts", name = 'PCA', ...){
  sce.runpca<-runPCA(sce.train, exprs_values = exprs_values, name = name, ...) 
  sce.rotation<-attr(reducedDim(sce.runpca, name), "rotation")
  sce.HVG<-rownames(sce.rotation)
  mu_hat<-rowMeans(assay(sce.train[sce.HVG, ], exprs_values))
  metadata(sce.runpca) = c(metadata(sce.runpca), list(pca_rotation = sce.rotation, 
                                                      pca_hvg = sce.HVG, mu_hat = mu_hat, exprs_values = exprs_values))
  return(sce.runpca)
}
#First apply runPCA() to return the same kind SingleCellExperiment object sce.runpca that contains the numeric
#matrix of every cell's coordinates in each of the PCs in sce.train in reducedDim(), then extract the rotation matrix
#as sce.rotation. Extracting row names from sce.rotation as highly variable genes, and generate row means for 
#sce.train(which is sce.HVG). Finally reset the names of attributes back to sce.runpca and set sce.runpca as a Raster
#object and return it. This function works as training the original input metadata as a new object encoding the 
#results from a PCA run.

is_trained_pca = function(sce){
  return(!is.null(metadata(sce)$mu_hat))
}
#This function checks all the valid values (not missing) in row means from the input sce. can also check if a SingleCellExperiment
#function starts by 'is_...' is trying to check if it's true or false

projectPCA<-function(sce.test, trainedPCA, name){
  stopifnot(is_trained_pca(trainedPCA))
  x=assay(sce.test, metadata(trainedPCA)$exprs_values)
  x_subset=x[metadata(trainedPCA)$pca_hvg, ]
  x_center=t(x_subset-metadata(trainedPCA)$mu_hat)
  reducedDim(sce.test, name) = x_center %*% as.matrix(metadata(trainedPCA)$pca_rotation)
  # consider adding some sort of hash/heuristic to indicate successful retraining.. 
  return(sce.test)
}
#First combine the sce.test (a SingleCellExperiment object) and the string indicating the part of trainedPCA contains 
#the expression values as x. Then extract the highly variable genes from trainedPCA as x_subset, center x by transposing
#the matrix of the differences between HVG and row means. Finally multiplicate x_center and the rotation matrix from
#trainedPCA together and save it into reducedDim() with a dim name and return. The function works as performing PCA
#on the test data that was projected by the trained model.

cluster_sce = function(sce, dimred = 'PCA', label_name = 'label', BLUSPARAM = NNGraphParam(cluster.fun = 'louvain'), ...) {
  cluster_ids <- bluster::clusterRows(reducedDim(sce, dimred), BLUSPARAM = BLUSPARAM, ...)
  colLabels(sce) = factor(cluster_ids)
  colData(sce)[[label_name]] =  colLabels(sce)
  sce
}
# what's the rows and columns for reducedDim(sce, dimred): rows are cells, columns are PCs

#Cluster the rows of sce using algorithm BLUSPARAM, factor the clustered rows as ids, and use these ids as
#column labels for sce. The function works as clustering rows of the input data and creating labels for them.

project_label<-function(sce_test, sce_train, dimred_train, dimred_test, label_name = 'label'){
  a<-reducedDim(sce_train, dimred_train) 
  b<-reducedDim(sce_test, dimred_test)
  # right now set k = 1, but eventually we'll want to see what happens
  # with k > 1 and report the majority vote and the concordance
  bNN = BiocNeighbors::queryKNN(a, b, k=1)
  if (is.null(colLabels(sce_train))) stop("Must Provide 'colLabels(sce_train)'") 
  blabel = colLabels(sce_train)[bNN$index]
  #overwrite labels:
  colLabels(sce_test) = blabel
  colData(sce_test)$knn_dist = bNN$distance
  colData(sce_test)[[label_name]] = blabel
  #return result:
  sce_test
}
#choose the right part of the reducedDim of testing/training data as two matrices with equal number and order of 
#columns, then identify 1 nearest neighbors in training data for each point in testing data. Check the column
#labels for training data. Then project the labels of training data to testing data corresponding to each of the
#column that points are. Finally paste the distance to the nearest neighbors to testing data and return the 
#modified testing data. The function works as querying the nearest neighbors for testing data from training data,
#and change the original labels of testing data to be the ones corresponding to neighbor points.

project_UMAP = function(sce_train, sce_test, dimred){
  # we need to make a wrapper around runUMAP so that ret_model = TRUE
  umap_test = umap_transform(reducedDim(sce_test, dimred), model = sce_train$train_umap)
  # need set the reduced dim with umap_test
}
#The function works as transforming the matrix of coordinates in testing data to the training data

align_PCA = function(query, reference){
  modi_query = query
  for(i in 1:ncol(modi_query)){
    if(cor(query[, i], reference[, i])<0){
      modi_query[ ,i] = modi_query[ ,i]*-1 
    }
  }
  print(modi_query)
}
#loop through each of the columns of query (all data) and reference (training data) data, check if the sign of 
#the diagonal of the correlation between query and reference is negative, then flip the sign of that column in 
#the query data. The function works as modifying query data in order to maximize the concordance between query 
#and reference data. 


split_and_label = function(sce, train_idx, dimred_test = 'PCA', trainPCA_args = list(exprs_values = "logcounts"), 
                            cluster_sce_args = list(label_name = 'label', BLUSPARAM = NNGraphParam())) {
  train <- sce[, train_idx]
  test <- sce[, -train_idx]
  ## This function defaults to setting reducedDim to "PCA", we could over-ride this...
  dimred_train = "PCA"
  train <- call_intercalate(trainPCA, sce.train = train, extra = trainPCA_args)
  train <- call_intercalate(cluster_sce, sce = train, dimred = dimred_train, extra = cluster_sce_args)
  
  if (length(train_idx) < ncol(sce)) {
    test <-
      projectPCA(test, train, name = dimred_test)
    clustered_test <-
      project_label(test, train, dimred_train, dimred_test)
  } else{
    clustered_test = NULL
  }
  list(test = clustered_test, train = train)
  #project_test<-projectPCA(zeisel.test_ori, clustered_test, name=name2)
}

vfold_resample = function(sce, fold_map, train_args = NULL, return_train = c('null', 'tbl', 'SingleCellExperiment'), return_test = c('tbl', 'SingleCellExperiment')){
  if(!inherits(fold_map, 'data.frame') || nrow(fold_map) != ncol(sce) || 
     !all(c('cell_idx', 'fold') %in% names(fold_map)) ) 
    stop('`fold_map` must be data.frame with fields `cell_idx` and `fold`.')
  return_train = match.arg(return_train)
  return_test = match.arg(return_test)
  fold_split = split(fold_map, fold_map[['fold']])
  fit = list()
  idx = seq_len(ncol(sce))
  pb = progress::progress_bar$new(total = length(fold_split), format = " processing fold :current [:bar]")
  for(i in seq_along(fold_split)){
    pb$tick()
    train_idx = setdiff(idx, fold_split[[i]][['cell_idx']])
    this_fit = call_intercalate(split_and_label, sce = sce, train_idx = train_idx, extra = train_args)
    fit[[i]] = list(test = switch(return_test, 
                                  tbl = tibble::as_tibble(as.data.frame(colData(this_fit$test), optional = TRUE)), # optional = TRUE: don't mangle column names
                                  SingleCellExperiment = this_fit$test),
                    train = switch(return_train,
                                   null = NULL,
                                   tbl = tibble::as_tibble(as.data.frame(colData(this_fit$train), optional = TRUE)),
                                   SingleCellExperiment = this_fit$train))
  }
  return(fit)
}

concordance_stats = function(ref, query, label = 'label', cell_id = 'cell_id'){
  stopifnot(all.equal(ref[[cell_id]], query[[cell_id]]))
  train_label =  ref[[label]]
  test_label = query[[label]]
    
  train_label = factor(train_label, levels = union(train_label, test_label))
  test_label = factor(test_label, levels = levels(train_label))
  cross_tab = table(train=train_label, test=test_label, exclude=NULL)
  remap = mclust::mapClass(as.character(train_label), as.character(test_label))
  # This maps table _indices_
  error = mclust::classError(as.character(train_label), as.character(remap$bTOa[as.character(test_label)]))
  accuracy = 1-error$errorRate
  permuted_tab = cross_tab
  permuted_tab[,names(remap$bTOa)] = permuted_tab[,as.character(remap$bTOa)]
  train_nclust = length(unique(train_label))
  test_nclust = length(unique(test_label))
  concordance_result = list(permuted_tab = permuted_tab, accuracy_rate=accuracy, 
                            adj_RandIndex=mclust::adjustedRandIndex(train_label, test_label), 
                            train_nclust = train_nclust,
                            test_nclust = test_nclust)
  return(concordance_result)
}


library(dplyr)
library(plyr)
library(stats)
library(utils)
library(caret)
library(reshape2)

# reading in the data
optdigits.train <- read.table('optdigits.train', sep =',')
optdigits.train <- tbl_df(optdigits.train) 
train_data <- select(optdigits.train, -V65) ## 3823 X 64
optdigits.test <- read.table('optdigits.test', sep = ',')
optdigits.test <- tbl_df(optdigits.test) 
test_data <- select(optdigits.test, -V65)


## The function f1 takes a training data and the number of clusters desired
## It implements the K-means clustering, 5 different times and returns a list
## of SSE, cluster centers and final clustering obtained at the end of each round
## of K-means clustering

f1 <- function(tr_data, K){
  set.seed(1080000)
  seed_vector <- sample(1:100000, 5, replace=F)
  seed_vector
  initial_centroid_list <- list()
  for(i in 1:length(seed_vector)){
    set.seed(seed_vector[i])
    randomRows <- c(sample(1:nrow(tr_data), K, replace = FALSE))
    ## print (randomRows)
    initial_centroid_list[[i]] <- tr_data %>% slice(randomRows)
  }
  SSE_VECTOR <- numeric(0)
  final_centroids <- list()
  final_cluster <- list()
  
  for(i in 1:length(initial_centroid_list)){
    initial_centroids <- initial_centroid_list[[i]]
    euc_dist <- as.data.frame(apply(initial_centroids, 1, function(x) 
      (apply(tr_data, 1, function(y) sqrt(sum((x-y)^2))))))
    
    clusters <- data.frame(Clusters = apply(euc_dist, 1, which.min))
    
    df <- cbind(tr_data, clusters)
    
    new_centroids <- list()
    for (j in 1:K){
      a <- filter(df, Clusters == j)
      new_centroids[[j]] <- apply(a[,-65], 2, mean)
    }
    new_centroids_df <- as.data.frame(do.call("rbind", new_centroids))
    
    
    while (identical(initial_centroids, new_centroids_df) == FALSE){
      df <- data.frame()
      initial_centroids <- data.frame()
      initial_centroids <- new_centroids_df
      euc_dist <- as.data.frame(apply(initial_centroids, 1, function(x) 
        (apply(tr_data, 1, function(y) sqrt(sum((x-y)^2))))))
      
      
      clusters <- data.frame(Clusters = apply(euc_dist, 1, which.min))
      
      df <- cbind(tr_data, clusters)
      new_cen <- list()
      new_centroids_df <- data.frame()
      for (j in 1:K){
        b <- filter(df, Clusters == j)
        new_cen[[j]] <- apply(b[,-65], 2, mean)
      }
      new_centroids_df <- as.data.frame(do.call("rbind", new_cen))  
    }
    ## SSE
    SSE <- numeric(0)
    for (k in 1:K){
      c <- filter(df, Clusters == k)
      dist <- apply(c[,-65], 1, function(x, y) sum((x -y)^2), y=new_centroids_df[k, 1:64])
      SSE <- append(SSE, dist)
    }
    SSE_final <- sum(SSE)  
    SSE_VECTOR <- append(SSE_VECTOR, SSE_final)
    ## final centroids
    final_centroids[[i]] <- new_centroids_df
    final_cluster[[i]] <- clusters
  }
  list(SSE_VECTOR = SSE_VECTOR, final_centroids=final_centroids, final_cluster=final_cluster)
}

## Result from function f1 
result1_10 <- f1(tr_data = train_data, K = 10)

vec_10 <- result1_10$SSE_VECTOR   ## vector of SSE
cen_list_10 <- result1_10$final_centroids  ## list of cluster centers
clust_list_10 <- result1_10$final_cluster   ## list of final clusters


## The function f2 finds the run with the lowest SSE
## It gives out the smallest SSE and the SSS for the clustering with the 
## lowest SSE
## It also gives the final cluster centers(lowest SSE) and the 
## final clustering(lowest SSE)

f2 <- function(vec, cen_list, clust_list) {
  SSE <- min(vec)
  fc <- clust_list[[which.min(vec)]]
  df <- cen_list[[which.min(vec)]]
  ind <- t(combn(nrow(df),2))
  out <- apply(ind, 1, function(x) sum((df[x[1],] - df[x[2],])^2))
  SSS <- sum(out)
  list(SSE=SSE, SSS=SSS, final_cluster_centers=df, final_clusters=fc)
}


result2_10 <- f2(vec=vec_10, cen_list=cen_list_10, clust_list=clust_list_10)

SSE_10 <- result2_10$SSE  ## SSE
SSE_10
SSS_10 <- result2_10$SSS  ## SSS
SSS_10
final_cluster_centers_10 <- result2_10$final_cluster_centers ##  final cluster centers
final_clusters_10 <- result2_10$final_clusters

## This function gives the mean entropy of clustering
f3 <- function(cluster, df, K){
  belong <- cbind(df[,65], cluster)
  num <- list()
  for (s in 1:K){
    z <- filter(belong, Clusters == s)
    num_vec <- numeric(0)
    for (t in 0:9){
      nv <- length(which(z$V65 == t))
      num_vec <- append(num_vec, nv)
    }
    num[[s]] <- num_vec
  }
  df3 <- as.data.frame(do.call("rbind", num))  
  total <- apply(df3, 1, sum)
  df3 <- cbind(df3, total)
  colnames(df3) <- c(0:9, "Total")  
  df4 <- t(apply(df3, 1, function(x) x/x[length(x)]))
  df5 <- t(apply(df4, 1, function(x) -x*log2(x)))
  entropy <- apply(df5, 1, sum, na.rm = TRUE)
  entropy
  fraction <- (df3$Total)/sum(df3$Total)
  df6 <- cbind(fraction, entropy)
  mean_entropy <- sum(apply(df6, 1, function(x) (x[1]*x[2])))
  list (mean_entropy = mean_entropy, cluster_df=df3)  
}

result3_10 <- f3(cluster=final_clusters_10, df=optdigits.train, K=10 )
mean_entropy_clustering_10 <- result3_10$mean_entropy ## mean entropy of clustering
mean_entropy_clustering_10  


## Associate each cluster center with the most frequent class

df_10 <- result3_10$cluster_df

max_class_10 <- numeric(0)
for (i in 1:10){
  name <- colnames(df_10[which.max(df_10[i, 1:10])])
  max_class_10[i] <- name
}
max_class_10


cluster_name_10 <- c('Clus1', 'Clus2', 'Clus3', 'Clus4', 'Clus5', 'Clus6', 
                  'Clus7', 'Clus8', 'Clus9', 'Clus10')

cluster_class_10 <- cbind(cluster_name_10, as.numeric(max_class_10))
cluster_class_10

translator_function_10 = function(element) {
  switch(element,
         '1' = 5,
         '2' = 9,
         '3' = 6,
         '4'= 0,
         '5' = 7,
         '6' = 1,
         '7' = 3,
         '8' = 2,
         '9' = 8,
         '10' = 4)
}

true_pred_10 <- cbind(optdigits.train[,65], final_clusters_10)

switch_apply_10 = unlist(sapply(true_pred_10$Clusters, translator_function_10))

true_pred_10 <- cbind(true_pred_10, switch_apply_10)

confMat_train_10 <- confusionMatrix(true_pred_10$V65, true_pred_10$switch_apply)

## Running the K-means algorithm on the test data set
euc_dist_test_10 <- as.data.frame(apply(final_cluster_centers_10, 1, function(x) 
  (apply(test_data, 1, function(y) sqrt(sum((x-y)^2))))))

clusters_test_10 <- data.frame(Clusters = apply(euc_dist_test_10, 1, which.min))

true_pred_test_10 <- cbind(optdigits.test[,65], clusters_test_10)

switch_apply_test_10 = unlist(sapply(true_pred_test_10$Clusters, translator_function_10))

true_pred_test_10 <- cbind(true_pred_test_10, switch_apply_test_10)

confMat_test_10 <- confusionMatrix(true_pred_test_10$V65, true_pred_test_10$switch_apply_test)


## Experiment-2--------------------------------------------------------------------##
## K-means clustering with 30 clusters

result1_30 <- f1(tr_data = train_data, K = 30)

vec_30 <- result1_30$SSE_VECTOR   ## vector of SSE
cen_list_30 <- result1_30$final_centroids  ## list of cluster centers
clust_list_30 <- result1_30$final_cluster   ## list of final clusters


result2_30 <- f2(vec=vec_30, cen_list=cen_list_30, clust_list=clust_list_30)

SSE_30 <- result2_30$SSE  ## SSE
SSE_30
SSS_30 <- result2_30$SSS  ## SSS
SSS_30
final_cluster_centers_30 <- result2_30$final_cluster_centers ##  final cluster centers
final_clusters_30 <- result2_30$final_clusters

result3_30 <- f3(cluster=final_clusters_30, df=optdigits.train, K=30 )
mean_entropy_clustering_30 <- result3_30$mean_entropy ## mean entropy of clustering
mean_entropy_clustering_30 

## Associate each cluster center with the most frequent class

df_30 <- result3_30$cluster_df

max_class_30 <- numeric(0)
for (i in 1:30){
  name <- colnames(df_30[which.max(df_30[i, 1:10])])
  max_class_30[i] <- name
}
max_class_30

cluster_name_30 <- c('Clus1', 'Clus2', 'Clus3', 'Clus4', 'Clus5', 'Clus6', 
                     'Clus7', 'Clus8', 'Clus9', 'Clus10', 'Clus11', 'Clus12', 'Clus13', 
                     'Clus14', 'Clus15', 'Clus16', 'Clus17', 'Clus18', 'Clus19', 'Clus20',
                     'Clus21', 'Clus22', 'Clus23', 'Clus24', 'Clus25',
                     'Clus26', 'Clus27', 'Clus28', 'Clus29', 'Clus30')

cluster_class_30 <- cbind(cluster_name_30, as.numeric(max_class_30))
cluster_class_30

translator_function_30 = function(element) {
  switch(element,
         '1' = 1,
         '2' = 7,
         '3' = 2,
         '4'= 6,
         '5' = 1,
         '6' = 7,
         '7' = 6,
         '8' = 9,
         '9' = 5,
         '10' = 7,
         '11' = 9,
         '12'= 1,
         '13' = 1,
         '14'= 3,
         '15'= 1,
         '16' = 2,
         '17' = 8,
         '18' = 0,
         '19' = 9,
         '20'= 2,
         '21'= 4,
         '22' = 4,
         '23' = 0,
         '24'= 4,
         '25' = 4,
         '26' = 3,
         '27' = 5,
         '28' = 8,
         '29' = 0,
         '30' = 3)
}

true_pred_30 <- cbind(optdigits.train[,65], final_clusters_30)

switch_apply_30 = unlist(sapply(true_pred_30$Clusters, translator_function_30))

true_pred_30 <- cbind(true_pred_30, switch_apply_30)

confMat_train_30 <- confusionMatrix(true_pred_30$V65, true_pred_30$switch_apply)

## Running the K-means (K = 30) algorithm on the test data set

euc_dist_test_30 <- as.data.frame(apply(final_cluster_centers_30, 1, function(x) 
  (apply(test_data, 1, function(y) sqrt(sum((x-y)^2))))))

clusters_test_30 <- data.frame(Clusters = apply(euc_dist_test_30, 1, which.min))

true_pred_test_30 <- cbind(optdigits.test[,65], clusters_test_30)

switch_apply_test_30 = unlist(sapply(true_pred_test_30$Clusters, translator_function_30))

true_pred_test_30 <- cbind(true_pred_test_30, switch_apply_test_30)

confMat_test_30 <- confusionMatrix(true_pred_test_30$V65, true_pred_test_30$switch_apply_test)


a <- as.numeric(final_cluster_centers_30[18, 1:64])
b <- matrix(a, nrow=8, ncol=8)
png("simpleIm_30_0.png")
par(mar = rep(0, 4))
image(b, axes = FALSE, col = grey(seq(0, 1, length = 256)))
dev.off()

c <- as.numeric(final_cluster_centers_10[9, 1:64])
d <- matrix(a, nrow=8, ncol=8)
png("simpleIm_10_8.png")
par(mar = rep(0, 4))
image(b, axes = FALSE, col = grey(seq(0, 1, length = 256)))
dev.off()

 



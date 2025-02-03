# Helper Functions

norm_vec <- function(x) sqrt(sum(x^2))

multivariate_Z_test <- function(X, linkage, K, k1, k2, sig) { 
  q <- ncol(X)
  hcl <- fastcluster::hclust(dist(X)^2, method=linkage)
  hcl_at_K <- cutree(hcl, K)
  
  diff_means <- colMeans(X[hcl_at_K == k1, , drop=F]) - colMeans(X[hcl_at_K == k2, , drop=F])
  stat <- norm_vec(diff_means) 
  n1 <- sum(hcl_at_K == k1)  
  n2 <- sum(hcl_at_K == k2) 
  squared_norm_nu <- 1/n1 + 1/n2
  scale_factor <- squared_norm_nu*sig^2
  
  pval <- 1 - pchisq(stat^2/scale_factor, df=q)
  return(list(stat=stat, pval=pval))
}

# Simulator Methods

## Wald test
make_multivariate_Z_test <- function(link, K) {
  new_method(name = sprintf("%s-Z-test-K-%s", link, K),
             label = sprintf("Multivariate Z test (%s, K=%s)", link, K),
             settings = list(link = link, K=K),
             method = function(model, draw, link, K) {
               X <- draw$data
               
               pair <- sample(1:K, 2)
               k1 <- pair[1] 
               k2 <- pair[2]
               
               results <- multivariate_Z_test(X, link, K, k1, k2, model$sig)
               
               hcl <- fastcluster::hclust(dist(X)^2, method=link)
               hcl_at_K <- cutree(hcl, K)
               n1 <- sum(hcl_at_K == k1)  
               n2 <- sum(hcl_at_K == k2) 
               
               true_means <- draw$means[draw$clusters, ]
               effect <- norm_vec(colMeans(true_means[hcl_at_K == k1, , drop=F]) - colMeans(true_means[hcl_at_K == k2, , drop=F]))/model$sig
               crosstab <- table(hcl_at_K, draw$clusters)
               recover <- ifelse(sum(crosstab != 0) == 3, 1, 0)
               npair <- c(n1, n2)
               
               return(list(results=results, effect=effect, recover=recover, npair=npair))
             })
}

## Their test with known variance
make_iso <- function(link, K) {
  new_method(name = sprintf("%s-iso-test-K-%s", link, K),
             label = sprintf("Isotropic (%s, K=%s)", link, K),
             settings = list(link = link, K=K),
             method = function(model, draw, link, K) {
               X <- draw$data
               
               pair <- sample(1:K, 2)
               k1 <- pair[1]
               k2 <- pair[2]
               
               library(fastcluster)
               hcl <- hclust(dist(X)^2, method=link)
               
               if(link == "complete") {
                 results <- test_complete_hier_clusters_approx(X, hcl, K, k1, k2, sig=model$sig)
               } else {
                 results <- test_hier_clusters_exact(X, link, hcl, K, k1, k2, sig=model$sig)
               }
               
               return(list(results=results))
             })
}

## Methods for hc linkage analysis

multivariate_Z_test_methods <- sapply(c("single", "average", "centroid", "complete"), 
                                      function(link) make_multivariate_Z_test(link, 3))

selective_methods <- sapply(c("single", "average", "centroid", "complete"),
                            function(link) make_iso(link, 3))
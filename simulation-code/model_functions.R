# ----------------
# Helper functions
# ----------------

## H0 1 cluster
generate_null_data <- function(n, q, sig) { 
  dat <- matrix(rnorm(n*q, sd=sig), n, q) 
  
  return(list(data=dat, means=matrix(0, 1, q), clusters=rep(1, n)))
}

## H1 2 clusters
generate_two_clusters_data <- function(n, q, len, sig, id=NULL) { 
  cl <- rep_len(1:2, length.out=n)
  
  a <- len/2
  mu <- rbind(c(-a, rep(0, q-1)), c(a, rep(0, q-1)))
  
  dat <- matrix(rnorm(n*q, sd=sig), n, q) + mu[cl, ]
  
  return(list(data=dat, means=mu, clusters=cl))
}

## H1 3 clusters
generate_three_clusters_data <- function(n, q, len, sig, id=NULL) { 
  cl <- rep_len(1:3, length.out=n)
  
  a <- len/2
  mu <- rbind(c(-a, rep(0, q-1)), c(rep(0, q-1), sqrt(3)*a), c(a, rep(0, q-1)))
  
  dat <- matrix(rnorm(n*q, sd=sig), n, q) + mu[cl, ]
  
  return(list(data=dat, means=mu, clusters=cl))
}

## H1 3 clusters unevenly spaced data
generate_three_unevenly_spaced_clusters_data <- 
  function(n, q, len, sig, id = NULL) { 
  cl <- rep_len(1:3, length.out = n)
  
  a <- len / 2
  # Uneven cluster centers.
  #   Cluster 1: (-a, 0)
  #   Cluster 2: (0, sqrt(3)*a)
  #   Cluster 3: (1.5*a, 0)
  mu <- rbind(
    c(-a, rep(0, q - 1)),
    c(0, sqrt(3) * a, rep(0, q - 2)),
    c(1.5 * a, rep(0, q - 1))
  )
  
  # Generate data points from a normal distribution (with standard deviation sig)
  # and shift them by the cluster center based on cl.
  dat <- matrix(rnorm(n * q, sd = sig), n, q) + mu[cl, ]
  
  return(list(data = dat, means = mu, clusters = cl))
  }

## H1 3 clusters with different variance across clusters
generate_three_clusters_diff_var_data <- function(n, q, len, sig, id = NULL) { 

  cl <- rep_len(1:3, length.out = n)
  
  a <- len / 2
  mu <- rbind(
    c(-a, rep(0, q - 1)),
    c(rep(0, q - 1), sqrt(3) * a),
    c(a, rep(0, q - 1))
  )
  
  dat <- matrix(0, n, q)
  
  for(i in 1:3) {
    idx <- which(cl == i)
    # Generate noise with the standard deviation for cluster i.
    noise <- matrix(rnorm(length(idx) * q, sd = sig[i]), nrow = length(idx), ncol = q)
    dat[idx, ] <- noise + matrix(rep(mu[i, ], each = length(idx)), ncol = q, byrow = FALSE)
  }
  
  return(list(data = dat, means = mu, clusters = cl))
}

# -------------------------
# Simulator model functions
# -------------------------

## 1 clusters
make_null_mod <- function(n, q, sig) { 
  new_model(name="null-mod", 
            label=sprintf("Global null data (n=%s, q=%s, sig=%s)", n, q,  sig), 
            params=list(n=n, q=q,  sig=sig), 
            simulate=function(nsim) { 
              lapply(1:nsim, function(x) generate_null_data(n, q, sig))
            })  
}

## 2 clusters
make_two_clusters_mod_with_id <- function(n, q, len, sig, id) { 
  new_model(name="two-clusters-mod-id", 
            label=sprintf("Two clusters data (n=%s, q=%s, len=%s, sig=%s, id=%s)", n, q, len, sig, id), 
            params=list(n=n, sig=sig, id=id), 
            simulate=function(nsim) { 
              lapply(1:nsim, function(x) generate_two_clusters_data(n, q, len, sig, id))
            })  
}

## 3 clusters, baseline
make_three_clusters_mod_with_id <- function(n, q, len, sig, id) { 
  new_model(name="three-clusters-mod-id", 
            label=sprintf("Three clusters data (n=%s, q=%s, len=%s, sig=%s, id=%s)", n, q, len, sig, id), 
            params=list(n=n, sig=sig, id=id), 
            simulate=function(nsim) { 
              lapply(1:nsim, function(x) generate_three_clusters_data(n, q, len, sig, id))
            })  
}

## 3 clusters, unevenly spaced means
make_three_unevenly_spaced_clusters_mod_with_id <- function(n, q, len, sig, id) { 
  new_model(name="three-clusters-mod-id", 
            label=sprintf("Three clusters data (n=%s, q=%s, len=%s, sig=%s, id=%s)", n, q, len, sig, id), 
            params=list(n=n, sig=sig, id=id), 
            simulate=function(nsim) { 
              lapply(1:nsim, function(x) generate_three_unevenly_spaced_clusters_data(n, q, len, sig, id))
            })  
}

## 3 clusters, different variance
make_three_clusters_diff_var_mod_with_id <- function(n, q, len, sig, id) { 
  new_model(name="three-clusters-mod-id", 
            label=sprintf("Three clusters data (n=%s, q=%s, len=%s, sig=%s, id=%s)", n, q, len, sig, id), 
            params=list(n=n, sig=sig, id=id), 
            simulate=function(nsim) { 
              lapply(1:nsim, function(x) generate_three_clusters_diff_var_data(n, q, len, sig, id))
            })  
}



## @knitr init

# QQ-plots of the Wald test p-values and their p-values against the uniform.
# Under H0.

n <- 100
nfeat <- 2
sig <- 1

source("~/Documents/RDirectory/clusterpval-assumptions/simulation-code/model_functions.R")
source("~/Documents/RDirectory/clusterpval-assumptions/simulation-code/type1_method_functions.R")
source("~/Documents/RDirectory/clusterpval-assumptions/simulation-code/eval_functions.R")

name_of_simulation <- paste("naive-type1-n", n, "-q", nfeat, "-sig", sig,  sep="")
sim <- new_simulation(name=name_of_simulation, label=name_of_simulation)

sim <- sim %>% 
  generate_model(make_null_mod, n=n, q=nfeat, sig=sig)
sim <- sim %>% simulate_from_model(nsim=1000, index=1:10) 
sim <- sim %>% run_method(c(multivariate_Z_test_methods,
                            selective_methods), 
                          parallel=list(socket_names=10, libraries=c("clusterpval"))) 
sim <- sim %>% evaluate(c(stat, pval, rejects))
ev <- sim %>% evals %>% as.data.frame

save(ev, file=paste("~/Documents/RDirectory/clusterpval-assumptions/simulation-results", name_of_simulation, ".Rdata", sep=""))

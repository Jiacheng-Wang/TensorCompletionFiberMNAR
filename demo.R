#  #######################################################################
#       File-Name:      demo.R
#       Date:           Wed Apr 27 16:35:36 2022
#       Author:         JCW
#       Purpose:        
#       Input Files:    NONE
#       Output Files:   NONE
#       Data Output:    NONE
#       Dependencies:   NONE
#       Status:         In Progress
#  #######################################################################
set.seed(123)
d1 = d2= d3 = d4 = 10
r1 = r2 = r3 = r4 = 2
# generate the underlying true propensity score
core_tnsr <- as.tensor(array(runif(prod(r1,r2,r3)), dim = c(r1,r2,r3)))
factor1 <- matrix(runif(d1*r1), nrow = d1)
res <- qr(factor1)
factor1 <- qr.Q(res)
factor2 <- matrix(runif(d2*r2), nrow = d2)
res <- qr(factor2)
factor2 <- qr.Q(res)
factor3 <- matrix(runif(d3*r3), nrow = d3)
res <- qr(factor3)
factor3 <- qr.Q(res)

lizt <- list(factor1, factor2, factor3)
tnsr_s_prop <- ttl(core_tnsr,lizt,ms=c(1,2,3))

tnsr_propensity <- as.tensor(pnorm(tnsr_s_prop@data))
tnsr_mask <- array(runif(d1*d2*d3), dim = c(d1,d2,d3))
tnsr_mask <- as.tensor(ifelse(tnsr_propensity@data - tnsr_mask > 0, 1,0))

# generate the tensor observation with the 4th mode having missing fibers
core_tnsr <- as.tensor(array(runif(prod(r1,r2,r3,r4)), dim = c(r1,r2,r3,r4)))
factor1 <- matrix(runif(d1*r1), nrow = d1)
res <- qr(factor1)
factor1 <- qr.Q(res)
factor2 <- matrix(runif(d2*r2), nrow = d2)
res <- qr(factor2)
factor2 <- qr.Q(res)
factor3 <- matrix(runif(d3*r3), nrow = d3)
res <- qr(factor3)
factor3 <- qr.Q(res)
factor4 <- matrix(runif(d4*r4), nrow = d4)
res <- qr(factor4)
factor4 <- qr.Q(res)

lizt <- list(factor1, factor2, factor3, factor4)
tnsr_s <- ttl(core_tnsr,lizt,ms=c(1,2,3,4))
tnsr_noise <- tnsr_s + as.tensor(array(rnorm(d1*d2*d3*d4, mean = 0, sd = 1), dim = c(d1,d2,d3,d4)))
tnsr_obs <- tnsr_noise

miss_index <- which(tnsr_mask@data == 0, arr.ind = T)
d4_index <- matrix(1:d4, ncol = 1)
index_list <- list(miss_index,d4_index)
m.combs <- expand.grid(lapply(index_list, function(x) seq_len(nrow(x))))
total_miss_index <- do.call(cbind, Map(function(a, b) a[b, ], index_list, m.combs))
tnsr_obs@data[total_miss_index] <- NA

# perform the algorithm for tensor completion with missing fibers not at random
res <- tnsr_FMAR(tnsr_obs@data, tnsr_mask@data, missing_mode=4, link = 'probit', 
                 rank=c(r1,r2,r3,r4), rank_prop = c(r1,r2,r3), alpha = max(tnsr_s_prop@data), alpha_prop=max(tnsr_s@data))

# compute the mse for estimation
sum((res$theta - tnsr_s@data)^2)/prod(d1,d2,d3,d4) # mse = 0.003645591


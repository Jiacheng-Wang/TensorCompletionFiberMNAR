#  #######################################################################
#       File-Name:      TensorCompletion_FiberMNAR_main.R
#       Date:           Wed Apr 27 15:15:43 2022
#       Author:         JCW
#       Purpose:        Tensor completion with fibers missing not at random
#       Input Files:    NONE
#       Output Files:   NONE
#       Data Output:    NONE
#       Dependencies:   NONE
#       Status:         In Progress
#  #######################################################################
# load the necessary package
library(rTensor)
library(TensorComplete)
library(pracma)

# main function for tensor missing fiber imputation
#### Input: 
# ttnsr: order 4 array dataset with missing entries as NA
# omega_fiber: order 3 binary array indicating which fiber is observed (1) or missing (0)
# missing_mode: indicating fibers along which mode is missing not at random
# link: choice of link function, candidates are 'logit', 'probit', 'laplacian'
# rank: tucker rank choice for ttnsr
# rank_prop: tucker rank choice for propensity score
# alpha: maximum norm constrain on ttnsr
# alpha_prop: maximum norm constrain on propensity score
#### Output:
# C: core tensor
# A: factor matrices along four modes
# theta: estimate underlying true ttnsr
# iteration: number of iterations to converge
# cost: value for L_ips
tnsr_FMAR <- function(ttnsr, omega_fiber, missing_mode, link = 'logit',
                      rank, rank_prop, alpha = TRUE, alpha_prop = TRUE,epsilon = 1e-4){
  d = dim(ttnsr)
  A_1 = as.matrix(randortho(d[1])[,1:rank[1]])
  A_2 = as.matrix(randortho(d[2])[,1:rank[2]])
  A_3 = as.matrix(randortho(d[3])[,1:rank[3]])
  A_4 = as.matrix(randortho(d[4])[,1:rank[4]])
  C = rand_tensor(modes = rank)
  C = C*ifelse(is.logical(alpha),1,
               1/max(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data))*alpha)
  ## initial theta is in the interior
  if(is.logical(alpha)){
    alpha_minus=alpha_minus2=TRUE
  }else{
    alpha_minus=alpha-epsilon
    alpha_minus2=alpha-2*epsilon
    prevtheta <- ttl(C,list(A_1,A_2,A_3,A_4),ms=1:4)@data
    if(max(abs(prevtheta))>alpha){
      warning("the input tensor exceeds the magnitude bound. Perform rescaling on the core tensor...")
      C=C/max(abs(prevtheta))*alpha
    }
  }
  
  # estimate the propensity score
  est_s <- fit_ordinal(omega_fiber, r = c(rank_prop[1], rank_prop[2], rank_prop[3]),omega = 0, alpha = alpha_prop)
  if (link == 'logit'){
    est_p <- as.tensor(logit_link(est_s$theta))
  }else if(link == 'probit'){
    est_p <- as.tensor(probit_link(est_s$theta))
  }else if(linl == 'laplacian'){
    est_p <- as.tensor(laplacian_link(est_s$theta))
  }
  
  message('finish propensity score estimation!')
  # transfer the propensity estimate into an order 4 tensor
  ips <- array(rep(NA,prod(d[1],d[2],d[3],d[4])),dim = c(d[1],d[2],d[3],d[4]))
  for (i in 1:d[missing_mode]){
    if (missing_mode == 1){
        ips[i,,,] <- est_p@data
    }else if(missing_mode == 2){
      ips[,i,,] <- est_p@data
    }else if(missing_mode == 3){
      ips[,,i,] <- est_p@data
    }else{
      ips[,,,i] <- est_p@data
    }
  }
  
  result = list()
  error<- 3
  iter = 0
  cost=NULL
  while ((error > 10^-4)&(iter<100) ){
    iter = iter +1
    prevtheta <- ttl(C,list(A_1,A_2,A_3,A_4),ms=1:4)@data
    prev <- likelihood(ttnsr,prevtheta,ips)
    
    # update A_1
    W <- kronecker_list(list(A_4, A_3, A_2))%*%t(k_unfold(C,1)@data)
    A_1 <- comb(A_1,W,ttnsr,ips, k=1, alpha= alpha_minus2)
    # orthognalize A_1
    qr_res=qr(A_1)
    A_1=qr.Q(qr_res)
    C=ttm(C,qr.R(qr_res),1)
    
    # update A_2
    W <- kronecker_list(list(A_4, A_3, A_1))%*%t(k_unfold(C,2)@data)
    A_2 <- comb(A_2,W,ttnsr,ips, 2,alpha_minus)
    # orthognalize A_2
    qr_res=qr(A_2)
    A_2=qr.Q(qr_res)
    C=ttm(C,qr.R(qr_res),2)
    
    # update A_3
    W <- kronecker_list(list(A_4, A_2, A_1))%*%t(k_unfold(C,3)@data)
    A_3 <- comb(A_3,W,ttnsr,ips,3,alpha)
    # orthognalize A_3
    qr_res=qr(A_3)
    A_3=qr.Q(qr_res)
    C=ttm(C,qr.R(qr_res),3)
    
    # update A_4
    W <- kronecker_list(list(A_3, A_2, A_1))%*%t(k_unfold(C,4)@data)
    A_4 <- comb(A_4,W,ttnsr,ips,4,alpha)
    # orthognalize A_3
    qr_res=qr(A_4)
    A_4=qr.Q(qr_res)
    C=ttm(C,qr.R(qr_res),4)
    
    # update C
    C_update <- corecomb(A_1,A_2,A_3,A_4, C,ttnsr, ips)
    
    if(!is.logical(alpha) & max(abs(ttl(C_update,list(A_1,A_2,A_3,A_4),ms=1:4)@data))>=alpha_minus2){
      theta <- ttl(C,list(A_1,A_2,A_3,A_4),ms=1:4)@data
      new <- likelihood(ttnsr,theta,ips)
      cost = c(cost,new); break
    }else{
      C=C_update
      theta <- ttl(C,list(A_1,A_2,A_3,A_4),ms=1:4)@data
      new <- likelihood(ttnsr,theta,ips)
      cost = c(cost,new)
      (error <- abs((new-prev)/prev))
    }
    message(paste(iter,"-th  iteration -- L_ips value is", new ," -----------------"))
  }
  result$C <- C; result$A <- list(A_1,A_2,A_3,A_4)
  result$theta= theta
  result$iteration <- iter
  result$cost = cost
  return(result)
}



# dependent functions
# likelihood function
likelihood = function(ttnsr,theta,ips){
  # ttnsr is the tensor observation
  # theta is the true underlying/estimate tensor observation
  # ips is the estimate propensity score
  index=which(is.na(ttnsr)==F & is.na(theta)==F, arr.ind =  T)
  ttnsr=ttnsr[index]
  theta=theta[index]
  ips = ips[index]
  return(sum((ttnsr-theta/sqrt(ips))^2))
}
#### cost function based on unfolded matrix
h1 = function(A_1,W1,ttnsr,ips){
  theta_output =W1%*%c(A_1)
  l=sqrt(sum(((theta_output[is.na(ttnsr)==F]-ttnsr[is.na(ttnsr)==F])^2)/ips[is.na(ttnsr)==F]))
  return(l)
}
#### cost function based on Tucker representation
hc = function(A_1,A_2,A_3,A_4,C,ttnsr,ips){
  theta_output = c(ttl(C,list(A_1,A_2,A_3,A_4),ms=1:4)@data)
  l=sqrt(sum(((theta_output[is.na(ttnsr)==F]-ttnsr[is.na(ttnsr)==F])^2)/ips[is.na(ttnsr)==F]))
  return(l)
}
#### gradient function based on unfolded matrix
g1 = function(A_1,W,ttnsr,ips){
  theta_output =W%*%c(A_1)
  pretheta=theta_output-as.matrix(ttnsr)
  pretheta[is.na(pretheta)==T]=0
  l=(2*t(pretheta)/ips)%*%W
  return(l)
}
#### update one factor matrix at a time while holding others fixed
comb = function(A,W,ttnsr,ips, k, alpha=TRUE){
  nA = A
  tnsr1 <- k_unfold(as.tensor(ttnsr),k)@data
  ips1 <- k_unfold(as.tensor(ips),k)@data
  if (alpha==TRUE) {
    l <- lapply(1:nrow(A),function(i){optim(A[i,],
                                            function(x) h1(x,W,tnsr1[i,],ips1[i,]),
                                            function(x) g1(x,W,tnsr1[i,],ips1[i,]),method = "BFGS")$par})
    nA <- matrix(unlist(l),nrow = nrow(A),byrow = T)
  }else{
    l <- lapply(1:nrow(A),function(i){constrOptim(A[i,],
                                                  function(x) h1(x,W,tnsr1[i,],ips1[i,]),
                                                  function(x) g1(x,W,tnsr1[i,],ips1[i,]),
                                                  ui = as.matrix(rbind(W,-W)),ci = rep(-alpha,2*nrow(W)),method = "BFGS")$par})
    nA <- matrix(unlist(l),nrow = nrow(A),byrow = T)
  }
  return(nA)
}
#### gradient function with respect to theta
gradient_tensor=function(A_1,A_2,A_3,A_4,C,ttnsr,ips){
  theta_output = c(ttl(C,list(A_1,A_2,A_3,A_4),ms=1:4)@data)
  pretheta=theta_output-ttnsr
  pretheta[is.na(pretheta)==T]=0
  l=2*pretheta/ips
  return(l)
}
#### gradient function with respect to the core tensor (gradient not using kronecker product at all)
gc = function(A_1,A_2,A_3,A_4,C,ttnsr,ips){
  g = gradient_tensor(A_1,A_2,A_3,A_4, C,ttnsr,ips) ## first, take entrywise gradient w.r.t. theta
  g = ttl(as.tensor(g),list(t(A_1),t(A_2),t(A_3),t(A_4)),ms = 1:4)@data ## then, take entrywise gradient w.r.t. core tensor
  return(g)
}
#### update core tensor holding other factors fixed
corecomb = function(A_1,A_2,A_3,A_4, C,ttnsr,ips, alpha=TRUE){
  h <- function(x) hc(A_1,A_2,A_3,A_4, new("Tensor",C@num_modes,C@modes,data = x),ttnsr,ips)
  g <- function(x) c(gc(A_1,A_2,A_3,A_4, new("Tensor",C@num_modes,C@modes,data = x),ttnsr,ips))
  d <- optim(c(C@data),h,g,method="BFGS")
  C <- new("Tensor",C@num_modes,C@modes,data =d$par)
  return(C)
}
#### different choice of link function
logit_link <- function(a){
  return(1/(1+exp(-a)))
}
probit_link <- function(a){
  return(pnorm(a))
}
laplacian_link <- function(a){
  if (a<0){
    return(0.5*exp(a))
  }else{
    return(1-0.5*exp(-a))
  }
}


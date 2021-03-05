#==========================================================================================
# Dataset: HouseVotes data (5 replications)
# Distribution: g(m) \propto (\alpha^m) 
# Component: 30 (# of components can be adjusted with 'k' parameter)
# Mixture tuning (concentration) parameter: 20 (adjusted with '\rho')
# Iterations: 1e4 (adjusted with 'nreps')
#==========================================================================================
library(tidyverse);library(phyclust);library(Rcpp)
#==========================================================================================
setwd("/work/STAT/abhishek/HouseVotes")
train_original=readRDS("HouseVotes.rds")
#==========================================================================================
# Define the functions required for the algorithm.
#==========================================================================================
f=function(alpha,x,mu){(alpha^(sum(abs(x-mu))))/(((1-(alpha^(p+1)))/(1-alpha))*choose(p,sum(abs(x-mu))))}
repl=function(mu,n,dimen){s=sort(sample(1:dimen,n));replace(mu,s,1-mu[s])}
Jeffreys_prior=function(alpha,p)
{
  t1=((p+1)^2)*(alpha^p)
  t2=2*p*(p+2)*(alpha^(p+1))
  t3=((p+1)^2)*(alpha^(p+2))
  t4=alpha^((2*p)+2)
  t5=(1-(alpha^(p+1)))^2
  t6=(1-alpha)^2
  n=(1-t1+t2-t3+t4)
  num=ifelse(n<0,0,n)
  den=alpha*t5*t6
  prior=sqrt(num/den)
  return(prior)
}
func_mu=function(dat,mu_proposed,mu_current,theta)
{
  N=nrow(dat);m_proposed=rep(0,N);m_current=rep(0,N)
  for(i in 1:N){m_proposed[i]=sum(abs(dat[i,1:p]-mu_proposed));m_current[i]=sum(abs(dat[i,1:p]-mu_current))}
  return(log_post_density_ratio=((sum(m_proposed)-sum(m_current))*log(exp(theta)/(1+exp(theta))))-
           sum(log(choose(p,m_proposed)))+sum(log(choose(p,m_current))))
}
func_theta=function(dat,theta_proposed,theta_current,mu)
{
  N=nrow(dat);m=rep(0,N);for(i in 1:N){m[i]=sum(abs(dat[i,1:p]-mu))};M=sum(m)
  k1=exp(theta_current);k2=exp(theta_proposed)
  r1=k1/(1+k1);r2=k2/(1+k2)
  term1=M*log(r2/r1)
  term2=log((k2/((1+k2)^2))/(k1/((1+k1)^2)))
  term3=log(Jeffreys_prior(r2,p)/Jeffreys_prior(r1,p))
  term4=N*log(((1-(r2^(p+1)))/(1-r2))/((1-(r1^(p+1)))/(1-r1)))
  return(log_post_density_ratio=term1+term2+term3-term4)
}
update_mu=function(reps)
{
  for(j in 1:k)
  {
    index=(((j-1)*p)+1):(j*p)
    mu_star=repl(mu_mat[reps-1,index],l,p)
    mu_star_mat[reps,index]=mu_star
    log_acceptance_ratio=ifelse(nrow(subset(train,train[,(p+1)]==j))==0,-Inf,
                                min(0,func_mu(subset(train,train[,(p+1)]==j),mu_star,mu_mat[reps-1,index],theta_mat[reps-1,j])))
    u=runif(1,0,1)
    if(log(u)<log_acceptance_ratio){mu_mat[reps,index]=mu_star;acc_mu_mat[reps,j]=1}else{mu_mat[reps,index]=mu_mat[reps-1,index]}
  }
  return(list(mu_mat=mu_mat[reps,],mu_star_mat=mu_star_mat[reps,],acc_mu_mat=acc_mu_mat[reps,]))
}
update_theta=function(reps)
{
  for(j in 1:k)
  {
    index=(((j-1)*p)+1):(j*p)
    theta_star=rnorm(1,theta_mat[reps-1,j],sigma)
    theta_star_mat[reps,j]=theta_star
    log_acceptance_ratio=ifelse(nrow(subset(train,train[,(p+1)]==j))==0,-Inf,
                                min(0,ifelse(is.nan(func_theta(subset(train,train[,(p+1)]==j),
                                    theta_star,theta_mat[reps-1,j],mu_mat[reps,index]))==TRUE,0,
      func_theta(subset(train,train[,(p+1)]==j),theta_star,theta_mat[reps-1,j],mu_mat[reps,index]))))
    u=runif(1,0,1)
    if(log(u)<log_acceptance_ratio){theta_mat[reps,j]=theta_star;acc_theta_mat[reps,j]=1}else{theta_mat[reps,j]=theta_mat[reps-1,j]}
  }
  return(list(theta_mat=theta_mat[reps,],theta_star_mat=theta_star_mat[reps,],acc_theta_mat=acc_theta_mat[reps,]))
}
update_id=function(reps)
{
  alphas=exp(theta_mat[reps,])/(1+exp(theta_mat[reps,]))
  ids=idmat[,reps-1]
  for(i in 1:N)
  {
    fs=rep(0,k);props_num=rep(0,k);props=rep(0,k)
    for(j in 1:k)
    {
      index=(((j-1)*p)+1):(j*p)
      fs[j]=ifelse(is.nan(f(alphas[j],train[i,1:p],mu_mat[reps,index]))==TRUE,0,
                   f(alphas[j],train[i,1:p],mu_mat[reps,index]))
      props_num[j]=((sum(ids[-i]==j)+(rho/k))*fs[j])
    }
    if (sum(props_num)==0) {props=rep(1/k,k)} else {props=props_num/sum(props_num)}
    ids[i]=sample(c(1:k),1,prob=props)
  }
  return(ids)
}
replaceMissing=function(x,obsId,reps)
{
  j=idmat[obsId,reps];index=(((j-1)*p)+1):(j*p)
  alphas=exp(theta_mat[reps,])/(1+exp(theta_mat[reps,]))
  q=sum(is.na(x))
  missIndex=which(is.na(x)==TRUE);m=sum(abs(x[-missIndex]-mu_mat[reps,index][-missIndex]))
  s=sum((choose(q,0:q)*(alphas[j]^(0:q)))/choose(p,(0+m):(q+m)))
  prob=((choose(q,0:q)*(alphas[j]^(0:q)))/choose(p,(0+m):(q+m)))/s
  ind=sample(0:q,1,prob = prob)
  x[missIndex]=repl(mu_mat[reps,index][missIndex],ind,q)
  return(x)
}
cppFunction('NumericMatrix arr(IntegerMatrix m, int nreps, Function f){
            NumericMatrix mat(nreps+1,nreps+1);
            for(int i=0; i<mat.nrow()-1; i++){
            for(int j=(i+1); j<mat.ncol(); j++){
            IntegerVector v1 = m(_,i), v2 = m(_,j);
            List out = f(v1,v2);
            mat(i,j)=out["Rand"];
            }
            }
            return mat;
            }')
#==========================================================================================
for(countData in 1:5)
{
  #==========================================================================================
  # Define the vectors and matrices required for storage.
  #==========================================================================================
  p=ncol(train_original);N=nrow(train_original);k=30
  sigma=0.5;rho=20;l=1                     
  nreps=1e4       #change as needed
  idmat=matrix(nrow=N,ncol=nreps+1)
  mu_mat=matrix(nrow=nreps+1,ncol=p*k);mu_star_mat=matrix(nrow=nreps+1,ncol=p*k)
  theta_mat=matrix(nrow=nreps+1,ncol=k);theta_star_mat=matrix(nrow=nreps+1,ncol=k)
  acc_mu_mat=matrix(0,nrow=nreps+1,ncol=k);acc_theta_mat=matrix(0,nrow=nreps+1,ncol=k)
  #==========================================================================================
  # Initialize the vectors and matrices.
  #==========================================================================================
  idmat[,1]=sample(c(1:k),N,replace = TRUE,prob=rep(1/k,k))
  mu_mat[1,]=sample(c(0,1),p*k,replace = TRUE);mu_star_mat[1,]=mu_mat[1,]
  theta_mat[1,]=rnorm(k,0,sigma);theta_star_mat[1,]=theta_mat[1,]
  na_index=which(is.na(t(train_original)),TRUE)
  na_index=na_index[,c(2,1)];colnames(na_index)=colnames(na_index)[c(2,1)]
  unobs_mat=cbind(na_index,matrix(nrow=nrow(na_index),ncol=nreps+1))
  train=matrix(nrow=N,ncol=p)
  for(i in 1:N){train[i,]=as.vector(as.matrix(replaceMissing(train_original[i,],obsId=i,reps=1)))}
  for(i in 1:nrow(unobs_mat)){unobs_mat[i,3]=train[unobs_mat[i,1],unobs_mat[i,2]]}
  train=cbind(train,id=idmat[,1])
  #==========================================================================================
  # The sampler.
  #==========================================================================================
  for(reps in 2:(nreps+1))
  {
    mus=update_mu(reps)
    mu_mat[reps,]=mus$mu_mat;mu_star_mat[reps,]=mus$mu_star_mat;acc_mu_mat[reps,]=mus$acc_mu_mat
    
    thetas=update_theta(reps)
    theta_mat[reps,]=thetas$theta_mat;theta_star_mat[reps,]=thetas$theta_star_mat;acc_theta_mat[reps,]=thetas$acc_theta_mat
    
    idmat[,reps]=update_id(reps)
    
    train=matrix(nrow=N,ncol=p)
    for(i in 1:N){train[i,]=as.vector(as.matrix(replaceMissing(x=train_original[i,],obsId=i,reps)))}
    for(i in 1:nrow(unobs_mat)){unobs_mat[i,reps+2]=train[unobs_mat[i,1],unobs_mat[i,2]]}
    
    train=cbind(train,id=idmat[,reps])
  }
  #==========================================================================================
  # Obtain results for 'mu's.
  #==========================================================================================
  all_mus=mu_mat[,1:p];for(j in 2:k){index=(((j-1)*p)+1):(j*p);all_mus=rbind(all_mus,mu_mat[,index])}
  all_mus_dist=data.frame(all_mus) %>% group_by_all() %>% summarize(n=n()) %>% arrange(desc(n))
  unique_mus=all_mus_dist %>% select(-n) %>% as.matrix
  
  mus_dist=list()
  for(j in 1:k)
  {
    index=(((j-1)*p)+1):(j*p)
    assign(sprintf("mu%s",j),data.frame(mu_mat[,index]))
    assign(sprintf("mu%s_dist",j),eval(parse(text=sprintf("mu%s",j))) %>% group_by_all() %>% summarise(n=n()))
    mus_dist[[j]]=eval(parse(text=sprintf("mu%s_dist",j)))
  }
  #==========================================================================================
  # Obtain results for 'alpha's.
  #==========================================================================================
  alpha_mat=exp(theta_mat)/(1+exp(theta_mat));alphas=list();alphas_summary=list()
  for(i in 1:nrow(unique_mus))
  {
    a=alpha_mat[which(apply(mu_mat[,1:p],1,function(x) identical(x,as.vector(unique_mus[i,])))),1]
    for(j in 2:k)
    {
      index=(((j-1)*p)+1):(j*p)
      a=c(a,alpha_mat[which(apply(mu_mat[,index],1,function(x) identical(x,as.vector(unique_mus[i,])))),j])
    }
    alphas[[i]]=a;alphas_summary[[i]]=summary(a)
  }
  #==========================================================================================
  # Obtain results for predicted components and predicted vectors.
  #==========================================================================================
  unobs_pred=data.frame(cbind(unobs_mat[,1:2],apply(unobs_mat[,-(1:2)],1,mean)))
  unobs_pred$predicted=ifelse(unobs_pred$V3>0.5,1,0)
  
  train_pred=train_original
  for(i in 1:nrow(unobs_pred)){train_pred[unobs_pred$row[i],unobs_pred$col[i]]=unobs_pred$predicted[i]}
  obsNAindex=which(!(apply(train_original,1,function(x) sum(is.na(x)))==0));obsMatch=rep(NA,length(obsNAindex))
  #==========================================================================================
  # Obtain results for density estimates and compare them.
  #==========================================================================================
  data_dist=data.frame(train_original) %>% group_by_all() %>% summarize(n=n())
  data_dist_complete_index=which((apply(data_dist,1,function(x) sum(is.na(x)))==0))
  interested_xs=data_dist[data_dist_complete_index,] %>% arrange(desc(n)) %>% select(-c(n)) %>% as.matrix
  posterior_rel_freq_ids=matrix(nrow=nreps+1,ncol=k)
  for(i in 1:nrow(posterior_rel_freq_ids))
  {
    for(j in 1:ncol(posterior_rel_freq_ids))
    {
      posterior_rel_freq_ids[i,j]=sum(idmat[,i]==j)/N    
    }
  }
  B=0;estimated_probs=rep(NA,nrow(interested_xs))
  for(i in 1:nrow(interested_xs))
  {
    d_reps=rep(NA,(nreps+1)-(1+B)+1)
    for(reps in (1+B):(nreps+1))
    {
      d=rep(NA,k)
      for(j in 1:k){index=(((j-1)*p)+1):(j*p);d[j]=f(ifelse(alpha_mat[reps,j]==1,(1-1e-7),alpha_mat[reps,j]),
                                                     interested_xs[i,],mu_mat[reps,index])}
      d_reps[reps-B]=sum(d*posterior_rel_freq_ids[reps,],na.rm=TRUE)
    }
    estimated_probs[i]=mean(d_reps)
  }
  complete_data_dist=data_dist[data_dist_complete_index,] %>% arrange(desc(n)) %>% mutate(RelFreq=n/N)
  complete_data_dist$EstProb=estimated_probs
  #==========================================================================================
  mr=arr(idmat,nreps,RRand)
  mr[lower.tri(mr)]=t(mr)[lower.tri(mr)];diag(mr)=NA
  avgRI=apply(mr,2,function(x) mean(x,na.rm = TRUE))
  #==========================================================================================
  # Obtain matrix of posterior probabilities that observations come from the same component 
  # and corresponding heatmap
  #==========================================================================================
  library(mlbench)
  data("HouseVotes84")
  
  idmat_rep=idmat[which(HouseVotes84$Class=="republican"),]
  idmat_dem=idmat[which(HouseVotes84$Class=="democrat"),]
  idmat_all=rbind(idmat_rep,idmat_dem)
  coin_mat=matrix(NA,nrow=nrow(idmat_all),ncol=nrow(idmat_all))
  
  for(i in 1:nrow(coin_mat))
  {
    for(j in 1:nrow(coin_mat))
    {
      coin_mat[i,j]=sum(idmat_all[i,]==idmat_all[j,])/(nreps+1)
    }
  }
  
  library(RColorBrewer)
  heatmap(coin_mat,Colv = NA, Rowv = NA)
  #==========================================================================================
  rm(mr)
  save.image(file=sprintf("HouseVotes%s_v2.RData",countData))
  #==========================================================================================
  gdata::keep(train_original,countData,f,repl,Jeffreys_prior,func_mu,func_theta,
              update_mu,update_theta,update_id,replaceMissing,arr,sure=TRUE)
  #==========================================================================================
}
#==========================================================================================
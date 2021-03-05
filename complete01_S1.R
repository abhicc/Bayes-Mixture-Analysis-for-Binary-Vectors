#==========================================================================================
# Dataset: Simulated complete dataset . This is the MASTER code (1st sampler)
# Distribution: g(m) \propto (\alpha^m) 
# Component: k (# of components can be adjusted with 'k' parameter)
# Mixture tuning (concentration) parameter: \rho (adjusted with '\rho')
# Iterations: nreps (adjusted with 'nreps')
#==========================================================================================
library(tidyverse);library(LaplacesDemon)
#==========================================================================================
# Generate the simulated dataset. Change 'N','p','k_true','mu's, and 'alpha's as needed.
# NOTE: As 'k_true' changes, number of 'mu's and 'alpha's change accordingly.
#==========================================================================================
N=500;p=6;k_true=2     #change as needed, and change number of 'mu's and 'alpha's with k
f=function(alpha,m){(alpha^m)/((1-(alpha^(p+1)))/(1-alpha))}
repl=function(mu,n,dimen){s=sort(sample(1:dimen,n));replace(mu,s,1-mu[s])}
weight=rep(1/k_true,k_true);ids_true=sample(c(1:k_true),N,replace = TRUE,prob=weight)
mus_true=matrix(nrow=k_true,ncol=p)
for(i in 1:k_true){mus_true[i,]=sample(c(0,1),p,replace = TRUE);assign(sprintf("mu%s_true",i),mus_true[i,])}
alphas_true=round(runif(k_true,0,1),3);for(i in 1:k_true){assign(sprintf("alpha%s_true",i),alphas_true[i])}
dat=matrix(nrow=N,ncol=p+2);dat[,1]=ids_true
for(i in 1:N)
{
  dat[i,2]=sample(0:p,1,prob = f(eval(parse(text=sprintf("alpha%s_true",dat[i,1]))),0:p))
  dat[i,3:(p+2)]=repl(eval(parse(text=sprintf("mu%s_true",dat[i,1]))),dat[i,2],p)
}
train_original=dat[,3:(p+2)];dat=data.frame(dat);colnames(dat)[1:2]=c("id","m")
train_dist=data.frame(train_original) %>% group_by_all() %>% summarize(n=n())
m=rep(NA,nrow(train_dist))
for(i in 1:nrow(train_dist)){m[i]=sum(abs(as.vector(as.matrix(train_dist[i,1:p]))-as.vector(as.matrix(train_dist[which.max(train_dist$n),1:p]))))}
train_dist=cbind(train_dist,m=m)      
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
#==========================================================================================
# Define the vectors and matrices required for storage.
#==========================================================================================
p=ncol(train_original);N=nrow(train_original);k=5
sigma=0.5;rho=20;l=1
nreps=1000       #change as needed
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
train=train_original;train=cbind(train,id=idmat[,1])
#==========================================================================================
# The sampler.
#==========================================================================================
for(reps in 2:(nreps+1))
{
  mus=update_mu(reps)
  mu_mat[reps,]=mus$mu_mat;mu_star_mat[reps,]=mus$mu_star_mat;acc_mu_mat[reps,]=mus$acc_mu_mat
  
  thetas=update_theta(reps)
  theta_mat[reps,]=thetas$theta_mat;theta_star_mat[reps,]=thetas$theta_star_mat;acc_theta_mat[reps,]=thetas$acc_theta_mat
  
  train=subset(train,select=-id)
  
  idmat[,reps]=update_id(reps)

  train=cbind(train,id=idmat[,reps])
}
#==========================================================================================
# Obtain results for 'mu's.
#==========================================================================================
apply(acc_mu_mat,2,function(x) prop.table(table(x))[2]*100)

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
apply(acc_theta_mat,2,function(x) prop.table(table(x))[2]*100)

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
# Obtain results for the number of non-empty clusters.
#==========================================================================================
table(apply(idmat,2,function(x) length(unique(x))))
#==========================================================================================
# Obtain results for density estimates and compare them.
#==========================================================================================
n=as.vector(prop.table(table(dat[,1])));weight=rep(1/k_true,k_true)
interested_xs=train_dist %>% select(-c(n,m)) %>% as.matrix
f=function(alpha,x,mu){(alpha^(sum(abs(x-mu))))/(((1-(alpha^(p+1)))/(1-alpha))*choose(p,sum(abs(x-mu))))}
mixture_probs_n=rep(NA,nrow(interested_xs));mixture_probs_w=rep(NA,nrow(interested_xs))
for(i in 1:nrow(interested_xs))
{
  mp_n=rep(NA,k_true);mp_w=rep(NA,k_true)
  for(j in 1:k_true)
  {
    mp_n[j]=n[j]*f(eval(parse(text=sprintf("alpha%s_true",j))),interested_xs[i,],eval(parse(text=sprintf("mu%s_true",j))))
    mp_w[j]=weight[j]*f(eval(parse(text=sprintf("alpha%s_true",j))),interested_xs[i,],eval(parse(text=sprintf("mu%s_true",j))))
  }
  mixture_probs_n[i]=sum(mp_n);mixture_probs_w[i]=sum(mp_w)
}
mixture_probs=cbind(mixture_probs_n,mixture_probs_w)
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
    d_reps[reps-B]=sum(d*posterior_rel_freq_ids[reps,])
  }
  estimated_probs[i]=mean(d_reps)
}
cbind(estimated_probs=round(estimated_probs,3),round(mixture_probs,3))
sum(mixture_probs[,1]*log(mixture_probs[,1]/estimated_probs))    #KL divergence with relative frequencies
sum(mixture_probs[,2]*log(mixture_probs[,2]/estimated_probs))    #KL divergence with true weights
#==========================================================================================
# Obtain matrix of posterior probabilities that observations come from the same component 
# and corresponding heatmap
#==========================================================================================
idmat1=idmat[which(dat$id==1),]    # since we started with k_true=2
idmat2=idmat[which(dat$id==2),]
idmat_all=rbind(idmat1,idmat2)
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
save.image(file="complete01_S1.RData")
#==========================================================================================
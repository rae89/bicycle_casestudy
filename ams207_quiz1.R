#case study dataset
bikes<- c(16,9,10,13,19,20,18,17,35,55)
cars<- c(58,90,48,57,103,57,86,112,273,64)
tot<- bikes+cars
dat= cbind(bikes,cars,tot)
n=dat[,3]
y=dat[,1]


prop=y/n

########### part (a) #####################
loghyperpost = function(u, v)
{ 
  a = (exp(u) * exp(v))/(exp(u) + 1)
  b = exp(v)/(exp(u) + 1)
  
  J = length(dat[,1])
  x = J * (lgamma(a + b) - lgamma(a) - lgamma(b)) + log(a) + log(b) - 
    2.5 * log(a + b)
  for(i in (1:J)) {
    x <- x + lgamma(a + y[i]) + lgamma(b + n[i] - y[i]) -
      lgamma(a + b + n[i]) }
  x }
rpost=function(n){
  i=sample(1:length(pc),size=n,replace=T,prob=pc)
  iu=c(row(p))[i];iv=c(col(p))[i]
  u=u.grid[iu];v=v.grid[iv]
  return(cbind(jitter(u),jitter(v)))}

u.grid = seq(-2.3, -.5, length = 50) # log(alpha/beta)
v.grid = seq(0, 5, length = 50) # log(alpha+beta)

p=matrix(0,nrow=50,ncol=50) # initialize matrix for posterior
# eval’s over agrid x bgrid
for (i in 1:50){
  p[i,]=loghyperpost(u.grid[i],v.grid)
}
pc <- exp(c(p)) # all cell p’s in one vector
image(u.grid,v.grid,exp(p),xlab='log(alpha/beta)',ylab='log(alpha+beta)')
contour(u.grid,v.grid,exp(p),xlab='log(alpha/beta)',ylab='log(alpha+beta)',add=TRUE)

normp=exp(p)/sum(exp(p))
p.ugrid=apply(normp,1,sum)

plot(u.grid,p.ugrid,ylab="",xlab='log(alpha/beta)',type='l')
plot(v.grid,p.ugrid,ylab="",xlab='log(alpha+beta)',type='l')

# generate a posterior sample
hyper <- rpost(500) # simulate 500 draws from posterior
plot(hyper[,1],hyper[,2],xlab="log(alpha/beta)",ylab="log(alpha+beta",
     xlim=c(-2.3,-0.5),ylim=c(0,5),pch="o")
contour(u.grid,v.grid,normp,xlab="log(alpha/beta)",ylab="log(alpha+beta")


#transform back to alpha and beta
alpha = (exp(hyper[,1])*exp(hyper[,2]))/(exp(hyper[,1])+1)
beta = exp(hyper[,2])/(exp(hyper[,1])+1)

#simulate theta 
theta = c()
theta = rbeta(500,alpha+y,beta+n-y)

theta.ci= quantile(theta,c(.025,.975))
#     2.5%      97.5% 
#  0.08093801 0.45221807

plot(theta)
# simulate new block of 100 vehicles 
y.pred = rbinom(500,100,theta)
mean(y.pred)
#[1] 19.684
y.pred.ci= quantile(y.pred,c(.025,.975))
#   2.5%  97.5% 
#  7.000 46.525 

################### Part (b) #########################

# simualate theta for each block, which gives us a theta of n by 10,
# where n is the number of simulations. 
J = length(dat[,1])

theta2 = matrix(0, nrow=500, ncol=10)
for(j in 1:J){
  theta2[,j] = rbeta(500, alpha+y[j], beta+n[j]-y[j])
}

#posterior intervals for the thetas of each block
for(j in 1:J)
{
  print (quantile(theta2[,j],c(.025,.975)))
}
"
theat 1:              2.5%     97.5% 
                  0.1369373 0.3040224 
theta 2:              2.5%      97.5% 
                  0.05629672 0.17705251 
theta 3:              2.5%     97.5% 
                0.1059840 0.2746341 
theta 4:              2.5%     97.5% 
                0.1176816 0.2860838 
theta 5:              2.5%     97.5% 
                0.1048029 0.2201846 
theta 6:              2.5%     97.5% 
                0.1688393 0.3503182 
theta 7:              2.5%     97.5% 
                0.1150592 0.2440149 
theta 8:              2.5%      97.5% 
                0.08393317 0.20320685 
theta 9:              2.5%      97.5% 
                0.08547295 0.15529631 
theta 10:              2.5%     97.5% 
                0.3427029 0.5252146 
"
colMeans(theta2)
#  [1] 0.2129633 0.1065583 0.1794642 0.1912028 0.1606171 0.2506573 0.1779976 0.1393867
#  [9] 0.1187810 0.4292091


############## part(c) ################

#replications of model 1
y1.rep = rbinom(10,n,sample(theta,10,replace=F))

#replications of model 2
y2.rep = c()
for(j in 1:J){
  y2.rep[j] = rbinom(1,n[j],mean(theta2[,j]))
}
hist(y)
hist(y1.rep)
hist(y2.rep)

t.test(y1.rep, y, paired = TRUE) 
"
Paired t-test

data:  y1.rep and y
t = -0.5098, df = 9, p-value = 0.6224
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -14.679944   9.279944
sample estimates:
  mean of the differences 
-2.7 
"
t.test(y2.rep, y, paired = TRUE) 
"
  Paired t-test

data:  y2.rep and y
t = -0.3166, df = 9, p-value = 0.7588
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -4.072808  3.072808
sample estimates:
mean of the differences 
                   -0.5 
"


library(schoRsch)
library(xtable)

T.Test <- t.test(y1.rep, y, paired = TRUE) 
xtable(
  t_out(toutput=T.Test, n.equal = TRUE, welch.df.exact = TRUE, welch.n = NA,
        d.corr = TRUE, print = TRUE)
)

T.Test <- t.test(y2.rep, y, paired = TRUE) 
xtable(
  t_out(toutput=T.Test, n.equal = TRUE, welch.df.exact = TRUE, welch.n = NA,
        d.corr = TRUE, print = TRUE)
)






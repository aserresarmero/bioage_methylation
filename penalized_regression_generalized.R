library(ape);library(nlme);library(glmnet);library(reshape2);library(MASS);library(parallel)

####Implemented from: Mylona, K. and Goos, P. (2011) Penalized Generalized 
####Least Squares for Model Selection under Restricted Randomization.

### Most functions take in: X (matrix containing variables as rows)
### V (variance-covariance matrix)
### Y (response variable as vector)
### bet (penalized regression beta coefficients)
### lam (penalized regression lambdas)
### ... Or some name variation of these parameters

## Antoniadis (1997)
lasso2=function(bet,lam){
  lam^2-(abs(bet)-lam)^2*ifelse(abs(bet)<lam,1,0)
}

## Proper lasso
lasso=function(bet,lam){
  abs(bet)*lam
}

###SCAD penalty
lasso3=function(bbet,lam,alf=3.7){
  sapply(bbet,function(bet){
    if(abs(bet)<lam){abs(bet)*lam}else if(abs(bet)>lam & abs(bet)<lam*alf){
      ((alf^2-1)*(lam^2)-(abs(bet)-alf*lam)^2)/(2*alf-1)}else{
        1/2*(alf+1)*lam^2
      }
  })
}

### Make W matrix from paper 
makeW=function(bet,lam){
   length(bet)*diag(as.vector(lasso(bet,lam)/abs(bet)))
}

### Beta update from paper 
upbet=function(X,V,Y,W){
  ginv(t(X)%*%V%*%X+W)%*%t(X)%*%V%*%Y
}

### Cross validation error from paper
cob=function(X,V,Y,W){
  tm=ginv(t(X)%*%V%*%X+W)
  tm%*%t(X)%*%V%*%X%*%tm
}

### Cross validation error
cv=function(X,V,Y,W){
  ph=X%*%ginv(t(X)%*%V%*%X+W)%*%t(X)%*%V
  sum((Y-ph%*%Y)^2)/dim(V)[1]/(1-sum(diag(ph))/dim(V)[1])^2
}

#### Iterative upgrading of beta coefficients
optbet=function(X,V,Y,lam,plo=T){
  ##Initial beta
  bet=ginv(t(X)%*%V%*%X)%*%t(X)%*%V%*%Y
  W=1*diag(rep(1,length(bet)))
  m=c()
  posit=(1:length(bet))
  ## One iteration is enough for this proof example, we do 4
  for(i in 1:4){
    ## Upbet just adds the penalty as described in paper
    nb=upbet(X,V,Y,W)
    m=c(m,summary(lm(X%*%nb~Y))$"r.squared")
    ## Manually filter idle coefficients (bad idea for proper implementation)
    kee=abs(nb)>1/100
    nb=nb[kee]
    posit=posit[kee]
    ## Update W 
    W=makeW(nb,lam)
    X=X[,kee]
  }
  if(plo==T){plot(m)}
  bets=rep(0,length(bet));bets[posit]=nb
  return(bets)
}

#### Cross validation scheme
cvfun=function(X,V,Y,lamarr=10^(seq(-2,0.5,.3)),folds=4){
  mV=ginv(V)
  di=dim(mV)[1]/folds
  sizes=rep(floor(di),folds);sizes[1]=sizes[1]+dim(V)[1]%%folds
  ### Construct cross-validation folds
  smp=sample(unlist(sapply(1:length(sizes),function(x)rep(x,sizes[x]))))
  lX=tapply(1:dim(mV)[1],smp,function(x)X[x,])
  lV=tapply(1:dim(mV)[1],smp,function(x)ginv(V[x,x]))
  lVV=tapply(1:dim(mV)[1],smp,function(x)ginv(V[-x,-x]))
  lY=tapply(1:dim(mV)[1],smp,function(x)Y[x])
  ## Optimize for each lambda independently.
  op=mclapply(lamarr,function(l){
    tryCatch({
      ## For any ii fold 
      obb=optbet(X,V,Y,l,F)
      thing=sapply(1:folds,function(ii){
        ## Separate into trainning set
        a=do.call(rbind,lX[-ii]);b=lVV[[ii]];cz=do.call(c,lY[-ii]);
        ## Run the beta optimization for training set
        ob=optbet(a,b,cz,l,F)
        ## Cross validation error on test set together with betas
        c(cv(lX[[ii]][,ob!=0],lV[[ii]],lY[[ii]],makeW(ob[ob!=0],l)),obb)
      })
      ## Make CV interval and report best fold betas
      return(list(thing[1,],thing[-1,which.min(thing[1,])]))
    },error=function(e){return(list(rep(NA,folds),rep(0,dim(X)[2])))})
  })
  return(list(sapply(op,function(x)x[[2]]),sapply(op,function(x)x[[1]]),lamarr))
}


set.seed(123)
##Number of samples
n=300
##Number of sites per matrix
s=200
## Generate true regression weights
w=2*sample(c(-1,1),s,replace = T)*10**(seq(-5,1,length.out=s))

###Generate random phylogeny that will be contaminated by error 
rtr=rtree(n)

##Generate random phylogeny that won't be contaminated by error 
rtr3=rtree(n)

## Create matrix 1 from a phylogeny
x=mvrnorm(s,rep(0,n),vcv(rtr3))
## Create matrix 2 phylogenetically independent
x2=t(sapply(1:(s),function(x)rnorm(n)))
## Create matrix 3 (same vcv as error)
x3=mvrnorm(s,rep(0,n),vcv(rtr))
xf=rbind(x,x3,x2)

## Spread error across all matrices
e=mvrnorm(s*3,rep(0,n),vcv(rtr))

## Make response variable by adding up w*XF+e
y=apply(sweep(xf,1,rep(w,3),"*")+e,2,sum)

## Rescale variables (unnecessary)
X=cbind(rep(1,n),t(t((apply(xf,1,scale)))))
Y=scale(y)
V=vcv(rtr)

## Estimate regression weights
bets2=cvfun(X,V,Y,lamarr=exp(c(seq(-4,2,0.5))),folds=4)
boxplot(log(melt(bets2[[2]])[,3])~melt(bets2[[2]])[,2],axes=F,xlab="lambda",ylab="CV error")
axis(1,1:length(bets2[[3]]),labels = round(bets2[[3]],4));axis(2)
bets=bets2[[1]][,which.min(apply(bets2[[2]],2,median))]
cho=sort(which(!apply(bets2[[2]],2,mean)>min(apply(bets2[[2]],2,mean),na.rm=T)+0.1),decreasing = T)[1]
bets=bets2[[1]][,cho]
apply(bets2[[1]],2,function(x)sum(x!=0))
W=makeW(bets,bets2[[3]][cho])
W[is.na(W)]=0
err=diag(cob(X,ginv(V),Y,W))

### Repeat with glmnet
mod=cv.glmnet(t(xf),Y,alfa=1)
mi=mod$"lambda.1se"
W=makeW(coefficients(mod)[-1],mi)
W[is.na(W)]=0
err2=diag(cob(X[,-1],diag(rep(1,dim(V)[1])),Y,W))

plot(Y,predict(mod, t(xf), s = mi),xlab="True",ylab="Predicted")
plot(Y,X%*%bets,xlab="True",ylab="Predicted")
f=as.data.frame(t(xf));colnames(f)=paste0("V",1:dim(f)[2]);mod2=gls(Y~V111,correlation = corPagel(value = 1,phy=as.phylo(rtr),fixed=T),data=cbind.data.frame(f,Y))
res=recalc(mod2$modelStruct$corStruct,list(Xy = as.matrix((Y-X%*%bets)/sd(Y-X%*%bets))))$Xy[, 1]
qqnorm(Y-predict(mod,t(xf),s=mi));qqline(Y-predict(mod,t(xf),s=mi))
qqnorm(Y-t(xf)%*%bets[-1]);qqline(Y-t(xf)%*%bets[-1])

plot(coefficients(mod,s=mi)[-1],bets[-1],col=heat.colors(10,rev = T)[cut(abs(w),10)],xlab="regular",ylab="phylo")

plot(c(sqrt(abs(err2)),-sqrt(abs(err2))),c(1:(length(w)*3),1:(length(w)*3)),pch=20,cex=.5,xlim=c(-max(sqrt(abs(err2))),max(sqrt(abs(err2)))),col=adjustcolor("green",.3),axes=F,xlab="",ylab="");par(new=T)
plot(rep(w,3),1:(length(w)*3),axes=F,xlab="",ylab="",pch=20,col=adjustcolor("red",.2),xlim=c(-max(abs(w)),max(abs(w))));par(new=T)#;axis(1)
plot(coefficients(mod,s=mi)[-1],1:length(coefficients(mod,s=mi)[-1]),pch=21,cex=1.2,axes=F,xlab="",ylab="",xlim=c(-max(abs(coefficients(mod,s=mi)[-1])),max(abs(coefficients(mod,s=mi)[-1]))),bg=rep(c("orange","purple","black"),each=s));axis(1)
legend("bottomright",fill = c("green","red","white"),legend = c("Error","True","Beta"),cex=.7)

#Plot numeric implementation
pdf(file  = "Rplot2.pdf",height = 6,width = 10)
layout(cbind(matrix(rep(1,3),ncol=1),matrix(rep(2,6),ncol=2),3:5))
image(X,axes=F);abline(h=c(1:2/3),lty="dashed",lwd=4)
text(x=rep(.5,3),y=c(.15,.5,.85),c("C","B","A"),cex=2)
plot(c(sqrt(abs(err[-1])),-sqrt(abs(err[-1]))),c(1:(length(err[-1])),1:(length(err[-1]))),pch=20,cex=.5,xlim=c(-max(sqrt(abs(err[-1]))),max(sqrt(abs(err[-1])))),col=adjustcolor("green",.3),axes=F,xlab="Coefficients",ylab="");par(new=T)
plot(rep(w,3),1:(length(w)*3),axes=F,xlab="",ylab="",col=adjustcolor("red",.2),pch=20,cex=.5,xlim=c(-max(abs(w)),max(abs(w))));par(new=T)#;axis(1)
plot(bets[-1],1:length(bets[-1]),axes=F,xlab="",ylab="",pch=21,cex=1.2,xlim=c(-max(abs(bets[-1])),max(abs(bets[-1]))),bg=rep(c("orange","purple","black"),each=s));axis(1)
legend("bottomright",fill = c("green","red","white"),legend = c("Error","True","Beta"),cex=.7)
boxplot(log(melt(bets2[[2]])[,3])~melt(bets2[[2]])[,2],axes=F,xlab="lambda",ylab="CV error")
axis(1,1:length(bets2[[3]]),labels = round(bets2[[3]],4));axis(2)
plot(Y,X%*%bets,xlab="True",ylab="Predicted")
qqnorm(Y-t(xf)%*%bets[-1]);qqline(Y-t(xf)%*%bets[-1])
dev.off()


CC=ginv(t(chol(vcv(rtr))))
#CC=ginv(vcv(rtr))
mod=cv.glmnet(CC%*%t(xf),CC%*%Y,alpha=1)
mi=mod$lambda.1se
mi=0.08
W=makeW(coefficients(mod)[-1],lam=mi)
W[is.na(W)]=0
err3=diag(cob(X[,-1],ginv(V),Y,W))

### Plot glmnet implementation 
pdf(file  = "Rplot.pdf",height = 6,width = 10)
layout(cbind(matrix(rep(1,3),ncol=1),matrix(rep(2,6),ncol=2),3:5))
image(X,axes=F);abline(h=c(1:2/3),lty="dashed",lwd=4)
text(x=rep(.5,3),y=c(.15,.5,.85),c("C","B","A"),cex=2)
plot(c(sqrt(abs(err3)),-sqrt(abs(err3))),c(1:(length(w)*3),1:(length(w)*3)),pch=20,cex=.5,xlim=c(-max(sqrt(abs(err3))),max(sqrt(abs(err3)))),col=adjustcolor("green",.3),axes=F,xlab="Coefficients",ylab="");par(new=T)
plot(rep(w,3),1:(length(w)*3),axes=F,xlab="",ylab="",pch=20,col=adjustcolor("red",.2),xlim=c(-max(abs(w)),max(abs(w))));par(new=T)#;axis(1)
plot(coefficients(mod,s=mi)[-1],1:length(coefficients(mod,s=mi)[-1]),pch=21,cex=1.2,axes=F,xlab="",ylab="",xlim=c(-max(abs(coefficients(mod,s=mi)[-1])),max(abs(coefficients(mod,s=mi)[-1]))),bg=rep(c("orange","purple","black"),each=s));axis(1)
plot(mod)
plot(predict(mod,newx =t(xf),s=mi),Y,xlab="True",ylab="Predicted")
qqnorm(Y-predict(mod,t(xf),s=mi));qqline(Y-predict(mod,t(xf),s=mi))
dev.off()

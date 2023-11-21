library(DeconRNASeq)

### Read in and process methylation Athlas
me=read.csv("METHYLATION_ATLAS")
me=me[!duplicated(me[,1]),]
dec=read.csv("COMMON_CGs",header = T)
ad=read.csv("EXTRA_CGs",header=F);colnames(ad)=colnames(dec)
dec=dec[dec[,1]%in%me[,1],]
dec=rbind(dec,ad)
dec[(is.na(dec))]=0
tmp=unlist(dec[,1])
dec=as.data.frame(dec[,-1])
row.names(dec)=tmp

### Prepare Methylation array data
load("METHYLATION_ARRAY_MATRIX")
l=.csv("METHYLATION_ARRAY_METADATA", header=T)
extra=.table("EXTRA_CGs")
oth=read.csv("EXTRA_CGs2")
x=normalized_betas_sesame
x=x[,!colnames(x)%in%l$Basename[l$CanBeUsedForAgingStudies=="no"]]

### Match methylation array with atlas panel
dmat=x[match(rownames(dec),unlist(x[,1]),nomatch = 0),]

### Load Human methylation data
full=read.table("METHYLATION_MATRIX_HUMAN",header=T,sep="\t")
### Match Human methylation data with atlas panel
full=full[match(me[,1],full[,1],nomatch=0),]
me=me[match(full[,1],me[,1],nomatch=0),]
rful1=full[match(unlist(dmat[,1]),full[,1]),]
rful2=full[match(extra[match(ad[,1],extra[,1],nomatch = 0),3],full[,1]),]
rful2[,1]=extra[match(ad[,1],extra[,1],nomatch = 0),1];rful=rbind(rful1,rful2)
rful=rful[!is.na(rful[,1]),]
row.names(full)=full[,1];full=full[,-1];row.names(me)=me[,1];me=me[,-1];row.names(rful)=rful[,1];rful=rful[,-1]

dmat2=x[match(extra[match(ad[,1],extra[,1],nomatch = 0),2],unlist(x[,1])),]
dmat2[,1]=extra[match(ad[,1],extra[,1],nomatch = 0),1]
dmat=rbind(dmat,dmat2)
tmp=unlist(dmat[,1])
dmat=as.data.frame(dmat[,-1])
row.names(dmat)=tmp
## Deconvolution using standard method
m=DeconRNASeq(dmat,as.data.frame(dec),use.scale = T)

###Read in sheep data
sheep=read.table("METHYLATION_MATRIX_SHEEP",header=T)
## Match sheep data to atlas
dmatsh=sheep[match(rownames(dec),sheep[,1],nomatch = 0),]
dmat2sh=sheep[match(extra[match(ad[,1],extra[,1],nomatch = 0),2],sheep[,1]),]
dmat2sh[,1]=extra[match(ad[,1],extra[,1],nomatch = 0),1]
dmatsh=rbind(dmatsh,dmat2sh)
tmp=unlist(dmatsh[,1])
dmatsh=as.data.frame(dmatsh[,-1])
row.names(dmatsh)=tmp
msh=DeconRNASeq(dmatsh,as.data.frame(dec),use.scale = T)
shtab=read.table("METADATA_SHEEP",header=T)
prsh=prcomp(t(sheep[,-1]))
m2=as.data.frame(matrix(runif(prod(dim(dmat))),ncol=756))
colnames(m2)=colnames(dmat);row.names(m2)=row.names(dmat)

m2=DeconRNASeq(m2,as.data.frame(dec),use.scale = T)

library(quadprog)
## Quadprog function based on https://stackoverflow.com/questions/45577591/linear-regression-with-constraints-on-the-coefficients
decon_blues=function(Y,X,mod=1){
  Rinv <- solve(chol(t(X) %*% as.matrix(X)));
  C <- cbind(rep(1,length(X)), diag(length(X)))
  b <- c(1,rep(0,length(X)))
  sol=apply(Y,2,function(YY){
    d <- t(YY) %*% as.matrix(X)
    n=solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
    if(mod==1){n$solution}else{
      denom=sum(as.vector(Y[,1]-t(n$solution)%*%t(X))^2)/(length(X[,1])-1)
      num=sum((as.vector(t(n$solution)%*%t(X))-mean(t(n$solution)%*%t(X)))^2)/(length(X[,1])-1)
      se=sqrt(diag(mean(as.vector(Y[,1]-t(n$solution)%*%t(X))^2)*chol2inv(chol(t(X)%*%as.matrix(X)))))
      c(pt(-abs(n$solution/se),length(X[,1]))*2,1-pf(num/denom,length(X[1,])-1,length(X[,1])-1))
    }
  })
  if(mod==1){row.names(sol)=colnames(X);colnames(sol)=colnames(Y)}else{row.names(sol)=c(colnames(X),"FullReg");colnames(sol)=colnames(Y)}
  sol
}

m3=as.data.frame(matrix(runif(prod(dim(full))),ncol=490))
colnames(m3)=colnames(full);row.names(m3)=row.names(full)
m3=DeconRNASeq(m3,me)

nu=t(read.table("METHYLATION_MATRIX_LABS",header=T))
ke=read.csv("REFERENCE_ATLAS_AGAIN")
ke=ke[!duplicated(ke[,1]),]

tmp=nu[1,];nu=apply(nu[-1,],1,as.numeric)
row.names(nu)=tmp;nu=as.data.frame(t(nu))
ke=ke[match(row.names(nu),ke[,1]),]
row.names(ke)=ke[,1];ke=ke[,-1]
m4=DeconRNASeq(nu,ke)
fu=as.data.frame(matrix(runif(prod(dim(nu))),nrow=dim(nu)[1]));row.names(fu)=row.names(nu)
m5=DeconRNASeq(fu,ke)

full2=full[match(row.names(ke),row.names(full),nomatch = 1),]
ke=ke[match(row.names(full2),row.names(ke),nomatch = 0),]
nu=nu[match(row.names(nu),row.names(ke),nomatch = 0),]
m6=DeconRNASeq(full2,ke)

pdf("test_decon.pdf",15,10)
par(mfrow=c(2,2),cex.lab=2.5,cex.axis=1.8,cex.main=3,mar=c(4.1,6.1,4.1,2.1))
zzz=DeconRNASeq(full,me)
zzz2=DeconRNASeq(rful,as.data.frame(dec[match(row.names(rful),row.names(dec),nomatch = 0),]))
## Deconvolution using quadprog
sol2=decon_blues(full,as.data.frame(me),mod=2)
plot(-log(sapply(sol2[26,],function(y)ifelse(y<1E-20,1E-16,y)),10),ylim=c(0,16),xlab="Sample",ylab="-log Regression p-value",col=c("black","red")[(cutree(hclust(dist(zzz2$out.all)),k=13)%in%(3:5))+1],main="Human (N=5918)",pch=16,cex=.7);abline(h=-log(0.05),lty="dashed")
sol2=decon_blues(rful,as.data.frame(dec[match(row.names(rful),row.names(dec),nomatch = 0),]),mod=2)
plot(-log(sapply(sol2[26,],function(y)ifelse(y<1E-20,1E-16,y)),10),ylim=c(0,16),xlab="Sample",ylab="-log Regression p-value",main="Human (N=270)",pch=16,cex=.7,col=c("black","red")[(cutree(hclust(dist(zzz2$out.all)),k=13)%in%(3:5))+1]);abline(h=-log(0.05),lty="dashed")
sol2=decon_blues(dmat,dec,2)
plot(-log(sapply(sol2[26,],function(y)ifelse(y<1E-20,1E-16,y)),10),ylim=c(0,16),xlab="Sample",ylab="-log Regression p-value",main="Dog (N=270)",pch=16,col="red",cex=.7);abline(h=-log(0.05),lty="dashed")
sol2=decon_blues(nu,ke,2)
plot(-log(sapply(sol2[26,],function(y)ifelse(y<1E-20,1E-16,y)),10),ylim=c(0,16),xlab="Sample",ylab="-log Regression p-value",main="Dog Lab (N=702)",pch=16,col="red",cex=.7);abline(h=-log(0.05),lty="dashed")
mtext(LETTERS[1:4],side=1,line = rep(c(-58,-28),each=2),outer = T,at=rep(c(0.02,0.52),2),font=2,cex=3)
dev.off()

##Full

library(reshape2)
library(ggplot2)
library(RColorBrewer)


panc=c("Pancreatic_beta_cells","Pancreatic_acinar_cells","Pancreatic_duct_cells")
leuc=c("CD4T.cells_EPIC","NK.cells_EPIC","CD8T.cells_EPIC","B.cells_EPIC")
tis=sapply(strsplit(colnames(rful)[!duplicated(cutree(hclust(dist(m6$out.all)),k=13))],"_"),function(x)x[length(x)])
tis1=paste0(c(paste0(tis,unlist(sapply(rle(tis)$lengths,function(x)1:x))),"r_unif","Labs")," (N=",c(unname(table(cutree(hclust(dist(zzz2$out.all)),k=13))),"96","96"),")")
tis1=sub("colon4","hepatocytes2",tis1)
prop=rbind(do.call(rbind,by(m6$out.all,cutree(hclust(dist(m6$out.all)),k=13),function(x)apply(x,2,mean))),apply(m5$out.all,2,mean),apply(m4$out.all,2,mean))
prop2=rbind(do.call(rbind,by(m6$out.all,cutree(hclust(dist(m6$out.all)),k=13),function(x)apply(x,2,sd))),apply(m5$out.all,2,sd),apply(m4$out.all,2,sd))
prop=cbind.data.frame(prop,Leucocytes_EPIC=unname(apply(prop[,colnames(prop)%in%leuc],1,sum)));prop=prop[,!colnames(prop)%in%leuc]
prop=cbind.data.frame(prop,Pancreas=apply(prop[,colnames(prop)%in%panc],1,sum));prop=prop[,!colnames(prop)%in%panc]
prop2=cbind.data.frame(prop2,Leucocytes_EPIC=apply(prop2[,colnames(prop2)%in%leuc],1,max));prop2=prop2[,!colnames(prop2)%in%leuc]
prop2=cbind.data.frame(prop2,Pancreas=apply(prop2[,colnames(prop2)%in%panc],1,max));prop2=prop2[,!colnames(prop2)%in%panc]
prop=prop[,rev(c(1,2,19,3,8,5,6,4,9,7,13,10,11,14,15,12,18,16,17,20))];prop2=prop2[,rev(c(1,2,19,3,8,5,6,4,9,7,13,10,11,14,15,12,18,16,17,20))]

row.names(prop)=1:length(prop[,1]);prop=melt(as.matrix(prop));
row.names(prop2)=1:length(prop2[,1]);prop2=melt(as.matrix(prop2));
prop$Var1=tis1[as.numeric(prop$Var1)];prop2$Var1=tis1[as.numeric(prop2$Var1)]
prop$sd=as.numeric(prop2$value);prop$value=as.numeric(prop$value)
prop$Var1=factor(prop$Var1,levels=tis1)

ggplot(prop,aes(x=value*100,y=Var2,fill=Var2))+geom_bar(stat="identity")+facet_wrap(~Var1)+
  geom_errorbar(aes(xmin=value*100,x=value*100,xmax=(value+sd)*100),stat="identity")+theme_bw()+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.text=element_text(size=15),strip.text = element_text(size=15))+
  xlab('Cell %')+ylab('')+scale_fill_manual(name="Cell type",values=colorRampPalette(brewer.pal(10,"Paired"))(length(unique(prop$Var2))))

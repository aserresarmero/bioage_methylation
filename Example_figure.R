library(MASS);library(ape);library(nlme);library(phytools)

nn=600
l=5
cccex=c(1.8,2,3)+0.3
### Colorize black png silhouettes 
change_cols = function(replace_black, replace_white, theimg) {
  r_b = col2rgb(replace_black) / 255
  r_w = col2rgb(replace_white) / 255
  theimg[theimg == 1] <- 2
  for (i in 1:3) {
    theimg[,,i][theimg[,,i] < 1] <- r_b[i]
  }
  for (i in 1:3) {
    theimg[,,i][theimg[,,i] > 1] <- r_w[i]
  }
  return(theimg)
}

## Construct barplots signifying lifespan 
bplot=function(b,p,q,ttt){
  #offset from top
  whe=.2
  #Normalize lifespans
  b3=table(b)/max(table(b))*.4
  b=sort(unique(b))
  b2=b/max(b)
  par(mar=c(5.1,5.1,5.1,2.1))
  plot(0:(length(b)+1),0:(length(b)+1),ty="n",xlim=c(-0.25,1.25),axes=F,xlab="Chronological age (years)",ylab="",cex.axis=cccex[1],cex.lab=cccex[3]-0.9,cex.names=cccex[2],)
  mtext(text = ttt, side = 2,line = 0.6,at = 3,cex=cccex[3]-1.3)
  sapply(1:length(b2),function(x){
    im=change_cols(palette()[x+1],"transparent",p[[x]])
    im2=change_cols(palette()[x+1],"transparent",q[[x]])
    rasterImage(im,-0.4,x-0.6,0,x+0.4)
    rasterImage(im2, b2[x],x-0.6,b2[x]+.35,x+0.4)
    rect(0,x-b3[x],b2[x],x+b3[x],col=adjustcolor(x+1,.6),border=x+1,lwd=3)
    rect(0,x-b3[x],whe,x+b3[x],col=adjustcolor(x+1,.7),border = "transparent")
    text(whe/2,x,paste0(round(.2/b2[x]*100,0),"%"),cex=cccex[1])
  })
  axis(1,c(0,whe,b2),round(c(0,whe*max(b),b),0),cex.axis=cccex[1],cex.lab=cccex[2],cex.names=cccex[2])
  abline(v=whe,lty="dashed")
  points(whe,length(b2)+0.6,pch=25,bg="black",cex=cccex[3])
  text(whe,length(b2)+1,paste0(whe*max(b)," years"), cex=cccex[1])
}

## Build and write regression formula in panel
form=function(m,pos){
  mes=pos[1]/(50)
  mesy=pos[3]/(55)
  co=abs(c(log(2*(1-pt(abs(coefficients(summary(m))[-1,3][-4]),995))+1E-200,10)))
  co=sapply(co,function(x)min(x,150))
  pts=seq(pos[1],pos[1]*1.6,length.out=length(co))
  mno=40
  points(pts,rep(pos[2],length(co)),cex=co/mno)
  points(pts[1:(length(pts)-1)]*1.075,rep(pos[2],length(co)-1),pch=3,cex=.8)
  text(pts,rep(pos[2]+mesy*4,length(co)),row.names(coefficients(summary(m)))[-c(1,5)],cex=cccex[1])
}

## Generate confidence intervals for quadratic and linear regressions
cip=function(x,y,z,ra=1){
  by(cbind.data.frame(x,y=(y**2),as.numeric(z)),z,function(ff){
    mod=lm(y~x,data=ff);n=sort(sample(ff[,1],50,replace = T))
    p=predict(mod,data.frame(x=n),interval = "confidence",level=.99)
    lines(n,p[,2]**(ra),lty="dashed",type="l",xlim=c(0,max(x)),
          ylim=c(0,max(y**ra,na.rm = T)),main="",col=ff[1,3]+1,lwd=2)
    lines(n,p[,3]**(ra),lty="dashed",type="l",xlim=c(0,max(x)),
          ylim=c(0,max(y**ra,na.rm = T)),main="",col=ff[1,3]+1,lwd=2)
    polygon(c(rev(n), n), c(rev(p[ ,3]**ra), p[ ,2]**ra), col = adjustcolor(ff[1,3]+1,.45), border = NA)
    mmm=seq(min(x),max(x),length.out=100)
    lines(mmm,predict(mod,data.frame(x=mmm),interval = "confidence",level=.99)[,1]**ra,
          col=ff[1,3]+1, lty="dashed")
  })
}

## Generate the whole plot: Plot regression lines, formulas, CIs and trees
spiel=function(x1,x2,x3,x4,x4c,x5,y,e,ra=1){
  if(ra==1){r=2}else{r=1}
  names(x3)=paste0("t",1:nn)
  names(x4)=paste0("t",1:nn)
  d=as.phylo(hclust(as.dist(as.matrix(dist(x3)*1.3)^0.2)*4))
  ## Determine inner tree edge colors
  ns=(as.numeric(as.factor(round(c(x3,fastAnc(d,x3)),0)))+1)
  d$tip.label=rep(".", length(d$tip.label))
  ## This formula is contaimanted by a tree from interaction variable x3
  mod=lm(y~Age+Lsp+Int+Tree+Sex,
    data=cbind.data.frame(Age=x1,Lsp=x2,Tree=x1/x3,Int=x1/x2,y=y**r,Sex=x5))
  ## This formula is not so contaimanted by a tree bc x4 is permuted
  mod2=lm(y~Age+Lsp+Int+Tree+Sex,
    data=cbind.data.frame(Age=x1,Lsp=x2,Tree=x1/x4,Int=x1/x2,y=y**r,Sex=x5))
  plot(x1,y**r,bg=adjustcolor(as.numeric(factor(x2))+1,.8),col=adjustcolor("black",.5),
    xlab="Chronological age (years)",ylab="Methylation age (arbitrary)",pch=c(21,22)[x5],
    bty="n",ylim=c(min(y**r),max(y**r)),cex.axis=cccex[1],cex.lab=cccex[2],cex.names=cccex[2])
  form(mod,c(max(x1)*.45,(max(y**r)-min(y**r))*0.35+min(y**r),max(y**r)*0.55))
  form(mod2,c(max(x1)*.45,(max(y**r)-min(y**r))*0.15+min(y**r),max(y**r)*0.55))
  legend(1,max(y**2)*1,c("Female","Male"),pch = c(21,22),bty = "n",cex=cccex[2])
  cip(x=x1,y=y,z=factor(x2),ra=ra)
  par(new=T);plot.phylo(d,type="fan",edge.color=ns[match(d$edge[,2],1:length(ns),nomatch = 0)],
    tip.color="transparent",x.lim = c(-29,1)+2,y.lim = c(-6,9)-0.5, edge.width = 0.5)
  ring((x2/max(x2)*30+rnorm(nn,0,4))/100*1.5,offset = 10/100,d,col=palette()[as.numeric(cut(x2,l))+1])
  par(new=T);plot.phylo(d,type="fan",edge.color=ns[match(d$edge[,2],1:length(ns),nomatch = 0)],
    tip.color="transparent",x.lim = c(-29,1)+2,y.lim = c(-1,14)+1.3,edge.width = 0.5)
  ring((as.numeric(x4)/max(x4)*30+rnorm(nn,0,7))/80*1.5,offset = 10/100,
    phy = d,col=palette()[x4c[as.integer(cut(x4,l))]+1])
}


set.seed(11)
library(ape)
library(phytools)
library(png)

pic<-readPNG("~/Documents/doggo/Slide1.png")
chi=readPNG("~/Documents/doggo/chihuahua.png")
whwt<-readPNG("~/Documents/doggo/whwt.png")
bcol=readPNG("~/Documents/doggo/border_collie.png")
gshp=readPNG("~/Documents/doggo/german_shepherd.png")
grd<-readPNG("~/Documents/doggo/greatdane.png")

baby<-readPNG("~/Documents/baby.png")
bones<-readPNG("~/Documents/bones.png")
bones2<-readPNG("~/Documents/doggo/bones2.png")

p=list(grd,gshp,bcol,whwt,chi)
bab=rep(list(baby),5)
bon=rep(list(bones),5)
bon2=rep(list(bones2),5)

pdf("~/test_examp.pdf",width=15, height=9)
#tiff("~/Figures_methylation_paper/Figure1.tiff",width=15*72, height=9*72)
layout(t(sapply(c(0,3),function(x)x+c(rep(1,2),rep(2:3,each=3)))))

## Do spiel for human and dog simulated lifespans: the coefficients are arbitrary but 
## they are the same for human and dog

## Beware!!
## Effectively this is not a generailzed model, but it plots much nicer.
## I'm attaching a similar script below using generalized errors. The p-values are
## similar but the regression lines are a lot messier.
for (x in 1:2){
    if(x==1){
      ##Sample human
      x2=sample(seq(65,100,l=l),nn,replace = T,
      prob=1.3^(1:length(seq(65,100,l=l))*(rev(1:length(seq(65,100,l=l))))))
    }else{
      ## sample dog
      x2=sample(seq(8,17,l=l),nn,replace = T)
    }
    ## make a variable x3 that follows x2 
    x3=5+abs((1:length(unique(x2)))[as.integer(factor(x2))]+rnorm(nn,0,.1))
    ## make a variable x4 that permutes x2 but keeps assignations 
    x4c=sample(1:length(unique(x2)))
    x4=5+x4c[as.integer(factor(x2))]+rnorm(nn,0,.4)
    ## create ages 
    x1=sapply(x2,function(x)runif(1,0,x))
    ## add sex for good measure
    x5=sample(1:2,nn,replace=T)
    ## make an error (if I had scaled the variables, then it would't be so nasty)
    e=abs(rnorm(nn,4,3))
    ## Because the variables are not scaled I need to give everything arbitrary weights
    y=sqrt(4*x1+x1/x2*50+x5*2+3*e)
    par(mar=c(8,2,8,2))
    if(x==1){
      bplot(x2,bab,bon,ttt="HUMAN")
    }else{bplot(x2,p,bon2,ttt="DOG")}
    par(mar=c(4.5,5,4,1))
    spiel(x1,x2,x3,x4,x4c,x5,y,e,1/2)
    spiel(x1,x2,x3,x4,x4c,x5,y,e,1)
  }
mtext(c("A","B","C","D","E","F"),side = 1,line = rep(c(-60,-27),each=3),at = rep(c(-94,-68,-30),2),font=2,cex=cccex[3]-0.9)
dev.off()

## Build model likelihood 
# lik=function(X,y,cova){
#   #like=function(X,C,y,bet){(-0.5*t(y-as.vector(X%*%bet))%*%C%*%(y-as.vector(X%*%bet))-0.5*deter)[1]-n/2*log(2*pi)}
#   C=cova
#   carro=ginv(t(X)%*%C%*%X)
#   bet=carro%*%t(X)%*%C%*%y
#   err=as.numeric(t(y-X%*%bet)%*%C%*%(y-X%*%bet)/(dim(C)[1]-length(bet)-1))
#   ne=as.vector(t(bet)/sqrt(abs(diag(err*carro))))
#   f=cbind.data.frame(beta=bet,err=sqrt(abs(diag(err*carro))),tval=ne,pval=2*pt(abs(ne),dim(C)[1]-length(bet)-1,lower.tail=F))
#   row.names(f)=colnames(X)
#   f
# }

# nn=600
# sapply(1:2,function(x){
#   spiel=function(x1,x2,x3,x4,x5,x5c,cod,y1,y,e){
#     cod$edge.length=cod$edge.length/median(cod$edge.length)/1000
#     d=cod
#     ns=as.numeric(as.factor(round(c(rep(0,length(x2)),fastAnc(d,x2)),0)))+1
#     d$tip.label=rep(".", length(d$tip.label))
#     mod=gls(y~Age+Lsp+Int+Sex,data=cbind.data.frame(Age=x1,Lsp=x2,Int=x1/x2,y=y1**r,Sex=x3,Miau=x4),correlation = corPagel(value=1, phy=co,fixed = T))
#     #mod2=gls(y~Age+Lsp+Int+Sex,data=cbind.data.frame(Age=x1,Lsp=x5,Int=x1/x5,y=y**r,Sex=x3),correlation = corPagel(value=1, phy=co,fixed = T))
#     mod2=gls(y~Age+Lsp+Int+Sex,data=cbind.data.frame(Age=x1,Lsp=x5,Int=x1/x5,y=y**r,Sex=x3,Miau=x4),correlation = corPagel(value=1, phy=co,fixed = T))
#     plot(x1,ty**r,bg=adjustcolor(as.numeric(factor(x2))+1,.8),col=adjustcolor("black",.5),xlab="Chronological age (years)",ylab="Methylation age (arbitrary)",pch=c(21,22)[x3], bty="n",ylim=c(min(ty**r),max(ty**r)),cex.axis=1.3,cex.lab=1.5)
#     form(m=mod,pos=c(max(x1)*.45,(max(ty**r)-min(ty**r))*0.35+min(ty**r),max(ty**r)*0.55))
#     form(m=mod2,pos=c(max(x1)*.45,(max(ty**r)-min(ty**r))*0.15+min(ty**r),max(ty**r)*0.55))
#     legend(1,max(ty**2)*1,c("Female","Male"),pch = c(21,22),bty = "n")
#     cip(x1,ty,factor(x2))
#     par(new=T);plot.phylo(d,type="fan",edge.color=ns[match(d$edge[,2],1:length(ns),nomatch = 0)],tip.color="transparent",x.lim = c(-30,3)+1,y.lim = c(-7,10)-0.5)
#     ring((x4/max(x4)*30+rnorm(nn,0,4))/100,offset = 1/100,d,col=palette()[as.numeric(as.factor(x2))+1])
#     par(new=T);plot.phylo(d,type="fan",edge.color=ns[match(d$edge[,2],1:length(ns),nomatch = 0)],tip.color="transparent",x.lim = c(-30,3)+1,y.lim = c(-2,15)+0.1)
#     ring(x5/20,offset = 1/100,d,col=palette()[x5c[as.integer(factor(x2))]+1])
#   }
#   set.seed(2)
#   l=6
#   if(x==1){
#     x2=sample(seq(65,100,l=l),nn,replace = T,prob=1.5^(1:length(seq(65,100,l=l))*(rev(1:length(seq(65,100,l=l))))))
#   }else{x2=sample(seq(8,17,l=l),nn,replace = T)}
#   x2t=x2
#   x22=as.factor(x2t)
#   i=T
#   while(i==T){
#     levels(x22)=sample(levels(x22))
#     if(sum(levels(x22)==levels(factor(x2t)))<length(levels(x22))/3){i=F}
#   }
#   x22=as.numeric(levels(x22)[as.numeric(x22)])
#   names(x2)=paste0("t",1:nn)
#   names(x22)=paste0("t",1:nn)
#   x22t=x22
#   if(x==1){x2=x2+rnorm(nn,0,2);x22=x22+rnorm(nn,0,6)}else{x2=x2+rnorm(nn,0,2)/5;x22=x22+rnorm(nn,0,6)/5}
#   co=sapply(norm2(x2),function(x)sapply(norm2(x2),function(y)abs(x-y)))
#   co=as.matrix(nearPD(cov.from.dist(co))$mat)
#   x1=sapply(x2t,function(y)abs(runif(1,0,y)))
#   x3=sample(1:2,nn,replace=T)
#   x11=sapply(x22t,function(y)abs(runif(1,0,y)))
#   err=co%*%rnorm(nn)
#   cod=as.phylo(hclust(as.dist(dist.from.cov(co))))
#   y1=2*norm2(x1)+10*norm2(norm2(x1)*norm2(x2))+0.3*norm2(x3)+3*err
#   y2=2*norm2(x11)+10*norm2(norm2(x11)*norm2(x22))+0.3*norm2(x3)+3*err
#   X=cbind.data.frame(x1=norm2(x1),x2=norm2(x2),x3=norm2(x3),x11=norm2(x11),x22=norm2(x22),y1,y2)
#   plot(x1,y1,bg=adjustcolor(as.numeric(factor(x2t))+1,.8),col=adjustcolor("black",.5),xlab="Chronological age (years)",ylab="Methylation age (arbitrary)",pch=c(21,22)[x3], bty="n",cex.axis=1.3,cex.lab=1.5)
#   cip(x1,y1,factor(x2t))
#   plot(x11,y2,bg=adjustcolor(as.numeric(factor(x22t))+1,.8),col=adjustcolor("black",.5),xlab="Chronological age (years)",ylab="Methylation age (arbitrary)",pch=c(21,22)[x3], bty="n",cex.axis=1.3,cex.lab=1.5)
#   cip(x11,y2,factor(x22t))
#   
#   m1=lik(as.matrix(cbind.data.frame(cept=rep(1,dim(X)[1]),X[,1:3],Int=norm2(X[,1]*X[,2]))),y=as.vector(X[,6]),ginv(co))
#   m2=lik(as.matrix(cbind.data.frame(cept=rep(1,dim(X)[1]),x1=X[,4],x2=X[,5],x3=X[,3],Int=norm2(X[,4]*X[,5]))),y=as.vector(X[,7]),ginv(co))
#   summary(gls(model = y1~x1+x2+x3+x1:x2,data=X,correlation = corBrownian(1,cod)))$tTable
#   summary(gls(model = y2~x11+x22+x3+x11:x22,data=X,correlation = corBrownian(1,cod)))$tTable
#   
#   x11(height=10,width=10)
#   par(mar=c(0,0,0,0))
#   contMap( cod, norm2(x2),ftype="off",col="transparent",type="fan",res = 20, legend=F,xlim=c(-10,10),ylim=c(-10,10))
#   q=contMap( cod, norm2(x2),ftype="off",col="transparent",type="fan",res = 20, legend=F,plot=F)
#   ring(norm2(x2)/9.95,cod,col="black",offset=0)
#   ring(norm2(x2)/10,cod,col=q$cols[floor(seq(1,length(q$cols),l=10))][cut(x2,10)],offset=0)
#   contMap( cod, norm2(x2),col="transparent",ftype="off",type="fan",res = 20, legend=F)
#   ring(norm2(x22)/9.95,cod,col="black",offset=0)
#   ring(norm2(x22)/10,cod,col=q$cols[floor(seq(1,length(q$cols),l=20))][as.numeric(cut(norm2(x22)*5,20))],offset=0)
#   #contMap( cod, x22,tip.color="transparent", offset=0)
#   #plot(as.phylo(cod),type="fan",edge.color=ns[match(d$edge[,2],1:length(ns),nomatch = 0)])
#   par(mar=c(8,2,8,2))
#   # if(x==1){
#   #   bplot(x2,bab,bon)
#   # }else{bplot(x2,p,bon2)}
#   par(mar=c(4,4,4,2))
#   spiel(x1,x2,x3,x4,x5,x5c,co,y1,y,e,1/2)
#   spiel(x1,x2,x3,x4,x5,x5c,co,y1,y,e,1)
# })
# dev.off()

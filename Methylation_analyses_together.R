##### Dependencies and functions:
set.seed(7)
cccex=c(1.5,1.8,2.5)
library(msir)
library(MASS)
library(ape)
library(stringdist)
library(RColorBrewer)
library(glmnet)
library(raster)
library(reshape2)
library(grDevices)
library(gridExtra)
library(nlme)
library(parallel)
library(rwc)
library(Matrix)
library(kinship2)
library(phytools)


norm2=function(bla){(bla - min(bla)) / ( max(bla) - min(bla) )+0.0001}


##### GWAS

#### Generalized regression implementation featuring inverse instead of GLS package 
### Takes in methylation matrix X, vcv matrix C and response variable y 
genreg_inv=function(X,C,y){
  carro=ginv(t(X)%*%C%*%X)
  bet=carro%*%t(X)%*%C%*%y
  err=as.numeric(t(y-X%*%bet)%*%C%*%(y-X%*%bet)/(dim(C)[1]-length(bet)-1))
  ne=as.vector(t(bet)/sqrt(abs(diag(err*carro))))
  cbind.data.frame(beta=bet,tvl=ne,pval=2*pt(abs(ne),dim(C)[1]-length(bet)-1,lower.tail=F))
}

#### Wrapper call to GWAS over datasets
# Takes in methylation matrix met to be treated as response (will be linearized) to loop over,
# vcv matrix cova, covariate matrix data (including age and moderator) and quot (whether to treat the moderator as quotient T)
gwas<-function(met,cova,data,quot=T){
  C=ginv(cova)
  deter=determinant(2*pi*cova,log=T)[[1]][1]
  ### Jitter matrix in case of degenerateness 
  if(any(is.infinite(deter))){
    R=diag(0.005,dim(C)[1])
    R[upper.tri(R)]=rnorm(dim(R)[1]*(dim(R)[1]-1)/2,0,0.01)
    R=R+t(R)
    C=solve(as.matrix(nearPD((cova)**(1+R))[[1]]))
    deter=determinant(2*pi*nearPD((cova)**(1+R))[[1]],log=T)[[1]][1]
  }
  ##Likelihood of BLUP coefficients calculation from https://web.vu.lt/mif/a.buteikis/wp-content/uploads/PE_Book/4-6-Multiple-GLS.html
  like=function(X,C,y,bet,err){(-0.5/err*t(y-as.vector(X%*%bet))%*%C%*%(y-as.vector(X%*%bet))-0.5*deter)[1]-dim(cova)[1]*log(err)}
  apply(met, 1,function(x){
    la=tryCatch({bc=boxcox(x~ag+.,data=data,plotit = F);bc$"x"[which.max(bc$y)]},error=function(e)1)
    if(la>0.05 | la<(-0.05)){y=(x^la-1)/la}else{y=log(x)}
    #X=as.matrix(cbind.data.frame(Intercept=rep(1,dim(data)[1]),met=x,data[,-which(colnames(data)=="ag")],row.names = NULL))
    X=as.matrix(cbind.data.frame(Intercept=rep(1,dim(data)[1]),data,row.names = NULL))
    m1=genreg_inv(X,C,y)
    ## Error variance estimator from https://web.vu.lt/mif/a.buteikis/wp-content/uploads/PE_Book/4-6-Multiple-GLS.html
    err=as.numeric(t(y-X%*%m1[,1])%*%C%*%(y-X%*%m1[,1])/(dim(C)[1]-length(m1[,1])-1))
    l1=like(X,C,y,m1[,1],err)
    row.names(m1)=c("Intercept",colnames(data))
    if(quot==T){
      la=tryCatch({bc=boxcox(x~ag+.+I(ag/lsp),data=data,plotit = F);bc$"x"[which.max(bc$y)]},error=function(e)1)
      if(la>0.05 | la<(-0.05)){y=(x^la-1)/la}else{y=log(x)}
      X=as.matrix(cbind.data.frame(Intercept=rep(1,dim(data)[1]),data,Int=data$ag/data$lsp))
    }else{
      la=tryCatch({bc=boxcox(x~ag+.+I(ag*lsp),data=data,plotit = F);bc$"x"[which.max(bc$y)]},error=function(e)1)
      if(la>0.05 | la<(-0.05)){y=(x^la-1)/la}else{y=log(x)}
      X=as.matrix(cbind.data.frame(Intercept=rep(1,dim(data)[1]),data,Int=data$ag*data$lsp))
    }
    m2=genreg_inv(X,C,y)
    err=as.numeric(t(y-X%*%m2[,1])%*%C%*%(y-X%*%m2[,1])/(dim(C)[1]-length(m2[,1])-1))
    l2=like(X,C,y,m2[,1],err)
    row.names(m2)=c("Intercept",colnames(data),"Int")
    X=as.matrix(cbind.data.frame(Intercept=rep(1,dim(data)[1]),rev(data)))
    m3=genreg_inv(X,C,x)
    row.names(m3)=c("Intercept",colnames(rev(data)))
    #list(m1,m2,m3,pchisq(2*abs(l2-l1),1,lower.tail = F))
    ## Report main effect p-values, change in likelihood from adding interaction p-value,
    ## moderator p-value, and interaction effect p-value respectively
    res=cbind(m1[2,dim(m1)[2]],ifelse(any(is.infinite(c(l1,l2))),Inf,
      pchisq(2*abs(l2-l1),1,lower.tail = F)),m3[2,dim(m3)[2]],m2[dim(m2)[1],dim(m2)[2]])
  })
}

#### Plot GWAS dual density p-values
## Uses values from GWAS function. Takes in naked moderator p-values (plif_no_int), 
## moderator interaction p-values(pag_int), and main effect p-values (pag_no_int)
## bio and chron plot the active coefficients (sites) of bioage and chronage glmnet models
pplot=function(plif_no_int,pag_int,pag_no_int,bi=1200,bio,chron){
  #ss="lambda.min"
  ## Max number of coefficients in regression
  mno=20
  ss=chron$lambda[chron$nzero>mno][1];sss=bio$lambda[bio$nzero>mno][1];
  ## Normalize p-values for axes
  nplif=norm2(-log(plif_no_int,10))
  npai=norm2(-log(pag_int,10))
  npan=norm2(-log(pag_no_int,10))
  ## Number of bins to report in legend
  lll=4
  ## Number of samples from gradient
  k=30
  li=c(-0.1,1.1)
  for (i in 1:2){
    par(xpd=T)
    if(i==1){ppp=npan;par(mar=c(0,5,2,4));nam="age";tmp=pag_no_int}
      else{ppp=nplif;par(mar=c(2,5,0,4));nam="lifespan"; tmp=plif_no_int}
    #make 2d density and single axis densities
    kd=kde2d(npai,ppp,n = bi,lims=c(li,li))
    d1=density(npai)
    d2=density(ppp)
    kdn=as.data.frame(raster(kd),xy=T)
    ## make colors less sharp 
    kdn$layer=(kdn$layer+min(kdn$layer[kdn$layer>0]))^(1/5)
    #Plot 2d density 
    #par(mar=c(2.1,7.1,2.1,2.1))
    plot(kdn[,1:2],bg=adjustcolor(c("transparent",colorRampPalette(RColorBrewer::brewer.pal(8,"Paired"))(k*2)[1:(k-1)])[cut(kdn$layer,k)],1),
         cex=1.1,pch=22,col="transparent", xlim=li,ylim=li,bty="n", axes=F,xlab ="",
         ylab=bquote(paste('-log'['10']*' (', .(nam),')')))
    par(new=T)
    #Plot 2d density contours
    contour(kd,frame.plot = F, xlim=li,ylim=li,xlab="",ylab="", axes=F)
    yintr=abs(log10(median(tmp)))/abs(log10(0.5))*log10(length(tmp)/0.01)
    xintr=abs(log10(median(npai)))/abs(log10(0.5))*log10(length(tmp)/0.01)
    par(xpd=F)
    ### Make ""significance"" line, kind of meaningless
    liny=function(x){abs(yintr/log10(min(tmp)))+-abs(yintr/log10(min(tmp))/xintr*log10(min(pag_int)))*x}
    segments(x0 = 0.001,y0 = liny(0.001), x1 = 1,y1=0,lty="dashed")
    #abline(v=0.001,h=0.001, lty="dashed",lwd=2)
    par(xpd=T)
    ## Plot active coefficients bioage
    points(npai[abs(coefficients(bio,s=sss)[-1])>0],ppp[abs(coefficients(bio,s=sss)[-1])>0],
           col="black",bg="gold",pch=21,cex=1.7,xlim=li,ylim=li)
    ## Plot active coefficients chronage
    points(npai[abs(coefficients(chron,s=ss)[-1])>0],ppp[abs(coefficients(chron,s=ss)[-1])>0],
           col="black",bg="lightgray",pch=21,cex=1.7,xlim=li,ylim=li)
    ## Plot margin densities
    if(i==1){
      legend(li[2]/2,li[2]+0.25,legend = c("Chronological","Biological"),fill=c("lightgray","gold"),title = "Clock CpGs", bg="transparent", cex=1.8,bty='n')
      polygon(d1$x,d1$y/sum(d1$y)*1.6+li[2]*1.015,col ="black")
    }
    polygon(d2$y/sum(d2$y)*1.2+li[2]*1.015,d2$x,col ="black")
    axis(2,at=seq(li[1],li[2],length.out=8),labels = -round(seq(0,-max(-log(tmp,10)),length.out=8),1))
  }  
  axis(1,at=seq(li[1],li[2],length.out=8),labels = -round(seq(0,-max(-log(pag_int,10)),length.out=8),1))
  legend_image <- as.raster(matrix(colorRampPalette(RColorBrewer::brewer.pal(8/2,"Paired"))(k*2), ncol=1))
  text(x=li[2]/1.8, y = 1, labels =expression("-log"[10]~"(count)"),cex=1.8)
  text(x=li[2]/1.7, y = seq(li[2]/2,li[2]/2*1.5,l=lll), labels = round(rev(seq(0,log(length(plif_no_int),10),l=lll)),0),cex=1.5)
  rasterImage(legend_image, li[2]/2, li[2]/2, li[2]/2*1.1,li[2]/2*1.5)
}

#### Plot Global regression
## Takes in all clock objects belonging to a dataset.
## met is the methylation matrix
## Additional option sparse=T limits number of coefficients to less than 25
## Additional option pgls=T plots GLS predictions instead of OLS predictions
## age is a vector of ages belonging to the samples and vari are the moderators
plreg<-function(met,chron,chron2,bio,bio2,pace,pace2,age,vari,vcvv=1,scAge,scAge2,sparse=T,pgls=T,quot=F){  ## Make confidence intervals function
  cip=function(x,y,z,ra=1){
    by(cbind.data.frame(x,y=y,as.numeric(z)),z,function(ff){
      mod=lm(y~x,data=ff);n=sort(sample(ff[,1],50,replace = T))
      p=predict(mod,data.frame(x=n),interval = "confidence",level=.99)
      lines(n,p[,2]**(ra),lty="dashed",type="l",xlim=c(0,max(x)),ylim=c(0,max(y**ra,na.rm = T)),main="",col=ff[1,3]+1,lwd=2)
      lines(n,p[,3]**(ra),lty="dashed",type="l",xlim=c(0,max(x)),ylim=c(0,max(y**ra,na.rm = T)),main="",col=ff[1,3]+1,lwd=2)
      polygon(c(rev(n), n), c(rev(p[ ,3]**ra), p[ ,2]**ra), col = adjustcolor(ff[1,3]+1,.45), border = NA)
      mmmm=seq(min(x),max(x),length.out=100)
      lines(mmmm,predict(mod,data.frame(x=mmmm),interval = "confidence",level=.99)[,1]**ra,col=ff[1,3]+1, lty="dashed")
    })
  }
  if(is.null(dim(vcvv))){vcvv=diag(rep(1,dim(met)[2]))}
  # Max Number of coefficients allowed
  mno=20
  ## adapt glmnet parameter search lambda parameter "s" for each regression object.
  ## Function got a little ugly as more features were added
  if(sparse==T){
    s=chron$lambda[chron$nzero>mno][1];ss=chron2$lambda[chron2$nzero>mno][1];
    sss=bio$lambda[bio$nzero>mno][1];ssss=bio2$lambda[bio2$nzero>mno][1];
  }else{
    s="lambda.min";ss="lambda.min";sss="lambda.min";ssss="lambda.min"
  }
  tab=data.frame(matrix(NA, nrow = 4, ncol = 4),row.names = c("Mean","Age","Mod","Int"))
  tab2=data.frame(matrix(NA, nrow = 4, ncol = 4),row.names = c("Mean","Age","Mod","Int"))
  qqq=c()
  qqq2=c()
  ### Create predictions for bioage and chronage, hanling the different number of coefficients
  ### Loop over different clocks
  for (i in 1:4){
    if(i==1){
      if(sparse==T){
          ppp=as.vector(predict(chron,newx = t(met),s=s));pppp=as.vector(predict(chron2,newx = t(met),s=ss));
          nam=paste0("CC (N = ",chron$nzero[chron$nzero>mno][1],")")
          }else{ppp=as.vector(predict(chron,newx = t(met),s=s));
          pppp=as.vector(predict(chron2,newx = t(met),s=ss));nam=paste0("CC (N = ",chron$nzero[chron$index[1]],")")
        }}
    else if(i==2){
      if(sparse==T){
        ppp=as.vector(predict(bio,newx = t(met),s=sss));nam=paste0("BC (N = ",bio$nzero[bio$nzero>mno][1],")");pppp=as.vector(predict(bio2,newx = t(met),s=ssss));
        }else{ppp=as.vector(predict(bio,newx = t(met),s=sss));nam=paste0("BC (N = ",bio$nzero[bio$index[1]],")");pppp=as.vector(predict(bio2,newx = t(met),s=ssss));
      }}
    else if(i==3){ppp=as.vector(pace);pppp=as.vector(pace2);nam="PM"}else if(i==4){ppp=as.vector(scAge);pppp=as.vector(scAge2);nam="SC"};
    ## Linearize coefficients with box-cox
    la1=tryCatch({bc=boxcox(I(ppp)~age,plotit = F);bc$"x"[which.max(bc$y)]},error=function(e)1)
    la2=tryCatch({bc=boxcox(I(pppp)~age,plotit = F);bc$"x"[which.max(bc$y)]},error=function(e)1)
    if(la1>0.2 | la1<(-0.2)){resp=(ppp^la1-1)/la1}else{resp=log(ppp)}
    if(la2>0.2 | la2<(-0.2)){resp2=(pppp^la2-1)/la1}else{resp2=log(pppp)}
    ## Plot clocks. If scage or pacemaker represent OLS instead of GLS
    if(pgls==T | i>2){
      plot(age,resp,bg=as.numeric(cut(vari,4))+1,pch=21,cex=1.9,main=nam,xlab="",ylab="",cex.axis=2.4)
      cip(age,resp,factor(as.numeric(cut(vari,4))+1),ra = 1)
    }else{
      plot(age,resp2,bg=as.numeric(cut(vari,4))+1,pch=21,cex=1.9,main=nam,xlab="",ylab="",cex.axis=2.4)
      cip(age,resp2,factor(as.numeric(cut(vari,4))+1),ra = 1)
    }
    ### Extract model p-values for representation
    mm=lm(Met~Age+Mod+Int,cbind.data.frame(Met=resp,Age=age,Mod=vari,Int=((if(quot==F){vari}else{1/vari})*age)))
    mmm=lm(Met~Age+Mod+Int,cbind.data.frame(Met=resp2,Age=age,Mod=vari,Int=((if(quot==F){vari}else{1/vari})*age)))
    qqq=c(qqq,2*summary(mm)$"r.squared")
    if(i>2){
      tab[,i]=as.vector(2*pt(-abs(summary(mm)$coefficients[,3]),df=dim(met)[2]))
      tab2[,i]=as.vector(2*pt(-abs(summary(mmm)$coefficients[,3]),df=dim(met)[2]))
    }else{
      tab2[,i]=as.vector(2*pt(-abs(summary(mm)$coefficients[,3]),df=dim(met)[2]))
      tab[,i]=as.vector(2*pt(-abs(summary(mmm)$coefficients[,3]),df=dim(met)[2]))
    }
  }
  ## Represent p-values as circles
  tab=melt(as.matrix(tab[length(tab[,1]):1,]))
  tab2=melt(as.matrix(tab2[length(tab2[,1]):1,]))
  par(mar=c(4,4,3,2))
  plot(as.numeric(as.factor(tab[,2]))-rep(c(0.15,0),each=length(as.numeric(as.factor(tab2[,2])))/2),
       as.numeric(as.factor(tab[,1])),cex=1.6*sqrt(-log(tab[,3]+1E-300,10)/70+0.2+ifelse(-log(tab[,3],10)>5,0.3,0))*2,
       axes=F,xlab="",ylab="",pch=20,col=ifelse(tab[,3]<0.01,"red","gray"),ylim=c(0,length(unique(tab[,1]))),
       xlim=c(0,max(as.numeric(as.factor(tab[,2])))+0.3))
  points((as.numeric(as.factor(tab2[,2]))+0.15)[1:(length(as.numeric(as.factor(tab2[,2])))/2)],
         as.numeric(as.factor(tab2[,1]))[1:(length(as.numeric(as.factor(tab2[,2])))/2)],
         cex=1.6*sqrt(-log(tab2[,3]+1E-300,10)/70+0.2+ifelse(-log(tab2[,3],10)>5,0.3,0))*2,pch=20,
         col=ifelse(tab2[,3]<0.01,"blue","gray"),ylim=c(0,length(unique(tab2[,1]))))
  points(1:length(unique(tab[,2])),rep(0,length(unique(tab[,2]))),cex=2*1.8,pch=20,col="black")
  points(1:length(unique(tab[,2])),rep(0,length(unique(tab[,2]))),cex=qqq*1.8,pch=20,col="orange")
  axis(3,at = 1:length(unique(tab[,2])),labels = c("CC","BC","PM","SC"),lwd = 0,cex=1.6,gap.axis=-10000)
  axis(2,at = 0:length(unique(tab[,1])),labels = c(expression("R"^2~" of 1"),
        as.character(unique(tab[,1]))),las=1,lwd = 0,cex=1.6,hadj = 0)
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

####Test
plreg2<-function(met,chron,chron2,bio,bio2,abio,abio2,pace,pace2,apace,apace2,age,vari,vari2,vcvv=1,scAge,scAge2,ascAge,ascAge2,sparse=T,pgls=T){  
  ## Make confidence intervals function
  cip=function(x,y,z,ra=1){
    by(cbind.data.frame(x,y=y,as.numeric(z)),z,function(ff){
      mod=lm(y~x,data=ff);n=sort(sample(ff[,1],50,replace = T))
      p=predict(mod,data.frame(x=n),interval = "confidence",level=.99)
      lines(n,p[,2]**(ra),lty="dashed",type="l",xlim=c(0,max(x)),ylim=c(0,max(y**ra,na.rm = T)),main="",col=ff[1,3]+1,lwd=2)
      lines(n,p[,3]**(ra),lty="dashed",type="l",xlim=c(0,max(x)),ylim=c(0,max(y**ra,na.rm = T)),main="",col=ff[1,3]+1,lwd=2)
      polygon(c(rev(n), n), c(rev(p[ ,3]**ra), p[ ,2]**ra), col = adjustcolor(ff[1,3]+1,.45), border = NA)
      mmmm=seq(min(x),max(x),length.out=100)
      lines(mmmm,predict(mod,data.frame(x=mmmm),interval = "confidence",level=.99)[,1]**ra,col=ff[1,3]+1, lty="dashed")
    })
  }
  if(is.null(dim(vcvv))){vcvv=diag(rep(1,dim(met)[2]))}
  # Max Number of coefficients allowed
  mno=20
  ## adapt glmnet parameter search lambda parameter "s" for each regression object.
  ## Function got a little ugly as more features were added
  if(sparse==T){
    s=chron$lambda[chron$nzero>mno][1];ss=chron2$lambda[chron2$nzero>mno][1];
    sss=bio$lambda[bio$nzero>mno][1];ssss=bio2$lambda[bio2$nzero>mno][1];
    asss=abio$lambda[abio$nzero>mno][1];assss=abio2$lambda[abio2$nzero>mno][1]
  }else{
    s="lambda.min";ss="lambda.min";sss="lambda.min";ssss="lambda.min";asss="lambda.min";assss="lambda.min"
  }
  tab=data.frame(matrix(NA, nrow = 4, ncol = 4),row.names = c("Mean","Age","Mod","Int"))
  tab2=data.frame(matrix(NA, nrow = 4, ncol = 4),row.names = c("Mean","Age","Mod","Int"))
  tab3=data.frame(matrix(NA, nrow = 4, ncol = 4),row.names = c("Mean","Age","Mod","Int"))
  tab4=data.frame(matrix(NA, nrow = 4, ncol = 4),row.names = c("Mean","Age","Mod","Int"))
  qqq=c()
  qqq2=c()
  ### Create predictions for bioage and chronage, hanling the different number of coefficients 
  ### Loop over different clocks
  for (i in 1:4){
    if(i==1){
      if(sparse==T){
        ppp=as.vector(predict(chron,newx = t(met),s=s));pppp=as.vector(predict(chron2,newx = t(met),s=ss));
        nam=paste0("CC (N = ",chron$nzero[chron$nzero>mno][1],")")
      }else{ppp=as.vector(predict(chron,newx = t(met),s=s));
      pppp=as.vector(predict(chron2,newx = t(met),s=ss));nam=paste0("CC (N = ",chron$nzero[chron$index[1]],")")
      }}
    else if(i==2){
      if(sparse==T){
        ppp=as.vector(predict(bio,newx = t(met),s=sss));nam=paste0("BC (N = ",bio$nzero[bio$nzero>mno][1],"), Weight");pppp=as.vector(predict(bio2,newx = t(met),s=ssss));
        appp=as.vector(predict(abio,newx = t(met),s=asss));anam=paste0("BC (N = ",abio$nzero[abio$nzero>mno][1],"), Lifespan");apppp=as.vector(predict(abio2,newx = t(met),s=assss));
      }else{
        ppp=as.vector(predict(bio,newx = t(met),s=sss));nam=paste0("BC (N = ",bio$nzero[bio$index[1]],"), Weight");pppp=as.vector(predict(bio2,newx = t(met),s=ssss));
        appp=as.vector(predict(abio,newx = t(met),s=asss));anam=paste0("BC (N = ",abio$nzero[abio$index[1]],"), Lifespan");apppp=as.vector(predict(abio2,newx = t(met),s=assss));
      }}
    else if(i==3){ppp=as.vector(pace);pppp=as.vector(pace2);nam="PM, Weight";appp=as.vector(apace);apppp=as.vector(apace2);anam="PM, Lifespan"
    }else if(i==4){ppp=as.vector(scAge);pppp=as.vector(scAge2);nam="SC, Weight";appp=as.vector(ascAge);apppp=as.vector(ascAge2);anam="SC, Lifespan"};
    ## Linearize coefficients with box-cox
    la1=tryCatch({bc=boxcox(I(ppp)~age,plotit = F);bc$"x"[which.max(bc$y)]},error=function(e)1)
    la2=tryCatch({bc=boxcox(I(pppp)~age,plotit = F);bc$"x"[which.max(bc$y)]},error=function(e)1)
    if(la1>0.2 | la1<(-0.2)){resp=(ppp^la1-1)/la1}else{resp=log(ppp)}
    if(la2>0.2 | la2<(-0.2)){resp2=(pppp^la2-1)/la1}else{resp2=log(pppp)}
    if(i>1){ala1=tryCatch({bc=boxcox(I(appp)~age,plotit = F);bc$"x"[which.max(bc$y)]},error=function(e)1)
      ala2=tryCatch({bc=boxcox(I(apppp)~age,plotit = F);bc$"x"[which.max(bc$y)]},error=function(e)1)
      if(ala1>0.2 | ala1<(-0.2)){aresp=(appp^ala1-1)/ala1}else{aresp=log(appp)}
      if(ala2>0.2 | ala2<(-0.2)){aresp2=(apppp^ala2-1)/ala1}else{aresp2=log(apppp)}
    }
    ## Plot clocks. If scage or pacemaker represent OLS instead of GLS
    if(pgls==T | i>2){
      plot(age,resp,bg=as.numeric(cut(vari,4))+1,pch=21,cex=1.9,main=nam,xlab="",ylab="",cex.axis=2.4)
      cip(age,resp,factor(as.numeric(cut(vari,4))+1),ra = 1)
      if(i>1){plot(age,aresp,bg=as.numeric(cut(vari2,4))+1,pch=22,cex=1.9,main=anam,xlab="",ylab="",cex.axis=2.4)
      cip(age,aresp,factor(as.numeric(cut(vari2,4))+1),ra = 1)}
    }else{
      plot(age,resp2,bg=as.numeric(cut(vari,4))+1,pch=21,cex=1.9,main=nam,xlab="",ylab="",cex.axis=2.4)
      cip(age,resp2,factor(as.numeric(cut(vari,4))+1),ra = 1)
      if(i>1){plot(age,aresp2,bg=as.numeric(cut(vari2,4))+1,pch=22,cex=1.9,main=anam,xlab="",ylab="",cex.axis=2.4)
        cip(age,aresp2,factor(as.numeric(cut(vari2,4))+1),ra = 1)}
    }
    ### Extract model p-values for representation
    mm=lm(Met~Age+Mod+Int,cbind.data.frame(Met=resp,Age=age,Mod=vari,Int=(vari*age)))
    mmm=lm(Met~Age+Mod+Int,cbind.data.frame(Met=resp2,Age=age,Mod=vari,Int=(vari*age)))
    amm=lm(Met~Age+Mod+Int,cbind.data.frame(Met=resp,Age=age,Mod=vari,Int=((1/vari2)*age)))
    ammm=lm(Met~Age+Mod+Int,cbind.data.frame(Met=resp2,Age=age,Mod=vari,Int=((1/vari2)*age)))
    qqq=c(qqq,2*summary(mm)$"r.squared")
    qqq2=c(qqq2,2*summary(amm)$"r.squared")
    if(i>2){
      tab[,i]=as.vector(2*pt(-abs(summary(mm)$coefficients[,3]),df=dim(met)[2]))
      tab2[,i]=as.vector(2*pt(-abs(summary(mmm)$coefficients[,3]),df=dim(met)[2]))
      tab3[,i]=as.vector(2*pt(-abs(summary(amm)$coefficients[,3]),df=dim(met)[2]))
      tab4[,i]=as.vector(2*pt(-abs(summary(ammm)$coefficients[,3]),df=dim(met)[2]))
    }else{
      tab4[,i]=as.vector(2*pt(-abs(summary(amm)$coefficients[,3]),df=dim(met)[2]))
      tab3[,i]=as.vector(2*pt(-abs(summary(ammm)$coefficients[,3]),df=dim(met)[2]))
      tab2[,i]=as.vector(2*pt(-abs(summary(mm)$coefficients[,3]),df=dim(met)[2]))
      tab[,i]=as.vector(2*pt(-abs(summary(mmm)$coefficients[,3]),df=dim(met)[2]))
    }
  }
  ## Represent p-values as circles
  tab=melt(as.matrix(tab[length(tab[,1]):1,]))
  tab2=melt(as.matrix(tab2[length(tab2[,1]):1,]))
  tab3=melt(as.matrix(tab3[length(tab3[,1]):1,]))
  tab4=melt(as.matrix(tab4[length(tab4[,1]):1,]))
  par(mar=c(4,4,3,2),xpd=NA)
  plot(as.numeric(as.factor(tab[,2]))-rep(c(0.15,0),each=length(as.numeric(as.factor(tab[,2])))/2),
       as.numeric(as.factor(tab[,1]))+0.15,cex=1.6*sqrt(-log(tab[,3]+1E-300,10)/70+0.2+ifelse(-log(tab[,3],10)>5,0.3,0))*2,
       axes=F,xlab="",ylab="",pch=20,col=ifelse(tab[,3]<0.01,"red","gray"),ylim=c(0,length(unique(tab[,1]))),
       xlim=c(0,max(as.numeric(as.factor(tab[,2])))+0.3))
  points((as.numeric(as.factor(tab2[,2]))+0.15)[1:(length(as.numeric(as.factor(tab2[,2])))/2)],
         as.numeric(as.factor(tab2[,1]))[1:(length(as.numeric(as.factor(tab2[,2])))/2)]+0.15,
         cex=1.6*sqrt(-log(tab2[,3]+1E-300,10)/70+0.2+ifelse(-log(tab2[,3],10)>5,0.3,0))*2,pch=20,
         col=ifelse(tab2[,3]<0.01,"blue","gray"),ylim=c(0,length(unique(tab2[,1]))))
  points(as.numeric(as.factor(tab3[,2]))-rep(c(0.15,0),each=length(as.numeric(as.factor(tab3[,2])))/2),
         as.numeric(as.factor(tab3[,1]))-0.15,cex=1.4*sqrt(-log(tab3[,3]+1E-300,10)/70+0.2+ifelse(-log(tab3[,3],10)>5,0.3,0))*2,
         pch=15,col=ifelse(tab3[,3]<0.01,"red","gray"),ylim=c(0,length(unique(tab3[,1]))))
  points((as.numeric(as.factor(tab4[,2]))+0.15)[1:(length(as.numeric(as.factor(tab4[,2])))/2)],
         as.numeric(as.factor(tab4[,1]))[1:(length(as.numeric(as.factor(tab4[,2])))/2)]-0.15,
         cex=1.4*sqrt(-log(tab4[,3]+1E-300,10)/70+0.2+ifelse(-log(tab4[,3],10)>5,0.3,0))*2,pch=15,
         col=ifelse(tab4[,3]<0.01,"blue","gray"),ylim=c(0,length(unique(tab4[,1]))))
  points(1:length(unique(tab[,2]))-0.2,rep(0,length(unique(tab[,2]))),cex=2*1.8,pch=20,col="black")
  points(1:length(unique(tab[,2]))-0.2,rep(0,length(unique(tab[,2]))),cex=qqq*1.8,pch=20,col="orange")
  points(1:length(unique(tab[,2]))+0.2,rep(0,length(unique(tab[,2]))),cex=2*1.6,pch=15,col="black")
  points(1:length(unique(tab[,2]))+0.2,rep(0,length(unique(tab[,2]))),cex=qqq2*1.6,pch=15,col="orange")
  axis(3,at = 1:length(unique(tab[,2])),labels = c("CC","BC","PM","SC"),lwd = 0,cex=1.6,gap.axis=-10000)
  axis(2,at = 0:length(unique(tab[,1])),labels = c(expression("R"^2~" of 1"),
        as.character(unique(tab[,1]))),las=1,lwd = 0,cex=1.6,hadj = 0)
  par(mar=c(5.1, 4.1, 4.1, 2.1), xpd=F)
}

#### correct methylation data for PCs, sex and non-linearities
## Linearize methylation matrix with box-cox
cmet<-function(met,data){
  t(apply(met,1,function(xxx){
    p=xxx
    la=tryCatch({bc=boxcox(p~I(data$ag),plotit = F);bc$"x"[which.max(bc$y)]},error=function(e)1)
    if(la>0.05 | la<(-0.05)){(p^la-1)/la}else{log(p)}
}))}

#### Pacemaker
### Implementation of Pacemaker from https://github.com/NuttyLogic/EpigeneticPacemaker
## Takes in a methylation matrix met and a response variable vector "age" corresponding to the samples
## co 0.6 filters sites with correlation below value 0.6. Folds are the number of cross-validation folds
pacemaker<-function(met,age,co=0.6,folds=10){
  ## Filter sites below correlation or containing NAs
  ia=abs(apply(met,1,function(v)cor(as.numeric(v),age)))>co
  ib=sapply(ia,function(v)!any(is.na(v)))
  pacem=met[ia & ib,]
  agn=c()
  ## Initialize error differential
  nex=0
  di=length(age)/folds
  ## Classify samples that will belong to each fold
  sizes=rep(floor(di),folds);sizes[1]=sizes[1]+length(age)%%folds
  smp=sample(unlist(sapply(1:length(sizes),function(x)rep(x,sizes[x]))))
  ## Break methylation matrix and response variable into folds
  lX=as.list(tapply(1:length(age),smp,function(x)pacem[,x]))
  lY=tapply(1:length(age),smp,function(x)age[x])
  ## Loop over fold "e"
  for(e in 1:length(lX)){
    ## Separate fold into test (itself) and training (rest)
    tesX=as.matrix(lX[[e]]);tesY=lY[[e]];trX=do.call(cbind,lX[-e]);trY=unlist(lY[-e])
    ## Calculate mean methylation per site
    means=apply(trX,1,mean)
    ## For fold in test set, optimize regression up to 50 times or if error differential drops
    for(i in 1:50){
      ## Copied from script function find_solution_direct
      bets=apply(sweep(sweep(trX,1,means,"-"),2,trY-mean(trY),"*"),1,sum)/sum((trY-mean(trY))**2)
      ints=means-bets*mean(trY)
      ## Copied from script function update_states
      agu=(t(bets)%*%as.matrix(tesX)-sum(bets*ints))/sum(bets**2)
      ## Copied from script function calc_error
      ex=sqrt(sum((tesY-t(matrix(agu,ncol=1)%*%matrix(bets,nrow=1)))**2)/prod(dim(tesX)))
      ## Append ages calculated in test set
      if(abs(ex-nex)<0.0001){agn=c(agn,agu);break}
      nex=ex
    }
  }
  ## Return predicted ages in the same order they were introduced 
  return(as.vector(agn)[order(order(smp))])
}

### From lajoyce et al. unpublished
## Takes in a methylation matrix "dat" and a response variable vector "vari" corresponding to the samples
## co parameter filters sites with correlation below value co.

scAge=function(dat,vari,co,folds){
  ## Handle missing sites and filter by correlation value
  dat=dat[apply(dat,1,function(y){v=cor(y,vari);if(!is.na(v)){return(abs(v)>co & abs(sd(y)/mean(y))>0.05)}else{return(F)}}),]
  if(dim(dat)[1]==0){warning("No CpGs above rho");return(NA)}
  di=length(vari)/folds
  ## Create folds for data and response variable
  sizes=rep(floor(di),folds);sizes[1]=sizes[1]+length(vari)%%folds
  smp=sample(unlist(sapply(1:length(sizes),function(x)rep(x,sizes[x]))))
  lX=as.list(tapply(1:length(vari),smp,function(x)dat[,x]))
  lY=as.list(tapply(1:length(vari),smp,function(x)vari[x]))
  ## Number of samples of pseudo-probability normal 
  lo=100
  ## For any fold "fu" in parallel (if that works, probably not)
  ## Mclapply will return a list of folds of equal size, we cbind them together into a dataframe
  mi=do.call(cbind,mclapply(1:length(lX),function(fu){
    ## Separate data and response into test (fold) and training (rest)  
    tesX=as.matrix(lX[[fu]]);tesY=lY[[fu]];trX=do.call(cbind,lX[-fu]);trY=unlist(lY[-fu])
    ## Initialize result matrix with zeros for every data point in test set and the 100 samples
    m=matrix(rep(0,lo*length(tesY)),nrow=lo)
    ## For a site i in the training matrix X 
    for(i in 1:dim(trX)[1]){
      ## Calculate loess curve on training set
      ju=loess.sd(trX[i,],trY)
      ## Interpolate mean and SD age from test methylation value
      intm=approx(ju$x,ju$y,tesX[i,])$y
      intsd=approx(ju$x,ju$sd,tesX[i,])$y
      ## For every test data point mean and sd interpolated in the loess curve
      mu=do.call(cbind,mapply(function(intm,intsd){
        ## Calculate 100 normal samples with mean age  
        dd=dnorm(seq(min(vari)/1.5,max(vari)*1.5,length.out=lo),intm,intsd,log=T);
        ## Bound any infinite values and make NAs unlikely
        dd[is.na(dd)|dd<(-10)]=-10; dd
      },intm,intsd,SIMPLIFY = F))
      ## Add log site log probabilities together
      m=m+mu
    }
    m
  }))
  ## Return sample age likelihoods in the same order they were introduced
  ## (Before rearranging them for the folds)
  mi=mi[,order(order(smp))]
  ## Report maximum likelihood age
  seq(min(vari)/1.5,max(vari)*1.5,length.out=lo)[apply(mi,2,which.max)]
}

####Summary plot for axes of variation and variables
## dis is the distance matrix from methylation data
## age is a vector containing sample ages
## female is a vector containing 1s and 2s for male and female
## breed is the sample's breeds
## lsp is the moderator 
## phyl is the genetic distance matrix

## must be called in conjunction with a very specific layout outside the function
## bc I didn't know that  ape tree layouts and newplots had to be controlled from inside contMap/plot.phylo
trpl<-function(dis,age,female,breed,lsp,phyl){
  par(bg="transparent",new=T)
  ## make everything into ape format 
  ph=as.phylo(hclust(dis))
  phyl=as.phylo(hclust(as.dist(dist.from.cov(phyl))))
  ## normalize branch lengths
  ph$edge.length=norm2(ph$edge.length)
  phyl$edge.length=norm2(phyl$edge.length)
  names(age)=ph$tip.label
  names(lsp)=ph$tip.label
  ## Function contmap takes care of branch coloring
  q=contMap( ph, age,ftype="off",col="transparent",type="fan",res = 20, legend=F,plot=F)
  ## Change colors to match palette
  q$cols[1:length(q$cols)]=colorRampPalette(palette()[-1][1:5])(length(q$cols))
  par(xpd=NA)
  ## plot methylation tree
  plot(q,ftype="off",col="transparent",type="fan",res = 20, legend=T,mar=rep(4,4),leg.txt="Chron age",fsize=c(0.1,2.5))
  ## plot breed 
  ring(0.1,ph,col=RColorBrewer::brewer.pal(10,"Paired")[as.numeric(as.factor(breed))],offset=0.05)
  ## plot sex
  ring(0.1,ph,col=c("brown","orange")[female],offset=0.16)
  ## plot moderator
  ring(norm2(lsp)/3+0.02,ph,col=colorRampPalette(palette()[-1][1:5])(40)[cut(lsp,40)],offset = 0.28)
  par(xpd=NA)
  ## Make 
  q=contMap( phyl, age,ftype="off",col="transparent",type="fan",res = 20, legend=F,plot=F)
  ## Change colors to match palette
  q$cols[1:length(q$cols)]=colorRampPalette(palette()[-1][1:5])(length(q$cols))
  ## plot genetic tree
  plot(q,ftype="off",col="transparent",type="fan",res = 20, legend=F,mar=rep(5,4))
  ## plot breed
  ring(0.1,phyl,col=RColorBrewer::brewer.pal(10,"Paired")[as.numeric(as.factor(breed))],offset=0.05)
  ## plot sex
  ring(0.1,phyl,col=c("brown","orange")[female],offset=0.17)
  ## plot moderator
  ring(norm2(lsp)/2+0.02,phyl,col=colorRampPalette(palette()[-1][1:5])(40)[cut(lsp,40)],offset = 0.28)
}

####Process data

#################################
## Saliva capture (Pellegrini).##
#################################

## Read in methylation
wo=read.table("CAPTURE_DATASET",header=F)
metadata=read.table("CAPTURE_DATASET_METADATA",header=T)
metadata$sex=tolower(metadata$sex);metadata$spayed.neutered=tolower(metadata$spayed.neutered)
metadata$wieght=as.numeric(gsub("[^0-9,\\.]","",metadata$wieght,perl = T))
metadata$breed=tolower(metadata$breed)
### Change young to adult weights 
metadata$wieght[!is.na(metadata$spayed.neutered) & metadata$age<1 & metadata$breed=="french_bulldog"]=22
metadata$wieght[!is.na(metadata$spayed.neutered) & metadata$age<1 & metadata$breed=="pug"]=16
metadata$wieght[!is.na(metadata$spayed.neutered) & metadata$age<1 & metadata$breed=="rottweiler"]=90
metadata$wieght[!is.na(metadata$spayed.neutered) & metadata$age<1 & metadata$breed=="dachshund"]=24
metadata$wieght[!is.na(metadata$spayed.neutered) & metadata$age<1 & metadata$breed=="yorki_terrier"]=5
metadata$wieght[!is.na(metadata$spayed.neutered) & metadata$age<1 & metadata$breed=="caucasian_shepherd"]=105
wo=t(wo)
colnames(wo)=wo[1,];wo=wo[-1,]
wo=as.data.frame(wo)
nam=wo[,1:2]

## Read in SNPs
pelsnp=read.csv("SNPS_CAPTURE",header = T)
pelsnp=pelsnp[,apply(pelsnp,2,function(x)sum(x=="NaN"))/dim(pelsnp)[1]*100<26]
pelsnp=pelsnp[,match(colnames(wo),paste(gsub(".bam.txt","",pelsnp[1,]),".CGmap",sep = ""),nomatch = 0)]
wo=wo[,match(paste(gsub(".bam.txt","",pelsnp[1,]),".CGmap",sep = ""),colnames(wo),nomatch = 0)]
wo=apply(wo,2,as.numeric)
wo=cbind(nam,wo)
wo=wo[,!colnames(wo)%in%metadata[metadata$age==0,1]]
## Remove missing sites
wo=wo[apply(wo[,-c(1,2)],1,function(x)sum(is.na(x))/length(x)<0.2),]
## Remove missing samples
wo=wo[,apply(wo[,-c(1,2)],2,function(x)sum(is.na(x))/length(x))<.15]
pelsnp=pelsnp[,match(colnames(wo[,-c(1,2)]),paste(gsub(".bam.txt","",pelsnp[1,]),".CGmap",sep = ""),nomatch = 0)]
colnames(pelsnp)=colnames(wo)[-(1:2)]
pelsnp=pelsnp[-1,]
pelsnp=apply(pelsnp,2,as.numeric)
###Dubious imputation for missing methylation (very few)
wo[,-(1:2)]=t(apply(wo[,-(1:2)],1,function(x){x[is.na(x)]=median(x,na.rm=T);x}))
pc=prcomp(t(wo[,-c(1,2)]))
wo[,-c(1,2)]<-sweep(t(pc$x[,-1]%*%t(pc$rotation[,-1])),1,-pc$center)
metadata=metadata[match(colnames(wo)[-(1:2)],metadata[,1]),]
pc=prcomp(t(wo[,-(1:2)]))
#Infer sex
metadata$sex[match(rownames(pc$x),metadata[,1])]=c("female","male")[(cmdscale(dist(t(wo[wo[,1]=="chrX",-1:-2])),k=2)[,1]>0)+1]
metadata$wieght[is.na(metadata$wieght)]=71.5
## Keep good SNPs (HW compliant)
psn=pelsnp[apply(pelsnp,1,function(x){x=x[!is.na(x)];p=sum(x==0);q=sum(x);
  n=2*length(x);hw=log((sum(x==1)/n)/(2*p*q/n^2));abs(hw)<1 & min(p,q)>1}),]
## Impute 1 missing age
pl=cv.glmnet(y=sqrt(metadata$age[match(colnames(wo[,-c(1,2,3)]),metadata[,1])]),x = t(wo[,-c(1,2,3)]),nfolds = 10)
mod=lm(sqrt(metadata$age[match(colnames(wo[,-c(1,2,3)]),metadata[,1])])~t(wo[coefficients(pl,s=0.1353)[-1]>0,-c(1,2,3)]))
metadata$age[1]=(sum(coefficients(mod)[-1]*wo[coefficients(pl,s=0.1353)[-1]>0,3])-coefficients(pl,s=0.1353)[1])^2
#Calculate vcv from SNPs
psn=as.dist(apply(psn,2,function(x)apply(psn,2,function(y){sum(abs(x-y),na.rm = T)/sum(!(is.na(x) | is.na(y)))})))
psn=vcv(nj(psn))
## Generate GWAS model matrix (idle code)
correc3=eigen(psn)
correc3=(sqrt(abs(diag(correc3$values)))%*%correc3$vectors)[,1:15]
row.names(correc3)=row.names(psn)
correc3=cbind.data.frame(metadata$age,correc3[match(metadata[,1],row.names(correc3)),],as.numeric(as.factor(metadata$sex)),norm2(metadata$wieght))
colnames(correc3)=c("ag",paste0("V",1:(dim(correc3)[2]-3)),"Fem","lsp")
#Correct nonlinearities
cmpell=cmet(wo[,-1:-2],correc3[,-2:-16])

## decorrelating matrix
C=ginv(t(chol(psn)))
## Note NO lifespan clocks were used because no lifespan measures available
## Biological age clock capture (phylogeny) (Note, these are not used bc no lifespan measures)
abiopell=cv.glmnet(y=C%*%(metadata$age/metadata$wieght),x = C%*%t(cmpell),nfolds = 10)
## Biological age clock capture (no phylogeny) (Note, these are not used bc no lifespan measures)
abiopell2=cv.glmnet(y=(metadata$age/metadata$wieght),x = t(cmpell),nfolds = 10)
## Biological weight clock capture (phylogeny)
wbiopell=cv.glmnet(y=C%*%(metadata$age*metadata$wieght),x = C%*%t(cmpell),nfolds = 10)
## Biological weight clock capture (no phylogeny)
wbiopell2=cv.glmnet(y=(metadata$age*metadata$wieght),x = t(cmpell),nfolds = 10)
## Scage age clock capture (phylogeny)
ascpell=scAge(t(C%*%t(as.matrix(wo[,-1:-2]))),vari = C%*%(metadata$age/metadata$wieght),co=.5,folds = 10)
## Scage age clock capture (no phylogeny)
ascpell2=scAge(wo[,-1:-2],vari = metadata$age/metadata$wieght,co=.5,folds = 10)
## Scage weight clock capture (phylogeny)
wscpell=scAge(t(C%*%t(as.matrix(wo[,-1:-2]))),vari = C%*%(metadata$age*metadata$wieght),co=.5,folds = 10)
## Scage weight clock capture (no phylogeny)
wscpell2=scAge(wo[,-1:-2],vari = metadata$age*metadata$wieght,co=.5,folds = 10)

## Chronological age clock capture (phylogeny)
chronpell=cv.glmnet(y=C%*%metadata$age,x = C%*%t(cmpell),nfolds = 10)
## Chronological age clock capture (no phylogeny)
chronpell2=cv.glmnet(y=metadata$age,x = t(cmpell),nfolds = 10)
## Pacemaker weight clock capture (phylogeny)
agpell=pacemaker(t(C%*%t(as.matrix(wo[,-1:-2]))),C%*%as.vector(metadata$age)*metadata$wieght,co = .5,folds = 3)
## Pacemaker weight  clock capture (no phylogeny)
agpell2=pacemaker(as.matrix(wo[,-1:-2]),as.vector(metadata$age)*metadata$wieght,co = .5,folds = 3)
## Pacemaker weight clock capture (phylogeny)
wgpell=pacemaker(t(C%*%t(as.matrix(wo[,-1:-2]))),C%*%as.vector(metadata$age)*metadata$wieght,co = .5,folds = 3)
## Pacemaker weight  clock capture (no phylogeny)
wgpell2=pacemaker(as.matrix(wo[,-1:-2]),as.vector(metadata$age)*metadata$wieght,co = .5,folds = 3)


##Weight age GWAS capture (phylogeny)
wppel=t(gwas(as.data.frame(wo[,-1:-2]),psn,as.data.frame(correc3[,-2:-16]),quot = F))
##Weight age GWAS capture (no phylogeny)
wppel2=t(gwas(as.data.frame(wo[,-1:-2]),diag(dim(psn)[1]),as.data.frame(correc3[,-2:-16]),quot = F))
correc3=eigen(psn)
correc3=(sqrt(abs(diag(correc3$values)))%*%correc3$vectors)[,1:15]
row.names(correc3)=row.names(psn)
correc3=cbind.data.frame(metadata$age,correc3[match(metadata[,1],row.names(correc3)),],as.numeric(as.factor(metadata$sex)),metadata$wieght)
colnames(correc3)=c("ag",paste0("V",1:(dim(correc3)[2]-3)),"Fem","lsp")
##Bioage age GWAS capture (phylogeny)
appel=t(gwas(met = as.data.frame(wo[,-1:-2]),cova = psn,data = as.data.frame(correc3[,-2:-16]),quot = T))

#############
## Horvath ##
#############
## Read in methylation data
gs=readRDS("METHYLATION_ARRAY_GOOD_PROBES")
gs=gs[!(is.na(gs$probeStart) | is.na(gs$probeEnd)),]
load("METHYLATION_ARRAY_DATA")
colo=read.csv("DOG_BREED_EQUIVALENCES")
colo$DogBreed[colo$DogBreed=="Cocker Spaniel"]="Cocker Spaniel (American)"
colo$DogBreed[colo$DogBreed=="Mastiff"]="Mastiff (English)"
colo2=read.csv("DOG_BREED_EQUIVALENCES2")
x=normalized_betas_sesame
x=x[unlist(x[,1])%in%gs[,1],]
l=read.csv("METHYLATION_ARRAY_METADATA", header=T)
l$DogBreed[l$DogBreed=="Cocker Spaniel"]="Cocker Spaniel (American)"
l$DogBreed[l$DogBreed=="Mastiff"]="Mastiff (English)"
l$DogBreed[l$DogBreed=="Manchester Terrier"]="Manchester Terrier (Standard)"
l$DogBreed=colo2$DogBreed[match(l$DogBreed,colo2$DogBreed.old)]
l$lifespan_category=cut(l$LifespanMedianHorvath,c(0,10,12,20))
pro=scan("PROBE_ORDER",what="char")
x=x[match(pro,unlist(x[,1]),nomatch = 0),]
x=x[,c("yes",l$CanBeUsedForAgingStudies)=="yes"]
l=l[l$CanBeUsedForAgingStudies=="yes",]
### Read in breed tree from parker et al. 2017
tre=read.nexus("TREE_PARKER_2017")
tre=cophenetic(tre)
labs=sapply(strsplit(colnames(tre),"_"),function(x)x[1])
##Create median distance matrix from breeds
ju=by(tre,labs,function(x)by(t(x),labs,function(y)median(unlist(y))))
ju=sapply(names(ju),function(x)ju[[x]])
bn=read.table("RRBS_BREEDS")
bn[bn[,2]=="WEIM",2]="WHPG"
bn[bn[,2]=="ESET",2]="ISET"
bn[bn[,2]=="JACK",2]="GLEN"
exch=t(sapply(gsub(" ","",unique(l$DogBreed)),function(x){p=sapply(bn[,1],
        function(y)stringdist(x,y,method = "lv"));c(x,bn[which.min(p),1])}))
l$DogTree=bn[match(exch[match(gsub(" ","",l$DogBreed),exch[,1]),2],bn[,1]),2]
##Build vcv and broadcast to all samples in methylation dataset by breed
ju2=vcv(as.phylo(nj(as.dist(ju))))
ju3=ju2[match(l$DogTree,rownames(ju2)),match(l$DogTree,rownames(ju2))]
## Generate model matrix (idle)
correc=eigen(ju2)
correc=(sqrt(abs(diag(correc$values)))%*%correc$vectors)[,1:15]
row.names(correc)=row.names(ju2)
correc=cbind.data.frame(ag=l$Age,correc[match(l$DogTree,rownames(correc)),],l$Female+1,l$Weight.kg.avg)
colnames(correc)=c("ag",paste0("V",1:(dim(correc)[2]-3)),"Fem","lsp")
##Jitter vcv to make it non-singular 
R=diag(0.005,dim(ju3)[1])
R[upper.tri(R)]=rnorm(dim(R)[1]*(dim(R)[1]-1)/2,0,0.01)
R=R+t(R)
C=as.matrix(nearPD((ju3)**(1+R))[[1]])
# Create decorrelation matrix
C=solve(t(chol(C)))
#Correct nonlinearities
cmhor=cmet(x[,-1],correc[,-2:-16])
##Biolgical age weight clock mammalian (phylogeny)
wbiohor=cv.glmnet(y=C%*%(norm2(l$Age)*norm2(l$Weight.kg.avg)),x = C%*%t(cmhor),nfolds = 10)
##Biolgical age weight clock mammalian (no phylogeny)
wbiohor2=cv.glmnet(y=(norm2(l$Age)*norm2(l$Weight.kg.avg)),x = t(cmhor),nfolds = 10)
##Biolgical age clock mammalian ( phylogeny)
abiohor=cv.glmnet(y=C%*%(l$Age/l$Lifespan.AvgAKC),x = C%*%t(cmhor),nfolds = 10)
##Biolgical age clock mammalian (no phylogeny)
abiohor2=cv.glmnet(y=(l$Age/l$Lifespan.AvgAKC),x = t(cmhor),nfolds = 10)
##Chronological age clock mammalian (phylogeny)
chronhor=cv.glmnet(y=C%*%l$Age,x = C%*%t(cmhor),nfolds = 10)
##Chronological age clock mammalian (no phylogeny)
chronhor2=cv.glmnet(y=l$Age,x = t(cmhor),nfolds = 10)
##Pacemaker weight clock mammalian (phylogeny)
wghor=pacemaker(t(C%*%t(as.matrix(x[,-1]))),C%*%(l$Age*l$Weight.kg.avg),co = .5,folds = 10)
##Pacemaker weight clock mammalian (no phylogeny)
wghor2=pacemaker(as.matrix(x[,-1]),l$Age*l$Weight.kg.avg,co = .5,folds = 10)
##Pacemaker lifespan clock mammalian (phylogeny)
aghor=pacemaker(t(C%*%t(as.matrix(x[,-1]))),C%*%(l$Age/l$LifespanMedianHorvath),co = .5,folds = 10)
##Pacemaker lifespan clock mammalian (no phylogeny)
aghor2=pacemaker(as.matrix(x[,-1]),l$Age/l$LifespanMedianHorvath,co = .5,folds = 10)
##scage clock mammalian (phylogeny)
aschor=scAge(t(C%*%t(as.matrix(x[,-1]))),vari = C%*%(l$Age/l$LifespanMedianHorvath),co=.5,folds = 10)
##scage clock mammalian (no phylogeny)
aschor2=scAge(x[,-1],vari = l$Age/l$LifespanMedianHorvath,co=.5,folds = 10)
##scage weight clock mammalian (phylogeny)
wschor=scAge(t(C%*%t(as.matrix(x[,-1]))),vari = C%*%(l$Age*l$Weight.kg.avg),co=.5,folds = 10)
##scage weight clock mammalian (phylogeny)
wschor2=scAge(x[,-1],vari = l$Age*l$Weight.kg.avg,co=.5,folds = 10)

##GWAS biological weight age (phylogeny)
wphor=t(gwas(x[,-1],ju3,correc[,-2:-16],quot = F))
wphor[is.na(wphor)]=1
##GWAS biological weight age (no phylogeny)
wphor2=t(gwas(x[,-1],diag(dim(ju3)[1]),correc[,-2:-16],quot = F))
wphor2[is.na(wphor2)]=1
## Create model matrix (idle)
correc=eigen(ju3)
correc=(sqrt(abs(diag(correc$values)))%*%correc$vectors)[,1:15]
row.names(correc)=row.names(ju3)
correc=cbind.data.frame(ag=l$Age,correc[match(l$DogTree,rownames(correc)),],l$Female+1,l$LifespanMedianHorvath)
colnames(correc)=c("ag",paste0("V",1:(dim(correc)[2]-3)),"Fem","lsp")
##GWAS biological age (phylogeny)
aphor=t(gwas(met = x[,-1],cova = ju3,data = correc[,-2:-16],quot = T))
aphor[is.na(aphor)]=1

##########
## RRBS ##
##########
## Read in SNP data 
d=read.table("RRBS_SNPS",header=T,comment.char=" ")
d=d[d[,1]!="chrX",]
b=read.table("RRBS_AGES")
wei=read.table("RRBS_WEIGHTS")
b=cbind(b,rep(NA,length(b[,1])),rep(NA,length(b[,1])),wei[,2]);colnames(b)=paste0("V",1:length(b[1,]))
#Read in methyaltion data
a=read.table("RRBS_METHYLATION_MATRIX",header=T,comment.char=" ")
a=a[,!colnames(a)%in%c("DINGO_1","DINGO_2","X5M")]
a=a[,c(T,!grepl("^X",colnames(a)[-1]))]
a1=a[,1:3];a=a[,-c(1:3)]
d1=d[,1:2];d=d[,-c(1:2)]
d=d[,match(colnames(a),colnames(d),nomatch=0)]
b=b[match(colnames(d),b[,1]),]
aff=read.table("RRBS_SNPS",header=T, comment.char=" ")
aff[aff=="1/1"]=2;aff[aff=="0/1"]=1;aff[aff=="0/0"]=0;aff[aff=="./."]=NA;aff1=aff[,1:2]
aff=aff[,match(colnames(d),colnames(aff))]
aff=apply(aff,2,as.numeric)
gd2=apply(aff,2,function(x)apply(aff,2,function(y){p=abs(x-y);sum(p,na.rm=T)/sum(!is.na(p))}))
#Change NAs to idle values (very few)
a=t(apply(a,1,function(x){x[is.na(x)]=median(x,na.rm=T);x}))
#Change NAs to idle values (very few)
d=t(apply(d,1,function(x){x[is.na(x)]=round(median(x,na.rm=T),0);x}))
#Make distance matrix SNPs
gd=dist(t(d),method="manhattan")
md=dist(t(a))
pc=prcomp(t(a))
# Determine sex from meyhylation data
na=sweep(t(pc$x[,-c(1:2)]%*%t(pc$rotation[,-c(1:2)])),1,-pc$center)
xa=a[a1[,1]=="chrX",];xmd=dist(t(xa))
ann=cutree(hclust(as.dist(xmd)),k=2)
b[,3]=ifelse(ann==1,"F","M")
nu=cv.glmnet(t(a[,!is.na(b[,5])]),b[!is.na(b[,5]),5],nfold=5,alpha=0.5)
b[,6]=b[,5];b[is.na(b[,5]),6]=predict(nu,newx=t(a[,is.na(b[,5])]),s="lambda.min")
b[,7]=cut(b[,4],c(0,10.99,12.99,20))
#Create model matrix (idle)
correc2=cmdscale(gd,k=15)
correc2=cbind.data.frame(ag=b[,6],correc2[match(b[,1],rownames(correc2)),],as.numeric(as.factor(b[,3])),b[,8])
colnames(correc2)=c("ag",paste0("V",1:(dim(correc2)[2]-3)),"Fem","lsp")
#Correct nonlinearities
cmprrb=cmet(a,correc2[,-2:-16])
# Create decorrelation matrix
C=solve(t(chol(vcv(nj(gd)))))

##Biolgical age weight clock rrbs (phylogeny)
wbiorrb=cv.glmnet(y=C%*%(b[,6]*b[,8]),x = C%*%t(cmprrb),nfolds = 10)
##Biolgical age weight clock rrbs (no phylogeny)
wbiorrb2=cv.glmnet(y=(b[,6]*b[,8]),x = t(cmprrb),nfolds = 10)
##Biolgical age clock rrbs (phylogeny)
abiorrb=cv.glmnet(y=C%*%(b[,6]/b[,4]),x = C%*%t(cmprrb),nfolds = 10)
##Biolgical age clock rrbs (no phylogeny)
abiorrb2=cv.glmnet(y=(b[,6]/b[,4]),x = t(cmprrb),nfolds = 10)
##Chronological age clock rrbs (phylogeny)
chronrrb=cv.glmnet(y=C%*%b[,6],x = C%*%t(cmprrb),nfolds = 10)
##Chronological age clock rrbs (no phylogeny)
chronrrb2=cv.glmnet(y=b[,6],x = t(cmprrb),nfolds = 10)
##Pacemaker age clock rrbs (phylogeny)
wgrrb=pacemaker(t(C%*%t(a)),C%*%(b[,6]*b[,8]),co=0.5,folds = dim(a)[2])
##Pacemaker age clock rrbs (no phylogeny)
wgrrb2=pacemaker(a,b[,6]*b[,8],co=0.5,folds = dim(a)[2])
##Pacemaker age clock rrbs (phylogeny)
agrrb=pacemaker(t(C%*%t(a)),C%*%(b[,6]*b[,4]),co=0.5,folds = dim(a)[2])
##Pacemaker age clock rrbs (no phylogeny)
agrrb2=pacemaker(a,b[,6]*b[,4],co=0.5,folds = dim(a)[2])
##Scage clock rrbs (phylogeny)
ascrrbs=scAge(t(C%*%t(a)),vari = C%*%(b[,6]/b[,4]),co=.5,folds = dim(a)[2])
##Scage clock rrbs (no phylogeny)
ascrrbs2=scAge(a,vari = b[,6]/b[,4],co=.5,folds = dim(a)[2])
##Scage weight clock rrbs (no phylogeny)
wscrrbs=scAge(t(C%*%t(a)),vari = C%*%(b[,6]*b[,8]),co=.5,folds = dim(a)[2])
##Scage weight clock rrbs (no phylogeny)
wscrrbs2=scAge(a,vari = b[,6]*b[,8],co=.5,folds = dim(a)[2])

##GWAS biological weight age rrbs (phylogeny)
wprrb=t(gwas(a,vcv(nj(gd)),correc2[,-2:-16],quot = F))
wprrb[is.na(wprrb)]=1
wprrb=wprrb[wprrb[,1]>1E-30 & wprrb[,3]>1E-30,]
correc2=cmdscale(gd,k=15)
correc2=cbind.data.frame(ag=b[,6],correc2[match(b[,1],rownames(correc2)),],as.numeric(as.factor(b[,3])),b[,4])
colnames(correc2)=c("ag",paste0("V",1:(dim(correc2)[2]-3)),"Fem","lsp")
##GWAS biological age rrbs (phylogeny)
aprrb=t(gwas(a,vcv(nj(gd)),correc2[,-2:-16],quot = T))
aprrb[is.na(aprrb)]=1
aprrb=aprrb[aprrb[,1]>1E-30 & aprrb[,3]>1E-30,]

###########
## Plots ##
###########
#cccex=c(1.5,1.8,2.5)+0.4
cccex=c(1.5,1.8,2.5)
### Figure 2, basic stats 
cairo_pdf(file = "stats.pdf",height=16,width=13,bg = "transparent")
layout(matrix(ncol=2,c(1:4,rep(5,4)),byrow = T))

### Counts of male and female
par(xpd=NA,mar=c(5.1,5.1,4.1,2.1),bg=NA)
bp=c(table(l$Female+1),table((metadata$sex=="female")+1),table((b[,3]=="F")+1))
names(bp)=rep(c("Male","Female"),3)
barplot(bp,ylab="Counts",col=adjustcolor(rep(c("red","green","purple"),each=2),0.3),cex.axis=cccex[1],cex.names=cccex[2],cex.lab=cccex[2])
segments(c(0.1,2.55,4.9),c(-50,-50,-50),c(2.45,4.8,7.3),c(-50,-50,-50))
text(c(1.3,3.8,6.2),-rep(80,3),c("Mammalian","Capture","RRBS"),cex=cccex[3])

### Density of weights 
par(xpd=NA,bg=NA)
plot(1,col="transparent",xlim=c(-2,90),ylim=c(0,0.04),xlab="Weight (Kg)",ylab="Density",cex.axis=cccex[1],cex.names=cccex[2],cex.lab=cccex[2])
par(xpd=F,bg=NA)
polygon(c(0,density(l$Weight.kg.avg)$x,0),c(0,density(l$Weight.kg.avg)$y,0),
        col = adjustcolor("red",0.3),border = adjustcolor("black",0.3))
par(new=T)
polygon(c(0,density(metadata$wieght*0.4535)$x,0),c(0,density(metadata$wieght*0.4535)$y,0),
        col = adjustcolor("green",0.3),border = adjustcolor("black",0.3))
par(new=T)
polygon(c(0,density(b[,8])$x,0),c(0,density(b[,8])$y,0),
        col = adjustcolor("purple",0.5),border = adjustcolor("black",0.3))
legend("topright",c("Mammalian","Capture","RRBS"),fill = adjustcolor(c("red","green","purple"),.3),cex = cccex[2],bty='n',bg="transparent")

### Density of ages 
par(xpd=NA,bg="transparent")
plot(1,col="transparent",xlim=c(-2,20),ylim=c(0,0.14),xlab="Age (yrs)",ylab="Density",cex.axis=cccex[1],cex.names=cccex[2],cex.lab=cccex[2])
par(xpd=F,bg="transparent")
polygon(c(0,density(l$Age)$x,0),c(0,density(l$Age)$y,0),
        col = adjustcolor("red",0.3),border = adjustcolor("black",0.3))
par(new=T)
polygon(c(0,density(metadata$age)$x,0),c(0,density(metadata$age)$y,0),
        col = adjustcolor("green",0.3),border = adjustcolor("black",0.3))
par(new=T)
polygon(c(0,density(b[,6])$x,0),c(0,density(b[,6])$y,0),
        col = adjustcolor("purple",0.3),border = adjustcolor("black",0.3))
legend("topright",c("Mammalian","Capture","RRBS"),fill = adjustcolor(c("red","green","purple"),.3),cex=cccex[2],bty='n',bg="transparent")
## Weight vs lifespan
par(xpd=NA,bg="transparent")
plot(l$Weight.kg.avg[!l$DogTree%in%gsub("[_,0-9]","",b[,1])],l$LifespanMedianHorvath[!l$DogTree%in%gsub("[_,0-9]","",b[,1])],
     ylab="Lifespan (yrs)",xlab="Weight (Kgs)",cex.axis=cccex[1],cex.names=cccex[2],cex.lab=cccex[2])
par(xpd=F,bg="transparent")
points(b[,8][gsub("[_,0-9]","",b[,1])%in%l$DogTree],b[,4][gsub("[_,0-9]","",b[,1])%in%l$DogTree],pch=22,bg="black")
points(b[,8][!gsub("[_,0-9]","",b[,1])%in%l$DogTree],b[,4][!gsub("[_,0-9]","",b[,1])%in%l$DogTree],pch=25,bg="grey")
legend("topright",c("Mammalian only","Mammalian & RRBS","RRBS only"),pt.bg = c("White","black","gray"),pch=c(21,22,25),cex=cccex[2],bty='n',bg="transparent")
mod=summary(lm(c(l$Weight.kg.avg,b[,8])~c(l$LifespanMedianHorvath,b[,4])))
text(x=20,y=8,paste0("rho = ",round(mod$"adj.r.squared",2),"\np-val = 4.29e-77"),cex=cccex[2])

## Tree of breeds 
par(bg="transparent")
bre=read.csv("PARKER_2017_TREE", header = T)
nju=as.phylo(hclust(as.dist(ju)))
nju=drop.tip(nju,nju$tip.label[!nju$tip.label%in%unique(l$DogTree)])
par(xpd=NA)
l$GROUP=bre$FCI_GROUP[match(l$DogTree,bre$BREED_ABBR)]
plot.phylo(nju,type="fan",font=2,cex=cccex[1],x.lim=c(-2500,2000),y.lim=c(-2500,2000),
   tip.color = RColorBrewer::brewer.pal(10,"Paired")[as.numeric(as.factor(l$GROUP[match(nju$tip.label,l$DogTree)]))])
ring(nju,offset=400,x=500*norm2(l$LifespanMedianHorvath[match(nju$tip.label,l$DogTree)])+5,
   col=as.numeric(cut(l$LifespanMedianHorvath[match(nju$tip.label,l$DogTree)],4))+1)
legend("bottomright",levels(cut(l$LifespanMedianHorvath[match(nju$tip.label,l$DogTree)],4)),
   fill=(1:length(levels(cut(l$LifespanMedianHorvath[match(nju$tip.label,l$DogTree)],4))))+1,
   title = "Lifespan (yrs)",cex=cccex[2],bty='n')
legend("topleft",gsub("[,()].*","",levels(as.factor(l$GROUP[match(nju$tip.label,l$DogTree)]))),
   fill=RColorBrewer::brewer.pal(10,"Paired")[as.numeric(as.factor(levels(as.factor(l$GROUP[match(nju$tip.label,l$DogTree)]))))],
   title = "FCI Group",cex=cccex[1],bg = "transparent",bty='n')
text(2000,2000,paste0("Pagel's λ = ",round(phylosig(nju,l$LifespanMedianHorvath[match(nju$tip.label,l$DogTree)],method = "lambda")$lambda,2)),cex=2)
mtext(text = LETTERS[1:5],line = c(-2,-2,-34,-34,-64),at = c(rep(c(0.01,0.51),2),0.01),outer=T,font=2,cex=2)
dev.off()
##### END Figure 2 stats plot #####

#### Start Supplementary figure colocalization ####
library(gtools)

png("test_coloc.png",width=800,height=1000)
par(mfrow=c(3,1))
## Length chromosomes canfam4
ls=c(123556469,84979418,92479059,89535178,89562946,78113029,81081596,76405709,61171909,70643054,74805798,72970719,64299765,61112200,64676183,60362399,65088165,56472973,55516201,58627490,51742555,61573679,53134997,48566227,51730745,39257614,46662488,41733330,42517134,40643782,39901454,40225481,32139216,42397973,28051305,31223415,30785915,24803098,124992030)
m=mixedorder(paste(wo[,1],wo[,2],sep = ":"))
m[!grepl("^chr[0-9,X]+",wo[,1])]=0
wo1=wo[m,1:2]
dd1=as.numeric(wo1[,2])+unlist(mapply(function(x,y)rep(x,y),c(0,cumsum(ls)[-length(ls)]),rle(wo1[,1])$length))
## Plot p-values chronological age
par(mar=c(5.1,6.1,4.1,2.1))
plot(dd1,-log10(wppel[m,1]), xlim=c(0,2353542698),col=c("black","grey")[unlist(mapply(function(x,y){rep(x,y)},rep(1:2, 100)[1:length(rle(wo1[,2])$length)],rle(wo1[,2])$length))],
  type="l",axes=F,xlab="",ylab="",main="CAPTURE",cex.lab=cccex[2]+0.2, cex.axis=cccex[2], cex.main=cccex[3])
axis(side=4,cex.lab=cccex[2]+0.2, cex.axis=cccex[2])
par(new=T)
## Plot p-values interaction
plot(dd1,-log10(wppel[m,4]),xlim=c(0,2353542698),col=c("black","grey")[unlist(mapply(function(x,y){rep(x,y)},rep(1:2, 100)[1:length(rle(wo1[,1])$length)],rle(wo1[,1])$length))],
  cex=.6,pch=20,ylab=expression("-log"[10]~"pval"),xlab="Position (bps)",cex.lab=cccex[2]+0.2, cex.axis=cccex[2])
abline(v=dd1[order(wppel[m,4])[1:10]],lty="dashed",col=adjustcolor("red",0.3),lwd=3)
loc1=dd1[order(wppel[m,4])[1:20]]
## Lengths chromosomes Canfam 3.1
ls=c(122678785,85426708,91889043,88276631,88915250,77573801,80974532,74330416,61074082,69331447,74389097,72498081,63241923,60966679,64190966,59632846,64289059,55844845,53741614,58134056,50858623,61439934,52294480,47698779,51628933,38964690,45876710,41182112,41845238,40214260,39895921,38810281,31377067,42124431,26524999,30810995,30902991,23914537,123869142)
m=mixedorder(paste(a1[,1],a1[,2],sep = ":"))
m[!grepl("^chr[0-9,X]+",a1[,1])]=0
a3=a1[m,1:2]
dd2=as.numeric(a3[,2])+unlist(mapply(function(x,y)rep(x,y),c(0,cumsum(ls)[-length(ls)]),rle(a3[,1])$length))
## Plot p-values chronological age
par(mar=c(5.1,6.1,4.1,2.1))
plot(dd2,-log10(wprrb[m,1]), xlim=c(0,2327339098),col=c("black","grey")[unlist(mapply(function(x,y){rep(x,y)},rep(1:2, 100)[1:length(rle(a3[,2])$length)],rle(a3[,2])$length))],
     type="l",axes=F,xlab="",ylab="",main="RRBS",cex.lab=cccex[2]+0.2, cex.axis=cccex[2],cex.main=cccex[3])
axis(side=4,cex.lab=cccex[2]+0.2, cex.axis=cccex[2])
par(new=T)
## Plot p-values interaction
plot(dd2,-log10(wprrb[m,4]),xlim=c(0,2327339098),col=c("black","grey")[unlist(mapply(function(x,y){rep(x,y)},rep(1:2, 100)[1:length(rle(a3[,1])$length)],rle(a3[,1])$length))],
     cex=.6,pch=20,ylab=expression("-log"[10]~"pval"),xlab="Position (bps)",cex.lab=cccex[2]+0.2, cex.axis=cccex[2])
loc2=dd2[order(wprrb[m,4])[1:20]]
abline(v=dd2[order(wprrb[m,4])[1:10]],lty="dashed",col=adjustcolor("red",0.3),lwd=3)
pro=read.table("ORDERED_PROBED_RRBS",header=T)
pro=pro[!is.na(pro[,3]),]
m=mixedorder(paste(pro[,2],pro[,3],sep=":"))
m[!grepl("^chr[0-9,X]+",pro[m,2])]=0
pro=pro[m,]
dd3=as.numeric(pro[,3])+unlist(mapply(function(x,y)rep(x,y),c(0,cumsum(ls)[-length(ls)]),rle(pro[,2])$length))
## Plot p-values chronological age
par(mar=c(5.1,6.1,4.1,2.1))
plot(dd3,-log10(wphor[match(pro[,1],unlist(x[,1])),1]), xlim=c(0,2327339098),col=c("black","grey")[unlist(mapply(function(x,y){rep(x,y)},rep(1:2, 100)[1:length(rle(pro[,2])$length)],rle(pro[,2])$length))],type="l",axes=F,xlab="",ylab="",main="MAMMALIAN",cex.lab=cccex[2]+0.2, cex.axis=cccex[2],cex.main=cccex[3])
axis(side=4,cex.lab=cccex[2]+0.2, cex.axis=cccex[2])
par(new=T)
loc3=dd3[order(wphor[match(pro[,1],unlist(x[,1])),4])[1:20]]
## Plot p-values interaction
plot(dd3,-log10(wphor[match(pro[,1],unlist(x[,1])),4]), xlim=c(0,2327339098),col=c("black","grey")[unlist(mapply(function(x,y){rep(x,y)},rep(1:2, 100)[1:length(rle(pro[,2])$length)],rle(pro[,2])$length))],cex=.6,pch=20,ylab=expression("-log"[10]~"pval"),xlab="Position (bps)",cex.lab=cccex[2]+0.2, cex.axis=cccex[2])
abline(v=dd3[order(wphor[match(pro[,1],unlist(x[,1])),4])[1:10]],lty="dashed",col=adjustcolor("red",0.3),lwd=3)
mtext(side=1,at=0.02,line = c(-33,-68,-103),text = c("C","B","A"),outer=T,font=2,cex=2)
dev.off()
#### END colocalization figure ####

#### Start Figure 5, dual p-values ####
cairo_pdf("test_rr.pdf",width=15,height=15,bg="transparent")
par(oma = c(5,7,3,2),cex.axis=cccex[2]+0.1,cex.lab=cccex[2]+0.5)
layout(cbind(matrix(1:6,nrow=6),matrix(c(7:10,11,11),nrow=6)))
## Plot Mammalian
wphor[wphor==0]=min(wphor[,3])
pplot(plif_no_int=wphor[,3],pag_int=wphor[,4],pag_no_int=wphor[,1],bi=600,bio=wbiohor,chron=chronhor)
## Plot rrbs
wprrb[wprrb==0]=min(wprrb[,3])
pplot(wprrb[,3],wprrb[,4],wprrb[,1],bio=wbiorrb,chron=chronrrb,bi=300)
## Plot capture
wppel[wppel==0]=min(wppel[,3])
pplot(wppel[,3],wppel[,4],wppel[,1],bio=wbiopell,chron=chronpell,bi=300)
## Plot Lifespan
aphor[aphor==0]=min(aphor[,3])
pplot(plif_no_int=aphor[,3],pag_int=aphor[,4],pag_no_int=aphor[,1],bi=600,bio=abiohor,chron=chronhor)
aprrb[aprrb==0]=min(aprrb[,3])
pplot(aprrb[,3],aprrb[,4],aprrb[,1],bio=abiorrb,chron=chronrrb,bi=300)
mtext(side = 1,line=2,bquote(paste('-log'['10']*' (interaction)')), cex=2,outer = TRUE)
mtext(line = 0,at=c(0.25,.75),text=c("Weight","Lifespan"), cex=2, outer=T, font=2)
mtext(side=2,line=2,at=c(0.16,0.52,0.86),text=c("CAPTURE","RRBS","MAMMALIAN"), cex=2, outer=T)
mtext(line = -1,at=c(0.017,.5),text=c("A","D"), cex=2, outer=T, font=2)
mtext(line = -36,at=c(0.017,.5),text=c("B","E"), cex=2, outer=T, font=2)
mtext(line = -72,at=c(0.017),text=c("C"), cex=2, outer=T, font=2)
dev.off()
### END Figure5, dual p-vals ###


#cccex=cccex+0.3

### Begin Figure 3, Regressions ###
pdf("test_3_test.pdf",width=25,height=25)
#Remake vcv for mammalian
R=diag(0.005,dim(ju3)[1])
R[upper.tri(R)]=rnorm(dim(R)[1]*(dim(R)[1]-1)/2,0,0.01)
R=R+t(R)
C=as.matrix(nearPD((ju3)**(1+R))[[1]])
lay=c(rep(1,3),rep(2:4,each=2),4)
lay1=matrix(rep(c(5,13,21),each=4),ncol=2,byrow=T)
lay2=cbind(rbind(rep(seq(6,11,2),each=2),rep(seq(7,11,2),each=2)),matrix(rep(12,4),ncol=2))
lay3=cbind(rbind(rep(seq(14,19,2),each=2),rep(seq(15,19,2),each=2)),matrix(rep(20,4),ncol=2))
lay4=t(matrix(rep(22:25,each=4),byrow=T,ncol=2))
#layout(rbind(c(rep(1,each=3),2,3),matrix(rep(4:(5*3+3),each=2),byrow = T,nrow=3)),heights = c(1,3,3,3))
layout(rbind(matrix(lay,nrow=1),cbind(lay1,rbind(lay2,lay3,lay4))),heights = c(1.2,rep(3,6)))
par(oma = c(5,10,1,2),cex.axis=cccex[2]+0.3,cex.lab=cccex[2]+0.7,xpd=NA,cex.main=cccex[3]-0.1)
## Make legends
plot(1.2,2,col="transparent",axes=F,xlab="",ylab="",xlim=c(1,2),ylim=c(0,2));
legend(1,7,legend = c("Low","Medium-Low","Medium-High","High"),fill = 2:5,horiz = T,title = "Moderator",cex=cccex[2]+0.5,bty='n',bg="transparent")
plot(1,1,col="transparent",axes=F,xlab="",ylab="",xlim=c(0,2),ylim=c(0,2));
legend(1.5,7,legend = c("Lifespan","Weight"),pch=c(15,20),col="black",pt.cex = 2/1.2*1.6,horiz = F,title = "Moderator",cex=cccex[2]+0.5,bty='n',bg="transparent")
plot(1,1,col="transparent",axes=F,xlab="",ylab="",xlim=c(0,2),ylim=c(0,2));
legend(1.5,7,legend = c("1e-5","1e-50","1e-100"),pch=21,col="black",pt.cex = c(1,2,3)/1.2*1.4,horiz = F,title = "p-value",cex=cccex[2]+0.5,bty='n',bg="transparent")
legend(1.35,7,legend = c("","",""),pch=22,col="black",pt.cex = c(1,2,3)/1.2*1.4,horiz = F,title = "",cex=cccex[2]+0.5,bty='n',bg="transparent")
plot(1,1,col="transparent",axes=F,xlab="",ylab="",xlim=c(0,2),ylim=c(0,2));
legend(1.2,7,legend = c("OLS","GLS"),fill = c("red","blue"),horiz = F,title = "Regression type",cex=cccex[2]+0.5,bty='n',bg="transparent")
par(xpd=F)
# ## Plot mammalian
# # plreg(met = cmhor,chron2 = chronhor2,chron = chronhor,bio2=wbiohor2,bio=wbiohor,
# #       pace2=aghor,pace=aghor2,age = l$Age,vari = l$Weight.kg.avg,vcvv=ginv(C),quot=F,
# #       scAge = wschor2,scAge2 = wschor,sparse=T,pgls=T)
# plreg2(met = cmhor,chron2 = chronhor2,chron = chronhor,bio2=wbiohor2,bio=wbiohor,
#       pace2=wghor,pace=wghor2,age = l$Age,vari = l$Weight.kg.avg,vari2 = l$LifespanMedianHorvath,vcvv=ginv(C),
#       scAge = wschor2,scAge2 = wschor,sparse=T,pgls=T,abio2=abiohor2,abio=abiohor,
#       apace2=aghor,apace=aghor2,ascAge = aschor2,ascAge2 = aschor)
# 
# # ## Plot rrbs
# # plreg(met = cmprrb,chron = chronrrb,chron2 = chronrrb2,bio=wbiorrb,bio2=wbiorrb2,
# #       pace=agrrb2,pace2=agrrb,age = b[,6],vari = b[,8],vcvv = ginv(vcv(nj(gd))),quot = F,
# #       scAge2 = wscrrbs,scAge = wscrrbs2,sparse=T,pgls=T)
# plreg2(met = cmprrb,chron = chronrrb,chron2 = chronrrb2,bio=wbiorrb,bio2=wbiorrb2,
#        abio=abiorrb,abio2=abiorrb2,apace=agrrb2,apace2=agrrb,ascAge2 = ascrrbs,ascAge = ascrrbs2,
#        pace=wgrrb2,pace2=wgrrb,age = b[,6],vari = b[,8],vari2=b[,4],vcvv = ginv(vcv(nj(gd)))
#        ,scAge2 = wscrrbs,scAge = wscrrbs2,sparse=T,pgls=T)
# ## Plot capture
# plreg(met = cmpell,chron = chronpell,chron2 = chronpell2,bio=wbiopell,bio2=wbiopell2,pace2=agpell,
#       pace=agpell2,age = metadata$age,vari = metadata$wieght,vcvv=ginv(psn),quot=F,
#       scAge2 = wscpell,scAge = wscpell2,sparse=T,pgls=T)
###########
# plreg(met = cmhor,chron2 = chronhor2,chron = chronhor,bio2=wbiohor2,bio=wbiohor,
#       pace2=aghor,pace=aghor2,age = l$Age,vari = l$Weight.kg.avg,vcvv=ginv(C),quot=F,
#       scAge = wschor2,scAge2 = wschor,sparse=F,pgls=F)
plreg2(met = cmhor,chron2 = chronhor2,chron = chronhor,bio2=wbiohor2,bio=wbiohor,
       pace2=wghor,pace=wghor2,age = l$Age,vari = l$Weight.kg.avg,vari2 = l$LifespanMedianHorvath,vcvv=ginv(C),
       scAge = wschor2,scAge2 = wschor,sparse=F,pgls=F,abio2=abiohor2,abio=abiohor,
       apace2=aghor,apace=aghor2,ascAge = aschor2,ascAge2 = aschor)

# plreg(met = cmprrb,chron = chronrrb,chron2 = chronrrb2,bio=wbiorrb,bio2=wbiorrb2,
#       pace=agrrb2,pace2=agrrb,age = b[,6],vari = b[,8],vcvv = ginv(vcv(nj(gd))),quot = F,
#       scAge2 = wscrrbs,scAge = wscrrbs2,sparse=F,pgls=F)
plreg2(met = cmprrb,chron = chronrrb,chron2 = chronrrb2,bio=wbiorrb,bio2=wbiorrb2,
       abio=abiorrb,abio2=abiorrb2,apace=agrrb2,apace2=agrrb,ascAge2 = ascrrbs,ascAge = ascrrbs2,
       pace=wgrrb2,pace2=wgrrb,age = b[,6],vari = b[,8],vari2=b[,4],vcvv = ginv(vcv(nj(gd)))
       ,scAge2 = wscrrbs,scAge = wscrrbs2,sparse=F,pgls=F)

plreg(met = cmpell,chron = chronpell,chron2 = chronpell2,bio=wbiopell,bio2=wbiopell2,pace2=agpell,
       pace=agpell2,age = metadata$age,vari = metadata$wieght,vcvv=ginv(psn),quot=F,
       scAge2 = wscpell,scAge = wscpell2,sparse=F,pgls=F)
mtext(side = 1, "Chronological Age", cex=cccex[2]+0.4,outer = TRUE)
mtext(side = 2, "Epigenetic Age (arbitrary units)", cex=cccex[2]+0.4,outer = TRUE,padj = -3)
mtext(side=2,line=1,at=c(0.18,0.5,0.78)+.03,text=c("CAPTURE","RRBS","MAMMALIAN"), cex=2, outer=T)
mtext(line = -13,at=c(0.015,0.815),text=c("A","B"), cex=cccex[2]+0.4, outer=T, font=2)
mtext(line = -70,at=c(0.015,0.815),text=c("C","D"), cex=cccex[2]+0.4, outer=T, font=2)
mtext(line = -128,at=c(0.015,0.815),text=c("E","F"), cex=cccex[2]+0.4, outer=T, font=2)
dev.off()
### END Figure 4, Regressions ###


### Begin Figure 3, Methylation and SNP trees ###
pdf("test_1.pdf",width=25,height=20)
#layout(matrix(c(1,1,4,4,7,7,2,3,5,6,8,9),ncol=2))
layout(matrix(c(1,1,4,4,2,3,5,6,7,7,10,10,8,9,10,10),nrow=4))

par(oma=rep(5,4))
### Spacer blank plot
plot(0,0,col="transparent",xlab="",ylab="",axes=F)
text(-1,1,"A",font=2,cex=5)
par(new=T)
colnames(C)=colnames(x[,-1]);rownames(C)=colnames(x[,-1])
## Plot mammalian
#trpl(dis = dist(t(x[,-1])),age = l$Age,female = l$Female+1,breed = l$GROUP,lsp=l$Weight.kg.avg,phyl=C)
trpl(dis = ddd,age = l$Age,female = l$Female+1,breed = l$GROUP,lsp=l$Weight.kg.avg,phyl=C)
plot(0,0,col="transparent",xlab="",ylab="",axes=F)
uu=unique(l$GROUP)
uu[is.na(uu)]="Mixed Breed or Unknown"
legend("center",c(gsub("[,()].*","",c(levels(as.factor(l$GROUP[match(nju$tip.label,l$DogTree)])))),"Unknown or Mixed breed"),
       fill=c(RColorBrewer::brewer.pal(10,"Paired")[as.numeric(as.factor(levels(as.factor(l$GROUP[match(nju$tip.label,l$DogTree)]))))],"white"),
       title = "FCI Group",cex=cccex[2]+0.5,bg="transparent",bty="n")
### Spacer blank plot
plot(0,0,col="transparent",xlab="",ylab="",axes=F)
text(-1,1,"B",font=2,cex=5)
par(new=T)
uu=sapply(metadata$breed,function(x){a=gsub(" ","_",tolower(l$DogBreed));l$GROUP[which.min(stringdist(a=x,b=a))]})
uu[grepl("[?,0-9+]|mix",metadata$breed)]=NA
### Plot capture
trpl(dis = dist(t(wo[,-1:-2])),age = metadata$age,female = 1+(metadata$sex=="female"),
     breed = uu,lsp=metadata$wieght*0.4536,phyl=psn)
plot(0,0,col="transparent",xlab="",ylab="",axes=F)
legend("center",legend = c("Male","Female"),fill=c("brown","orange"),title="Sex",
       cex=cccex[3]+0.6,bg="transparent",bty="n")
## Spacer blank plot
plot(0,0,col="transparent",xlab="",ylab="",axes=F)
par(new=T)
text(-1,1,"C",font=2,cex=5)
### Plot rrbs
trpl(dis = dist(t(a)),age = b[,6],female = 1+(b[,3]=="F"),
     breed = l$GROUP[match(gsub("[0-9]?_*","",b[,1]),l$DogTree)],lsp=b[,8],phyl=vcv(nj(gd)))
plot(0,0,col="transparent",xlab="",ylab="",axes=F)
legend("center",legend = sort(unique(cut(l$Weight.kg.avg,5))),
       fill=sort(unique(as.numeric(cut(c(l$Weight.kg.avg,metadata$wieght*0.4536,b[,8]),5))))+1,
       title="Weight (Kg)",cex=cccex[2]+0.6,bg="transparent",bty="n")
### Horrible  venn diagram made manually ###
par(xpd=NA,mar=rep(1,4))
plot(1,1,xlim=c(0.15,.85),ylim=c(0.15,.85),col="transparent",xlab="",ylab="",axes=F)
points(c(0.35,0.5,0.54),c(0.35,.6,0.25),cex=c(20,33,27)*3.4,pch=21,bg=adjustcolor(c("red","blue","green"),.4),col="transparent")
text(c(0.44,0.415,0.4,0.52,0.34,0.53,0.5),c(0.39,0.285,0.44,0.41,0.31,0.23,.55),c("200","2,442","4,771","32,492","30,930","244,334","~511,000"),col="black",font=2,cex=cccex[2]+0.9)
text(c(0.34,0.55,0.5),c(0.37,0.29,.6),c("Mammalian","RRBS + Capture","EPIC"),col="black",font=2,cex=cccex[3]+0.4)
text(0.2,0.79,"D",font=2,cex=5)
dev.off()
### END Figure 3, Methylation and SNP trees ###

### Begin Supplementary figures, cross-validation  and dissection of coefficients###
## Takes in chronological and bioage glment objects
cvplot=function(bio1,bio2,titl){
    par(mar=rep(5,4))
    col1=adjustcolor("blue",.5);col2=adjustcolor("red",.5)
    ##extract number of parameters and values from glmet object 1
    plot(bio1$nzero,bio1$cvm,col=col1,pch=20,axes=F,xlab="",ylab="",ylim=c(min(bio1$cvlo),max(bio1$cvup)))
    arrows(x0=bio1$nzero,x1=bio1$nzero,y1=bio1$cvup,y0=bio1$cvlo, code=3, angle=90, length=0.1,col=col1)
    mi1=round(seq(min(bio1$nzero),max(bio1$nzero),length.out=5),0)
    mi2=seq(min(bio1$cvlo),max(bio1$cvup),length.out=5)
    axis(1,at = mi1,labels = mi1,col=col1,col.axis=col1)
    axis(2,at = mi2,labels = format(mi2,scientific=T,digits=3),col=col1,col.axis=col1)
    abline(v = bio1$nzero[bio1$index],col=col1,lty="dashed",lwd=3)
    par(new=T)
    ##extract number of parameters and values from glmet object 2
    plot(bio2$nzero,bio2$cvm,col=col2,pch=20,axes=F,xlab="",ylab="",ylim=c(min(bio2$cvlo),max(bio2$cvup)))
    arrows(x0=bio2$nzero,x1=bio2$nzero,y1=bio2$cvup,y0=bio2$cvlo, code=3, angle=90, length=0.1,col=col2)
    mi1=round(seq(min(bio2$nzero),max(bio2$nzero),length.out=5),0)
    mi2=seq(min(bio2$cvlo),max(bio2$cvup),length.out=5)
    axis(3,at = mi1,labels = mi1,col=col2,col.axis=col2)
    axis(4,at = mi2,labels = format(mi2,scientific=T,digits=3),col=col2,col.axis=col2)
    abline(v = bio2$nzero[bio2$index],col=col2,lty="dashed",lwd=3)
    title(main = titl,line=2.5);
    mtext(text = c("# Coeff PGLS","MSE PGLS","# Coeff OGLS","MSE OGLS"),side=1:4,col = c(col1,col1,col2,col2),line = 3.5,cex=1.5)
}

#### Plot coefficient values against respective p-values and specific regression values
cairo_pdf("dissection.pdf",height=12,width=8)
par(mfrow=c(4,2),cex.axis=cccex[2]+0.1,cex.lab=cccex[2]+0.5,mar=c(5.1,6.1,4.1,2.1))
library(TeachingDemos)
lam="lambda.min"
coe=abs(coefficients(wbiohor,s=lam)[-1])>0
plot(abs(coefficients(wbiohor,s=lam)[-1])[coe],pch=20,
  -log10(wphor[coe,4]),xlab="coefficients PGLS weight",ylab=expression("-log"[10]~"(p-value"[interaction]~")"))   
mtext(LETTERS[1],side = 1,at = -0.1,line=-15,font=2,cex=2)
pts=c(4,94,207,sample((1:250)[-c(4,94,207)],4))
shadowtext(x = abs(coefficients(wbiohor,s=lam)[-1])[which(coe)[pts]],
    y = -log10(wphor[which(coe)[pts],4]),col="white",labels = LETTERS[2:(1+length(pts))],cex=1.5,font=2)
sapply(1:length(pts),function(x){
  p=summary(lm(unlist(cmhor[which(coe)[pts][x],])~l$Age*l$Weight.kg.avg))$coefficients
  p2=anova(lm(unlist(cmhor[which(coe)[pts][x],])~l$Age),
    lm(unlist(cmhor[which(coe)[pts][x],])~l$Age+l$Age*l$Weight.kg.avg))$"Pr(>F)"[2]
  plot(l$Age,unlist(cmhor[which(coe)[pts][x],]),col="black",bg=as.numeric(cut(l$Weight.kg.avg,4))+1,
    xlab="Age",ylab="Linearized methylation",pch=21)
  title(bquote("p"[int]~"=" ~ .(format(p[dim(p)[1],dim(p)[2]],scientific=T,digits=2)) ~"; p"[Δint]~"=" ~.(format(p2,scientific=T,digits=2))),cex=2.5)
  mtext(LETTERS[1+x],side = 1,at = -2,line=-15,font=2,cex=2)
})
dev.off()

###### Plot cross-validation plots
pdf("cross-validation.pdf",height=20,width=11)
par(mfrow=c(3,2),lwd=2,oma=c(2,5,5,2),cex.axis=cccex[2]+0.1,cex.lab=cccex[2]+0.5)
cvplot(wbiohor,wbiohor2,"");cvplot(abiohor,abiohor2,"")
cvplot(wbiorrb,wbiorrb2,"");cvplot(abiorrb,abiorrb2,"")
cvplot(wbiopell,wbiopell2,"")
plot.new();legend("center",fill = c("blue","red"),legend = c("PGLS","OGLS"),cex=3)
mtext(side=2,line=2,at=0.05+c(0.13,0.47,0.79),text = c("CAPTURE","RRBS","MAMMALIAN"),outer=T, cex=2)
mtext(line = 1,at=c(0.25,0.75),text=c("Weight","Lifespan"), cex=2, outer=T, font=2)
mtext(line = -2,at=c(0.015,0.5),text=c("A","B"), cex=2, outer=T, font=2)
mtext(line = -33,at=c(0.015,0.5),text=c("C","D"), cex=2, outer=T, font=2)
mtext(line = -63,at=c(0.015),text=c("E"), cex=2, outer=T, font=2)
dev.off()
### Begin Supplementary figure, cross-validation  and dissection of coefficients###


##### Analysis of GO enrichment
library(WebGestaltR)

HA=read.csv("METHYLATION_ARRAY_MANIFESTO")

g1=HA[HA[,2]%in%unlist(x[row.names(x)%in%as.character(which(order(norm2(wphor[,1])*norm2(wphor[,4]))%in%(1:100))),1]),3]
g1=HA[HA[,2]%in%unlist(x[row.names(x)%in%names(sort(wphor[head(order(wphor[,4]),1000),1]))[1:100],1]),3]
##Beware pellegrini in canfam 4!
g2=wo[row.names(wo)%in%names(sort(wppel[head(order(wppel[,4]),1000),1]))[1:100],1:2]
g2[,2]=as.numeric(g2[,2])-5000;g2=cbind(g2,g2[,2]+10000)
### Genes looked up in UCSC
g2=c("ASIC2","BTRC","C2CD2","C6H7orf50","CACNA2D3","CDK15","CEACAM23","CELF4","CLVS1","CTNNA2","CTNND2","DPYD","GPER1","ITGA11","KRBA2","LOC100686845","LOC100856635","LOC102151182","LOC102153193","LOC102153620","LOC102153911","LOC102154225","LOC102154968","LOC102155397","LOC102155449","LOC102155842","LOC102155920","LOC102156970","LOC106557691","LOC106559309","LOC106560027","LOC111091469","LOC111092136","LOC111092856","LOC111093083","LOC111094013","LOC111094349","LOC111095144","LOC111098954","LOC475957","LOC610636","LSAMP","MICAL3","MRPL10","OSBPL7","OTP","PAX9","PTPRG","RIN2","RIOK1","RPL26","SETD5","SLC25A21","SLC2A10","SLC38A6","SLF2","ST6GALNAC5","STAM","TMEM236","TRERF1","WDR41","ZBTB16")
#g2=c("ASIC2","BARHL1","C1H19orf12","C1QL3","C1QTNF7","CADPS","CBLN1","CDH22","CTNNA2","CTNND2","CUX2","EN1","FOXD3","FRMD4A","FZD7","GRIK3","HOXB7","HOXB8","HOXC10","HOXC11","HOXC4","HOXD10","IHH","INSM1","LARGE1","LHX6","LMO1","LOC102152454","LOC102152918","LOC102155469","LOC102156069","LOC111090286","LOC111090473","LOC111091843","LOC111092797","LOC111092969","LOC111093616","LOC111093850","LOC111094348","LOC119864218","LOC119865116","LOC119865451","LOC119865610","LOC119866157","LOC119867493","LOC119867649","LOC119867963","LOC119870328","LOC119870818","LOC119871239","LOC119871943","LOC119872371","LOC119872627","LOC119874319","LOC119876050","LOC119876746","LOC119877729","LRRTM1","MAF","MECOM","MIR137","MIR196A-2","NR5A1","NTRK2","OLIG2","OTP","PAX5","PBX3","PHOX2B","PITX2","PLD1","PRDM13","PRKG1","PRR5","PSD2","PTER","RSU1","RYR2","SALL1","SLC12A5","SOX2","TBX5","TFAP2A","TFAP2D","UNCX","URI1","VAX1","ZBTB20","ZIC3","ZNF235","ZZZ3")
a=data.frame(a);row.names(wprrb)=1:dim(wprrb)[1]
g3=a1[row.names(a)%in%names(sort(wprrb[head(order(wprrb[,4]),1000),1]))[1:100],1:3]
g3[,2]=g3[,2]-5000;g3[,3]=g3[,3]+5000;
g3=c("ARNT2","ASIC2","C12H6orf132","C1QL3","CACNA2D3","CELF4","CHD7","CTNNA2","CTNND2","DES","EN1","EPHA10","FEZF2","FOXD3","FOXG1","FZD7","HES7","HOXA11","HOXB3","HOXC11","HOXC4","HOXC5","HOXC6","HOXC8","HOXD10","HOXD9","IHH","INSM1","ISL1","KCNMA1","LBX1","LHFPL4","LHX6","LOC102151514","LOC102152853","LOC102153255","LOC102154721","LOC102154783","LOC102154791","LOC102156194","LOC106559786","LOC111091244","LOC111091843","LOC111094337","LOC111094348","LOC111095612","LOC111096999","LOC119864218","LOC119865610","LOC119867346","LOC119867649","LOC119867963","LOC119870328","LOC119871811","LOC119872371","LOC119872799","LOC119873653","LOC119874164","LOC119874319","LOC119876645","LOC119877729","LRRTM1","MEIS2","MIR137","NKX1-1","NR2F1","NR2F2","NRN1","OLIG2","OTP","PAX5","PHOX2B","POU3F3","PRDM13","PSD2","PTER","RSU1","SALL1","SKOR1","SLC12A5","TFAP2A","TFAP2D","UNC79","UNCX","ZBTB16","ZBTB20","ZIC1","ZIC4","ZNF235","ZZZ3")
#g3=c("ADCY1","CAB39L","CCDC191","CCDC88C","CDV3","CLIC6","COL5A2","COR12J1","CYSLTR2","DHCR7","EGFR","ELFN1","EMID1","ENC1","EVI5L","FAM193A","FAM53B","FAT3","FGFR2","FLYWCH1","GPHN","GRB7","GTF3C4","HES7","HOXA7","HOXA9","HTT","IL16","ITGA2B","ITGB2","ITIH1","KIAA1671","KSR2","LDB2","LOC102152259","LOC102152479","LOC102156374","LOC106558635","LOC106559534","LOC111090590","LOC111092337","LOC111094545","LOC111096724","LOC111097964","LOC488626","LOC488951","LOC491550","LOC607201","LOC612320","LRRC4C","LYPD3","MEGF8","MSANTD1","MSRA","MYT1L","NCAPD3","OBSCN","OR10J5","PIGQ","PKDCC","PPARGC1B","PRSS12","PSTPIP1","PYGO2","RAB40C","RAPGEF4","RTL1","SHC1","SLC25A47","SLC37A2","SLC9A3","SLCO2B1","TGFB2","TLE3","TMEM145","TOX2","TRIM35","TRIM37","VEPH1","VWA5B1","WARS","WDFY2","XAB2","XXYLT1","ZNF704")
WebGestaltR(enrichMethod="ORA", organism="hsapiens",enrichDatabase="geneontology_Biological_Process_noRedundant",
  interestGene = g1,referenceGene=unique(HA[,3]),interestGeneType="genesymbol",referenceGeneType="genesymbol")
WebGestaltR(enrichMethod="ORA", organism="hsapiens",enrichDatabase="geneontology_Biological_Process_noRedundant",
  interestGene = g2,referenceGene=unique(scan("GENE_SET_CAPTURE",what="c")),interestGeneType="genesymbol",referenceGeneType="genesymbol")
WebGestaltR(enrichMethod="ORA", organism="hsapiens",enrichDatabase="geneontology_Biological_Process_noRedundant",
  interestGene = g3,interestGeneType="genesymbol",referenceSet = "genome_protein-coding")

####Mboning et al example Figure
pdf("example_scage.pdf",9,6)
## number of samples
nnn=1000
## Make methylation 
xxx=rnorm(nnn,5,1)
## Make heteroskedastic error
eee=rnorm(nnn,0,seq(.1,2,length.out=nnn))
tes=c(2.5,5,7)
## Make non-linear y
yyy=abs(min(scale((xxx+eee)^2)))+scale((xxx+eee)^2)
layout(matrix(1:2,nrow = 1),widths = c(2,1))
## Loess 
loe=loess.sd(xxx,yyy)
par(mar=c(4.1,5.1,2.1,0))
### Plot!
plot(xxx/10,yyy,pch=20,xlab="Methylation",ylab="Chronological Age",frame.plot = F,ylim=c(min(yyy)-0.1,max(yyy)+0.1),cex.axis=1.1,cex.lab=1.3, main="Site X",cex.main=1.5)
lines(loe$x/10,loe$y+2*loe$sd,lwd=2)
lines(loe$x/10,loe$y-2*loe$sd,lwd=2)
lines(loe$x/10,loe$y,lwd=2,lty="dashed",col="red")
tesY=approx(loe$x,loe$y,xout = tes)$y
tesSD=approx(loe$x,2*loe$sd,xout = tes)$y
sapply(c(1,1.5,2),function(x)points(tes/10,tesY,col=rainbow(3),cex=x))
legend("topleft",c("Train",paste0("Test",1:3)),fill = c("black",rainbow(3)),bty='n')
legend("top",c("LOESS mean","LOESS sd"),col=c("red","black"),lty = c("dashed","solid"),bty='n')
segments(x0 = tes/10,x1=tes/10,y0=tesY-tesSD,y1=tesY+tesSD,col=rainbow(3),lwd=1.2)
dists=mapply(function(x,y)density(rnorm(10000,x,y)),tesY,tesSD/2,SIMPLIFY = F)
par(mar=c(4.1,0,2.1,2))
plot(0,0,col="transparent",ylim=c(min(yyy)-0.1,max(yyy)+0.1),xlim=c(0,.7),axes=F,xlab="Age likelihood",ylab="",cex.axis=1.1,cex.lab=1.3)
axis(4,cex.axis=1.1)
axis(1,cex.axis=1.1)
mapply(function(x,y)polygon(x =y$y,y=y$x,col=adjustcolor(x,.5),border = x),rainbow(3),dists,SIMPLIFY = F)
dev.off()
                                      
# ## Same as before, but lifespan
# pdf("test_2.pdf",width=15,height=20)
# R=diag(0.005,dim(ju3)[1])
# R[upper.tri(R)]=rnorm(dim(R)[1]*(dim(R)[1]-1)/2,0,0.01)
# R=R+t(R)
# C=as.matrix(nearPD((ju3)**(1+R))[[1]])
# layout(rbind(c(rep(1,each=3),2,3),matrix(rep(4:23,each=2),byrow = T,nrow=4)),heights = c(1,3,3,3,3))
# par(oma = c(5,16,1,5),cex.axis=cccex[2]+0.1,cex.lab=cccex[2]+0.5,xpd=NA,cex.main=cccex[3]-0.3)
# plot(1.2,1,col="transparent",axes=F,xlab="",ylab="",xlim=c(1,2),ylim=c(0,2));
# legend(1,1,legend = c("Low","Medium-Low","Medium-High","High"),fill = 2:5,horiz = T,title = "Moderator",cex=cccex[2],bty='n',bg="transparent")
# plot(1,1,col="transparent",axes=F,xlab="",ylab="",xlim=c(0,2),ylim=c(0,2));
# legend(1.2,4,legend = c("OLS","GLS"),fill = c("red","blue"),horiz = F,title = "Regression type",cex=cccex[2],bty='n',bg="transparent")
# plot(1,1,col="transparent",axes=F,xlab="",ylab="",xlim=c(0,2),ylim=c(0,2));
# legend(1.5,4,legend = c("1e-5","1e-50","1e-100"),pch=21,col="black",pt.cex = c(1,2,3)/1.2*1.4,horiz = F,title = "p-value",cex=cccex[2],bty='n',bg="transparent")
# par(xpd=F)
# plreg(met = cmhor,chron = chronhor,chron2 = chronhor2,bio=abiohor,bio2=abiohor2,
#       pace=aghor2,pace2=aghor,age = l$Age,vari = l$LifespanMedianHorvath,vcvv=ginv(C),
#       quot=T,scAge2 = aschor,scAge = aschor2,sparse=T,pgls=T)
# plreg(met = cmprrb,chron = chronrrb,chron2 = chronrrb2,bio=abiorrb,bio2=abiorrb2,
#       pace=agrrb2,pace2=agrrb,age = b[,6],vari = b[,4],vcvv = ginv(vcv(nj(gd))),
#       quot=T,scAge2 = ascrrbs,scAge = ascrrbs2,sparse=T,pgls=T)
# plreg(met = cmhor,chron = chronhor,chron2 = chronhor2,bio=abiohor,bio2=abiohor2,
#       pace=aghor2,pace2=aghor,age = l$Age,vari = l$LifespanMedianHorvath,vcvv=ginv(C),
#       quot=T,scAge2 = aschor,scAge = aschor2,sparse=F,pgls=F)
# plreg(met = cmprrb,chron = chronrrb,chron2 = chronrrb2,bio=abiorrb,bio2=abiorrb2,
#       pace=agrrb2,pace2=agrrb,age = b[,6],vari = b[,4],vcvv = ginv(vcv(nj(gd))),
#       quot=T,scAge2 = ascrrbs,scAge = ascrrbs2,sparse=F,pgls=F)
# mtext(side = 1, "Chronological Age", cex=2,outer = TRUE)
# mtext(side = 2, "Epigenetic Age (arbitrary units)", cex=2,outer = TRUE,padj = -3.5)
# mtext(side=2,line=1,at=c(0.15,0.35,0.58,0.83),text=c("RRBS","MAMMALIAN","RRBS","MAMMALIAN"), cex=2, outer=T,col=c("purple","purple","gold","gold"))
# mtext(line = -13,at=c(0.015,0.815),text=c("A","B"), cex=2, outer=T, font=2)
# mtext(line = -49,at=c(0.015,0.815),text=c("C","D"), cex=2, outer=T, font=2)
# mtext(line = -80,at=c(0.015,0.815),text=c("E","F"), cex=2, outer=T, font=2)
# mtext(line = -115,at=c(0.015,0.815),text=c("G","H"), cex=2, outer=T, font=2)
# mtext(side=2,line = 8,at=c(0.75),text=c("SPARSE"), cex=4, outer=T, font=2,col="gold")
# mtext(side=2,line = 8,at=c(0.25),text=c("NON-SPARSE"), cex=4, outer=T, font=2,col="purple")
# dev.off()
# 

##### Extra analyses not in paper 
## phylogenetic signal, PWD actual age of death, 

# #####Phylosig
# 
# #mvdnorm(unlist(x[1000,-1]),log = T,mean = 0,sigma = vcv(as.phylo(nj(as.dist(ju3)))))
# R=diag(0.005,dim(ju3)[1])
# R[upper.tri(R)]=rnorm(dim(R)[1]*(dim(R)[1]-1)/2,0,0.01)
# R=R+t(R)
# C=as.matrix(nearPD((ju3)**(1+R))[[1]])
# colnames(C)=colnames(x[,-1]);row.names(C)=colnames(x[,-1]);
# miau=nj(dist.from.cov(C))
# #maha=apply(x[,-1],1,function(x)mahalanobis(x,center=F,mean = 0,inverted=T,cov = miau))
# #plot(-log10(pchisq(maha,df = dim(miau)[1]/100,lower.tail = T)))
# maha1=apply(x[,-1],1,function(x)phylosig(tree = miau,x,test = F)[1])
# hist(maha1)
# 
# psn=pelsnp[apply(pelsnp,1,function(x){x=x[!is.na(x)];p=sum(x==0);q=sum(x);n=2*length(x);hw=log((sum(x==1)/n)/(2*p*q/n^2));abs(hw)<1 & min(p,q)>1}),]
# psn=as.dist(apply(psn,2,function(x)apply(psn,2,function(y){sum(abs(x-y),na.rm = T)/sum(!(is.na(x) | is.na(y)))})))
# #psn=vcv(nj(psn))
# psn=nj(psn)
# #miau=ginv(psn)
# miau=psn
# #maha=apply(wo[,-1:-2],1,function(x)mahalanobis(x,center=F,mean = 0,inverted=T,cov = miau))
# maha2=apply(wo[,-1:-2],1,function(x)phylosig(tree = miau,x,test=F)[[1]][1])
# hist(maha2)
# #plot(-log10(pchisq(maha,df = dim(miau)[1],lower.tail = T)))
# #hist(-pnorm(maha,mean = dim(miau)[1],sd =  2*dim(miau)[1],lower.tail = T,log.p = T))
# 
# #miau=solve(vcv(nj(gd)))
# miau=nj(gd)
# #maha=apply(a,1,function(x)mahalanobis(x,center=F,mean = 0,inverted=T,cov = miau))
# #plot(-log10(pchisq(maha,df = dim(miau)[1],lower.tail = T)))
# #hist(-pnorm(maha,mean = dim(miau)[1],sd =  2*dim(miau)[1],lower.tail = T,log.p = T))
# 
# maha3=apply(a,1,function(x)phylosig(tree = miau,x,test = F)[1])
# hist(maha3)

# #######
# ##PWD##
# #######
# ##Ages of protuguese water dogs not used
# 
# relate=read.csv("PWD_RELATIONS",header=F,sep = ":")
# #relate=read.csv("PWD_RELATIONS2",header=F,sep = ":")
# relate=relate[match(l$ExternalSampleID[match(colnames(x),l$Basename,nomatch = 0)],relate[,1],nomatch = 0),]
# relate=relate[apply(relate[,-1:-3],1,function(x)sum(x=="")<1),]
# 
# AOD=read.csv("PWD_AGES_OF_DEATH")
# AOD[,2]=as.numeric(gsub("#NUM!","NA",AOD[,2]))
# AOD=AOD[match(relate[,1],AOD[,1]),]
# l2=l[match(AOD[,1],l$ExternalSampleID),]
# l2$AOD=AOD[,2]
# x2=x[,match(l2$Basename,colnames(x),nomatch = 0)]
# l2$life=l2$Age/l2$AOD
# l2$life[l2$life>1]=1
# l3=l2[!is.na(AOD[,2]),]
# x3=x2[,match(l3$Basename,colnames(x2))]
# relate=relate[match(l3$ExternalSampleID[match(colnames(x3),l3$Basename,nomatch = 0)],relate[,1],nomatch = 0),]
# AOD2=AOD[!is.na(AOD[,2]),]
# 
# paroff=function(relate){
#   first=data.frame(gsub("\\+"," ",relate[,1]))
#   dat=relate[,-1]
#   ndat=relate[,-1]
#   i=1
#   res=data.frame()
#   while(i<=3){
#     ndat=dat[,((1:dim(dat)[2])%%(dim(dat)[2]/(2**i)))==1]
#     fu=do.call(rbind,lapply(0:(dim(ndat)[2]/2-1),function(x){miau=cbind(first[,x+1],ndat[,(2*x+1):(2*x+2)]);colnames(miau)=c("V1","V2","V3");miau}))
#     first=ndat
#     res=rbind(res,fu)
#     dat=dat[,!(((1:dim(dat)[2])%%(dim(dat)[2]/(2**i)))==1)]
#     i=i+1
#   }
#   fu=do.call(rbind,lapply(0:(dim(ndat)[2]/2-1),function(x){miau=cbind(first[,x+1],dat[,(2*x+1):(2*x+2)]);colnames(miau)=c("V1","V2","V3");miau}))
#   res=rbind(res,fu)
#   res=res[order(res[,1]),]
#   res=res[!duplicated(res[,1]),]
# }
# # miau=paroff(relate = relate[,-2])
# # miau=miau[!apply(miau,1,function(x)any(x=="")),]
# # m=makekinship(famid = rep(1,length(miau[,1])),id=miau[,1],father.id = miau[,2],mother.id = miau[,3])
# 
# rels=apply(relate[,-c(1:3)],1,function(x)apply(relate[,-c(1:3)],1,function(y)max(sum(x%in%y),sum(y%in%x))))
# colnames(rels)=gsub("\\+"," ",relate[,2]);row.names(rels)=gsub("\\+"," ",relate[,2])
# #m=m[match(colnames(rels),colnames(m),nomatch = 0),match(colnames(rels),colnames(m),nomatch = 0)]
# 
# CC=as.matrix(nearPD(rels)$mat)#C=diag(1,dim(AOD2)[1])
# #C=as.matrix(nearPD(m)$mat)#C=diag(1,dim(AOD2)[1])
# CC=solve(t(chol(CC)))
# 
# cmpwd=cmet(x3,cbind.data.frame(ag=l3$Age,Fem=as.numeric(as.factor(l3$Sex)),lsp=l3$LifespanMedianHorvath))
# abiopwd=cv.glmnet(y=CC%*%(l3$Age/AOD2[,2]),x = CC%*%t(cmpwd),nfolds = 10)
# abiopwd2=cv.glmnet(y=(l3$Age/AOD2[,2]),x = t(cmpwd),nfolds = 10)
# wbiopwd=cv.glmnet(y=CC%*%(l3$Age*AOD2[,2]),x = CC%*%t(cmpwd),nfolds = 10)
# wbiopwd2=cv.glmnet(y=(l3$Age*AOD2[,2]),x = t(cmpwd),nfolds = 10)
# chronpwd=cv.glmnet(y=CC%*%l3$Age,x = CC%*%t(cmpwd),nfolds = 10)
# chronpwd2=cv.glmnet(y=l3$Age,x = t(cmpwd),nfolds = 10)
# agpwd=pacemaker(t(CC%*%t(x3)),CC%*%(l3$Age/AOD2[,2]),co=0.5,folds = dim(x3)[2])
# agpwd2=pacemaker(x3,(l3$Age/AOD2[,2]),co=0.5,folds = dim(x3)[2])
# wgpwd=pacemaker(t(CC%*%t(x3)),CC%*%(l3$Age*AOD2[,2]),co=0.5,folds = dim(x3)[2])
# wgpwd2=pacemaker(x3,(l3$Age*AOD2[,2]),co=0.5,folds = dim(x3)[2])
# ascpwd=scAge(t(CC%*%t(x3)),CC%*%(l3$Age/AOD2[,2]),co=.5,folds = dim(x3)[2])
# ascpwd2=scAge(x3,(l3$Age/AOD2[,2]),co=.5,folds = dim(x3)[2])
# wscpwd=scAge(t(CC%*%t(x3)),CC%*%(l3$Age*AOD2[,2]),co=.5,folds = dim(x3)[2])
# wscpwd2=scAge(x3,(l3$Age*AOD2[,2]),co=.5,folds = dim(x3)[2])
# 
# appwd=t(gwas(met = x3,cova = as.matrix(nearPD(rels)$mat),data = cbind.data.frame(ag=l3$Age,Fem=as.numeric(as.factor(l3$Sex)),lsp=AOD2[,2]),quot = T))
# appwd[appwd==0]=min(appwd[,3])

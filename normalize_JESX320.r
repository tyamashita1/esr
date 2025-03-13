################################
#constants
bohr.mag <- 9.2740154E-24 #J/T
planck <- 6.6260755E-34 #J s
################################

dir="./NC3txt2/"

stdin <- commandArgs(trailingOnly = T)
lenstdin <- length(stdin)
fn <- c()
for(i in 1:lenstdin){
fn[i] = paste(dir,stdin[i],sep="")
}

colPal2 <- colorRampPalette(c("#8C36FF","#2CA7FF","#11D4AC","#BAE627","#FF8E19","#FF0E16"))

################################
# functions
# g-value 換算
gfun <- function(f,uF){ invisible( planck*uF*1E+6/f/bohr.mag*1000 ) }

# x0,y0を通る直線
linearbg <- function(s,x,x0,y0){
 c = y0-s*x0
 invisible(s*x+c) 
}

# 微分を取る
deriv1 <- function(x,y){
  len=length(x)
  ax = (x[2:(len)]+x[1:(len-1)])/2
  ay = (y[2:(len)]-y[1:(len-1)])/(x[2:(len)]-x[1:(len-1)])
  invisible(list(x=ax,y=ay))
}

# 周波数を取得してg-value換算
specfun <- function(fni,iwr){
  if(iwr==1){ print(fni,quote=F) }
  system(paste("sed -n 11,4106p ",fni," > tmp",sep=""),intern=F)
  da <- read.table("tmp",header=F,col.names=c("f","a"))
  system(paste("grep ' uF' ",fni," > uF",sep=""),intern=F)
  uf <- read.table("uF",header=F)
  write.table(uf$V3[1],file="uF2",quote=F,row.names=F,col.names=F)
  uf <- read.table("uF2",header=F,sep="F")
  freq = uf$V2[1]
  invisible(list(freq=freq,f=da$f,a=da$a,g=gfun(da$f,freq)))
}

# Mnマーカーのg-valueで横軸校正
specfun_cal <- function(fni,iwr){
  sp = specfun(fni,iwr)
  xap=seq(min(sp$g),max(sp$g),length=20000)
  apsp = approx(sp$g,sp$a,xout=xap)
  dv1 = deriv1(apsp$x,apsp$y)
# dv1 = deriv1(sp$g,sp$a)
  df = data.frame(g=dv1$x,d=abs(dv1$y))
  sf = subset(df,df$g<1.99)
  g1 = sf$g[which(sf$d==max(sf$d))]
  sf = subset(df,df$g>2.02)
  g2 = sf$g[which(sf$d==max(sf$d))]
  g1true = 1.9810
  g2true = 2.0340
  slope = (g2true-g1true)/(g2-g1)
  slope=mean(slope)
  print(slope)
  cut = g1true-slope*mean(g1)
  invisible(list(gc=slope*sp$g+cut,a=sp$a))
}

# スペクトル描画
drawser <- function(fn,iwr=1,itype=1){
  len=length(fn)
 #len=1
  cp=colPal2(len)
  for(i in 1:len){
  sp <- specfun(fn[i],iwr)
  
  if(itype==1){ lines(sp$f,sp$a,col=cp[i],lty=1,lwd=1.0) }
  if(itype==2){ lines(sp$g,sp$a,col=cp[i],lty=1,lwd=1.0) }
  #if(itype==2){ points(sp$g,sp$a,col=cp[i],pch=i%%3,cex=0.4) }
  if(itype==3){
    dv1 = deriv1(sp$g,sp$a)
    lines(dv1$x,dv1$y*1e-5,col=cp[i],lty=1,lwd=1.0)
   }
  }#i in len
}#drawser


# 積分
integ <- function(x,y,min=1.99,max=2.03){
  df = data.frame(x=x,y=y)
  sf = subset(df,df$x>min & df$x < max)
  itg = c()
  itg[1] = sf$y[1] * ( sf$x[2] - sf$x[1] )
  for(j in 2:length(sf$x)){
  itg[j] = sf$y[j] * ( sf$x[j] - sf$x[j-1] ) + itg[j-1]
  }
  invisible(list(x=sf$x,y=itg))
}


################################







#################################
#draw
pdf(paste("plot_",stdin[1],".pdf",sep=""), width = 32/2.54, height = 22/2.54)
xmax<- 342
xmin<- 331
ymax<- 2000
ymin<- -ymax

marc <- c(3.0, 4.0, 1, 1)
mgpc <- c(4, 0.5, 0)
omac <- c(1, 1, 1, 1)

split.screen(figs = c(2, 3))
par(oma=omac)
tck =  0.03
cexaxis=0.8

screen(1)
par(mar=marc)
par(mgp=mgpc)
plot(1,1,type='n',ylim=c(ymin,ymax),xlim=c(xmin,xmax),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab=NA,ylab=NA,bty="n")
#plot(1,1,type='n',ylim=c(ymin,ymax),xlim=c(xmin,xmax),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab=NA,ylab=NA,bty="n",log="y")

mtext("Magnetic field (mT)",side=1,line=1.5,cex=0.8)
mtext("Amplitude",side=2,line=3,cex=0.8)
mtext("Raw spectrum",side=3,line=0.5,cex=0.6)

if(1==1){
xv <- pretty(c(xmin,xmax))
xvs <- sprintf("%.1f",xv)
yv <- pretty(c(ymin,ymax))
yvs <- sprintf("%.0f",yv)
axis(side=1,las=F,tck=tck,cex.axis=cexaxis,at=xv,labels=xvs,family="sans")
axis(side=3,las=F,tck=tck,at=xv,labels=NA)
axis(side=2,las=T,tck=tck,cex.axis=cexaxis,at=yv,labels=yvs,family="sans")
axis(side=4,las=T,tck=tck,at=yv,labels=NA)
}

if(1==2){
ticks <- seq( 10, -10, by=-1)
tickValues <- 10**ticks
tickStrings <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(side=1,las=F,tck=tck,cex.axis=1.0,family="sans",at=tickValues,labels=tickStrings)
axis(side=2,las=T,tck=tck,cex.axis=1.0,family="sans",at=tickValues,labels=tickStrings)
t19 <- seq(1,9,by=1)
for(i in -10:10){ axis(side=1,las=F,tck=tck/2,at=t19*10**i,labels=NA) }
for(i in -10:10){ axis(side=2,las=T,tck=tck/2,at=t19*10**i,labels=NA) }
}

#legend("topright",col="black",legend=c(),cex=0.6,pch=21,pt.bg=cp,bg="white")

drawser(fn,iwr=1,itype=1)


screen(2)
xmin=1.97
xmax=2.04
par(mar=marc)
par(mgp=mgpc)
plot(1,1,type='n',ylim=c(ymin,ymax),xlim=c(xmin,xmax),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab=NA,ylab=NA,bty="n")
#plot(1,1,type='n',ylim=c(ymin,ymax),xlim=c(xmin,xmax),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab=NA,ylab=NA,bty="n",log="y")

mtext("Normalized magnetic field",side=1,line=1.5,cex=0.8)
mtext("Amplitude",side=2,line=3,cex=0.8)
mtext("X-axis converted",side=3,line=0.5,cex=0.6)

if(1==1){
xv <- pretty(c(xmin,xmax))
xvs <- sprintf("%.3f",xv)
yv <- pretty(c(ymin,ymax))
yvs <- sprintf("%.0f",yv)
axis(side=1,las=F,tck=tck,cex.axis=cexaxis,at=xv,labels=xvs,family="sans")
axis(side=3,las=F,tck=tck,at=xv,labels=NA)
axis(side=2,las=T,tck=tck,cex.axis=cexaxis,at=yv,labels=yvs,family="sans")
axis(side=4,las=T,tck=tck,at=yv,labels=NA)
}

if(1==2){
ticks <- seq( 10, -10, by=-1)
tickValues <- 10**ticks
tickStrings <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(side=1,las=F,tck=tck,cex.axis=1.0,family="sans",at=tickValues,labels=tickStrings)
axis(side=2,las=T,tck=tck,cex.axis=1.0,family="sans",at=tickValues,labels=tickStrings)
t19 <- seq(1,9,by=1)
for(i in -10:10){ axis(side=1,las=F,tck=tck/2,at=t19*10**i,labels=NA) }
for(i in -10:10){ axis(side=2,las=T,tck=tck/2,at=t19*10**i,labels=NA) }
}

#legend("topright",col="black",legend=c(),cex=0.6,pch=21,pt.bg=cp,bg="white")

drawser(fn,iwr=0,itype=2)

cp=colPal2(length(fn))
legend("topleft",col=cp,legend=fn[1:(length(fn))],lwd=1,bg="white",cex=0.7,pch=(1:(length(fn)))%%3)


screen(3)
ymax=1e+2
ymin=-ymax
xmin=1.97
xmax=2.04
par(mar=marc)
par(mgp=mgpc)
plot(1,1,type='n',ylim=c(ymin,ymax),xlim=c(xmin,xmax),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab=NA,ylab=NA,bty="n")
#plot(1,1,type='n',ylim=c(ymin,ymax),xlim=c(xmin,xmax),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab=NA,ylab=NA,bty="n",log="y")

mtext("Normalized magnetic field",side=1,line=1.5,cex=0.8)
mtext("Amplitude",side=2,line=3,cex=0.8)
mtext("Second derivative",side=3,line=0.5,cex=0.6)

if(1==1){
xv <- pretty(c(xmin,xmax))
xvs <- sprintf("%.3f",xv)
yv <- pretty(c(ymin,ymax))
yvs <- sprintf("%.1f",yv)
axis(side=1,las=F,tck=tck,cex.axis=cexaxis,at=xv,labels=xvs,family="sans")
axis(side=3,las=F,tck=tck,at=xv,labels=NA)
axis(side=2,las=T,tck=tck,cex.axis=cexaxis,at=yv,labels=yvs,family="sans")
axis(side=4,las=T,tck=tck,at=yv,labels=NA)
}

if(1==2){
ticks <- seq( 10, -10, by=-1)
tickValues <- 10**ticks
tickStrings <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(side=1,las=F,tck=tck,cex.axis=1.0,family="sans",at=tickValues,labels=tickStrings)
axis(side=2,las=T,tck=tck,cex.axis=1.0,family="sans",at=tickValues,labels=tickStrings)
t19 <- seq(1,9,by=1)
for(i in -10:10){ axis(side=1,las=F,tck=tck/2,at=t19*10**i,labels=NA) }
for(i in -10:10){ axis(side=2,las=T,tck=tck/2,at=t19*10**i,labels=NA) }
}

cp=colPal2(length(fn)-1)
#legend("topright",col=cp,legend=fn[1:(length(fn)-1)],lwd=1,bg="white",cex=0.5,pch=(1:(length(fn)-1))%%3)

drawser(fn,iwr=0,itype=3)

sp <- specfun(fn[1],iwr=0)
xap=seq(min(sp$g),max(sp$g),length=20000)
apsp = approx(sp$g,sp$a,xout=xap)
dv1 = deriv1(apsp$x,apsp$y)
df <- data.frame(g=dv1$x,d=abs(dv1$y))
sf <- subset(df,df$g<1.99)
g1 = sf$g[which(sf$d==max(sf$d))]
sf <- subset(df,df$g>2.02)
g2 = sf$g[which(sf$d==max(sf$d))]

abline(v=g1,lty=3)
abline(v=g2,lty=3)

g1true = 1.9810
g2true = 2.0340

slope <- (g2true-g1true)/(g2-g1)
slope=mean(slope)
cut <- g1true-slope*mean(g1)



screen(4)
ymax= 500
ymin=-ymax
xmin=1.97
xmax=2.04
par(mar=marc)
par(mgp=mgpc)
plot(1,1,type='n',ylim=c(ymin,ymax),xlim=c(xmin,xmax),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab=NA,ylab=NA,bty="n")
#plot(1,1,type='n',ylim=c(ymin,ymax),xlim=c(xmin,xmax),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab=NA,ylab=NA,bty="n",log="y")

mtext("Normalized magnetic field",side=1,line=1.5,cex=0.8)
mtext("Amplitude",side=2,line=3,cex=0.8)
mtext("Mn marker calibrated",side=3,line=0.5,cex=0.6)

if(1==1){
xv <- pretty(c(xmin,xmax))
xvs <- sprintf("%.3f",xv)
yv <- pretty(c(ymin,ymax))
yvs <- sprintf("%.1f",yv)
axis(side=1,las=F,tck=tck,cex.axis=cexaxis,at=xv,labels=xvs,family="sans")
axis(side=3,las=F,tck=tck,at=xv,labels=NA)
axis(side=2,las=T,tck=tck,cex.axis=cexaxis,at=yv,labels=yvs,family="sans")
axis(side=4,las=T,tck=tck,at=yv,labels=NA)
}

if(1==2){
ticks <- seq( 10, -10, by=-1)
tickValues <- 10**ticks
tickStrings <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(side=1,las=F,tck=tck,cex.axis=1.0,family="sans",at=tickValues,labels=tickStrings)
axis(side=2,las=T,tck=tck,cex.axis=1.0,family="sans",at=tickValues,labels=tickStrings)
t19 <- seq(1,9,by=1)
for(i in -10:10){ axis(side=1,las=F,tck=tck/2,at=t19*10**i,labels=NA) }
for(i in -10:10){ axis(side=2,las=T,tck=tck/2,at=t19*10**i,labels=NA) }
}

  len=length(fn)
  cp=colPal2(len)
 for(i in 1:len){
  spc <- specfun_cal(fn[i],iwr=0)
  lines(spc$gc,spc$a,col=cp[i],lty=1,lwd=1.0) 
 }#i in len




screen(5)
ymax= 0.6
ymin=-0.6
xmin=1.97
xmax=2.04
par(mar=marc)
par(mgp=mgpc)
plot(1,1,type='n',ylim=c(ymin,ymax),xlim=c(xmin,xmax),xaxt="n",yaxt="n",xaxs="i",yaxs="r",xlab=NA,ylab=NA,bty="n")
#plot(1,1,type='n',ylim=c(ymin,ymax),xlim=c(xmin,xmax),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab=NA,ylab=NA,bty="n",log="y")

mtext("Normalized magnetic field",side=1,line=1.5,cex=0.8)
mtext("Amplitude",side=2,line=3,cex=0.8)
mtext("Normalized by Mn marker amplitude",side=3,line=0.5,cex=0.6)

if(1==1){
xv <- pretty(c(xmin,xmax))
xvs <- sprintf("%.3f",xv)
yv <- pretty(c(ymin,ymax))
yvs <- sprintf("%.1f",yv)
axis(side=1,las=F,tck=tck,cex.axis=cexaxis,at=xv,labels=xvs,family="sans")
axis(side=3,las=F,tck=tck,at=xv,labels=NA)
axis(side=2,las=T,tck=tck,cex.axis=cexaxis,at=yv,labels=yvs,family="sans")
axis(side=4,las=T,tck=tck,at=yv,labels=NA)
}
box()

if(1==2){
ticks <- seq( 10, -10, by=-1)
tickValues <- 10**ticks
tickStrings <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(side=1,las=F,tck=tck,cex.axis=1.0,family="sans",at=tickValues,labels=tickStrings)
axis(side=2,las=T,tck=tck,cex.axis=1.0,family="sans",at=tickValues,labels=tickStrings)
t19 <- seq(1,9,by=1)
for(i in -10:10){ axis(side=1,las=F,tck=tck/2,at=t19*10**i,labels=NA) }
for(i in -10:10){ axis(side=2,las=T,tck=tck/2,at=t19*10**i,labels=NA) }
}

  len=length(fn)
  cp=colPal2(len)
 for(i in 1:len){
  spc <- specfun_cal(fn[i],iwr=0)
  df = data.frame(gc=spc$gc,a=spc$a)
  sf = subset(df,df$gc<1.99)
   amax1 = max(sf$a)
    gmax1 = sf$gc[which(sf$a==amax1)]
   amin1 = min(sf$a)
    gmin1 = sf$gc[which(sf$a==amin1)]
     height1 = amax1-amin1
  sf = subset(df,df$gc>2.03)
   amax2 = max(sf$a)
    gmax2 = sf$gc[which(sf$a==amax2)]
   amin2 = min(sf$a)
    gmin2 = sf$gc[which(sf$a==amin2)]
     height2 = amax2-amin2
  bgslope_up = (amax1-amax2)/(gmax1-gmax2)
  bgslope_lo = (amin1-amin2)/(gmin1-gmin2)
   bgslope=(bgslope_up+bgslope_lo)/2
   height=(height1+height2)/2
   spcbg = linearbg(bgslope,spc$gc,min(spc$gc),spc$a[which(spc$gc==min(spc$gc))])
  lines(spc$gc,spc$a/height,col=cp[i],lty=1,lwd=1.0) 
# lines(spc$gc,(spc$a - spcbg)/height,col=cp[i],lty=1,lwd=1.0) 
  lines(spc$gc,spcbg/height,col="black",lty=3,lwd=1.0) 
 }#i in len


screen(6)
ymax= 0.6
ymin=-0.6
xmin=1.97
xmax=2.04
par(mar=marc)
par(mgp=mgpc)
plot(1,1,type='n',ylim=c(ymin,ymax),xlim=c(xmin,xmax),xaxt="n",yaxt="n",xaxs="i",yaxs="r",xlab=NA,ylab=NA,bty="n")
#plot(1,1,type='n',ylim=c(ymin,ymax),xlim=c(xmin,xmax),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab=NA,ylab=NA,bty="n",log="y")

mtext("Normalized magnetic field",side=1,line=1.5,cex=0.8)
mtext("Amplitude",side=2,line=3,cex=0.8)
mtext("Height slope normalized, output to csv",side=3,line=0.5,cex=0.6)

if(1==1){
xv <- pretty(c(xmin,xmax))
xvs <- sprintf("%.3f",xv)
yv <- pretty(c(ymin,ymax))
yvs <- sprintf("%.1f",yv)
axis(side=1,las=F,tck=tck,cex.axis=cexaxis,at=xv,labels=xvs,family="sans")
axis(side=3,las=F,tck=tck,at=xv,labels=NA)
axis(side=2,las=T,tck=tck,cex.axis=cexaxis,at=yv,labels=yvs,family="sans")
axis(side=4,las=T,tck=tck,at=yv,labels=NA)
}
box()

if(1==2){
ticks <- seq( 10, -10, by=-1)
tickValues <- 10**ticks
tickStrings <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(side=1,las=F,tck=tck,cex.axis=1.0,family="sans",at=tickValues,labels=tickStrings)
axis(side=2,las=T,tck=tck,cex.axis=1.0,family="sans",at=tickValues,labels=tickStrings)
t19 <- seq(1,9,by=1)
for(i in -10:10){ axis(side=1,las=F,tck=tck/2,at=t19*10**i,labels=NA) }
for(i in -10:10){ axis(side=2,las=T,tck=tck/2,at=t19*10**i,labels=NA) }
}

  len=length(fn)
  cp=colPal2(len)
 for(i in 1:len){
  spc <- specfun_cal(fn[i],iwr=0)
  df = data.frame(gc=spc$gc,a=spc$a)
  sf = subset(df,df$gc<1.99)
   amax1 = max(sf$a)
    gmax1 = sf$gc[which(sf$a==amax1)]
   amin1 = min(sf$a)
    gmin1 = sf$gc[which(sf$a==amin1)]
     height1 = amax1-amin1
  sf = subset(df,df$gc>2.03)
   amax2 = max(sf$a)
    gmax2 = sf$gc[which(sf$a==amax2)]
   amin2 = min(sf$a)
    gmin2 = sf$gc[which(sf$a==amin2)]
     height2 = amax2-amin2
# bgslope_up = (amax1-amax2)/(gmax1-gmax2)
# bgslope_lo = (amin1-amin2)/(gmin1-gmin2)
#  bgslope=(bgslope_up+bgslope_lo)/2
#  height=(height1+height2)/2
#
   gg1=(gmax1+gmin1)/2
   gg2=(gmax2+gmin2)/2
   heightslope = (height2-height1)/(gg2-gg1)
   heightcut = height2 - heightslope*gg2
   spcheight = linearbg(heightslope,spc$gc,gg2,height2)
   print("heights")
   print(height1)
   print(height2)
   print(height )
#  spcbg = linearbg(bgslope,spc$gc,min(spc$gc),spc$a[which(spc$gc==min(spc$gc))])
  abline(h=0,lty=3)
  abline(h=0.5,lty=3)
  abline(h=-0.5,lty=3)
  lines(spc$gc,spc$a/spcheight,col=cp[i],lty=1,lwd=1.0) 
   wdf = data.frame(g=spc$gc,a=spc$a/spcheight)
  write.csv(wdf,file=paste(fn[i],"_calnr.csv",sep=""),row.names=F,quote=F)
   ap <- approx(x=wdf$g,y=wdf$a,xout=seq(1.980,2.035,length=5000))
   awdf = data.frame(g=ap$x,a=ap$y)
  write.csv(awdf,file=paste(fn[i],"_calnr_ap.csv",sep=""),row.names=F,quote=F)
 }#i in len




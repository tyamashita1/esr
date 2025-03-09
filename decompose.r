library(pracma)
library(latex2exp)

# location of normalized spectrum files
dir="data/"

#########
#constants
bohr.mag <- 9.2740154E-24 #J/T
planck <- 6.6260755E-34 #J s


#########

Yh <- function(h,hh,w){  invisible( exp(-(h-hh)**2/w**2) ) }

ah <- function(h1,h2,h3,h){ invisible( sqrt((h3**2-h2**2)*(h**2-h1**2)) ) }
bh <- function(h1,h2,h3,h){ invisible( sqrt((h2**2-h1**2)*(h3**2-h**2)) ) }

Xa <- function(h1,h2,h3,h,hh,w){ 
  m = (ah(h1,h2,h3,hh)/bh(h1,h2,h3,hh))**2
  r = ellipke(m)$k * Yh(h,hh,w) / bh(h1,h2,h3,hh) / hh**2 * h1*h2*h3*2/pi 
  invisible(r)
}
Xb <- function(h1,h2,h3,h,hh,w){ 
  m = (bh(h1,h2,h3,hh)/ah(h1,h2,h3,hh))**2
  r = ellipke(m)$k * Yh(h,hh,w) / ah(h1,h2,h3,hh) / hh**2 * h1*h2*h3*2/pi 
  invisible(r)
}


daikei_sekibun <- function(x,y){
  len = length(x)
  integ = ( y[1:(len-1)]+y[2:(len)] ) * ( x[2:(len)]-x[1:(len-1)] ) /2.0
  invisible(sum(integ))
}

deriv1 <- function(x,y){
  len=length(x)
  ax = (x[2:(len)]+x[1:(len-1)])/2
  ay = (y[2:(len)]-y[1:(len-1)])/(x[2:(len)]-x[1:(len-1)])
  invisible(list(x=ax,y=ay))
}

deriv3 <- function(x,y){
  len = length(x)
  xdrv = (x[1:(len-2)]+x[3:(len)])/2
  drv = ( -y[1:(len-2)]+y[3:(len)] ) / (-x[1:(len-2)]+x[3:(len)])
  invisible(list(x=xdrv,y=drv))
}

spcpowder <- function(hplt,h1,h2,h3,w){
  ep=1e-6
  n=5000
  aplt = hplt*0

  if(h1!=h2&h2!=h3){
    hva = seq(h1+ep,h2-ep,length=n)
    hvb = seq(h2+ep,h3-ep,length=n)
  }else if(h1==h2&h2!=h3){
    hva = c()
    hvb = seq(h2+ep,h3-ep,length=n)
  }else if(h1!=h2&h2==h3){
    hva = seq(h1+ep,h2-ep,length=n)
    hvb = c()
  }
  
  for(i in 1:length(hplt)){
  h=hplt[i]
   suma=0.0
   sumb=0.0
   if(length(hva)>0){ yva=Xa(h1,h2,h3,h,hva,w);suma=daikei_sekibun(hva,yva) }
   if(length(hvb)>0){ yvb=Xb(h1,h2,h3,h,hvb,w);sumb=daikei_sekibun(hvb,yvb) }
   aplt[i]=suma+sumb
  }#i nn
  
   invisible(list(hplt=hplt,aplt=aplt))
}#spcpowder

spcpowderderiv <- function(hplt,h1,h2,h3,w){
  if(h1<=h2&h2<=h3){ 
  spc0 = spcpowder(hplt,h1,h2,h3,w)
  dv = deriv3(spc0$hplt,spc0$aplt)
  invisible(list(x=dv$x,y=dv$y))
  }else{print("h1<=h2<=h3 required");stop()}
}

#########

derivSimple <- function(x,y){
  l=length(x)
  dy = ( y[2:(l)] - y[1:(l-1)] )/(x[2:(l)]-x[1:(l-1)])
  dx = (x[2:(l)]+x[1:(l-1)])/2
  invisible(list(x=dx,y=dy))
}

derivA <- function(x,y){
  n=5
  l=length(x)
  mn=n+1
  mx=l-n
  for(i in 1:n){
   if(i==1){ dy = ( y[(mn+i):(mx+i)] - y[(mn-i):(mx-i)] )/( x[(mn+i):(mx+i)] - x[(mn-i):(mx-i)] ) }
   if(i>=2){ dy = dy + ( y[(mn+i):(mx+i)] - y[(mn-i):(mx-i)] )/( x[(mn+i):(mx+i)] - x[(mn-i):(mx-i)] ) }
  }
  dy=dy/n
  dx = x[mn:mx]
  invisible(list(x=dx,y=dy))
}


Lor <- function(c,w,x){
 invisible( -1/pi * 2*w*( x-c ) / ( (x-c)**2+w**2 )**2  )
}

Gau <- function(c,w,x){
 invisible( -2/w**3/sqrt(pi) * ( x-c ) * exp(-(x-c)**2/w**2)  )
}

Gau6 <- function(c,w,x){
 invisible( -2/w/sqrt(pi) * ( x-c ) * exp(-(x-c)**2/w**2) )
}

GauLor <- function(c,w,x){
 y = -(2*exp(-(c-x)**2)*(x-c)*(2+c**2+w**2-2*c*x+x**2))/((c**2+w**2-2*c*x+x**2)*1e+5)**3
 invisible(y)
}

MnBGlinear <- function(x){
   gMn1=1.9810
   gMn2=2.0340
   wMn1=1.548e-04
   wMn2=1.537e-04
   slope=-3.643e-03
   cut=-6.779e-03
   aMn1=-5.710e-08
   aMn2=-5.743e-08
   invisible( aMn1 * Lor(1/gMn1, wMn1, x) + aMn2 * Lor(1/gMn2, wMn2, x) + slope * (x - 0.5) + cut )
}

MnMarker <- function(x,ww){
   gMn1=1.9810
   gMn2=2.0340
   wMn1=1.548e-04
   wMn2=1.537e-04
   aMn1=-5.710e-08
   aMn2=-5.743e-08
   invisible( aMn1 * Lor(1/gMn1, wMn1*ww, x) + aMn2 * Lor(1/gMn2, wMn2*ww, x) )
}

BGlinear <- function(x,slope=-3.643e-03,cut=0.0,shift=0.49165){
   invisible( slope * (x - shift) + cut )
}

Native <- function(x){
#aNat=c(0.07499440,0.12643628,0.00674372)
#wNat=c(0.00086581,0.00250166,0.00049291)
#ginvNat=c(0.49872650,0.49827217,0.49856562)
 aNat=c(0.08375975,0.24960170,0.01183407)
 wNat=c(0.00103601,0.00260525,0.00050936)
 ginvNat=c(0.49879748,0.49840057,0.49857366)
 model=x*0
 for(j in 1:length(aNat)){
   model = model + aNat[j]*Gau6(ginvNat[j],wNat[j],x)
 }
 invisible( model )
}


NativeW <- function(x,wscale){
#aNat=c(0.07499440,0.12643628,0.00674372)
#wNat=c(0.00086581,0.00250166,0.00049291)
#ginvNat=c(0.49872650,0.49827217,0.49856562)
 aNat=c(0.08375975,0.24960170,0.01183407)
 wNat=c(0.00103601,0.00260525,0.00050936)
 ginvNat=c(0.49879748,0.49840057,0.49857366)
 model=x*0
 for(j in 1:length(aNat)){
   model = model + aNat[j]*Gau6(ginvNat[j],wNat[j]*wscale,x)
 }
 invisible( model )
}



NativeFit <- function(x,aNat,wNat,ginvNat){
 model=x*0
 for(j in 1:length(aNat)){
   model = model + aNat[j]*GauLor(ginvNat[j],wNat[j],x)
 }
 invisible( model )
}



pmrunif <- function(len){ 
  r = runif(min=-1,max=+1,n=len)
  if(r>=0){ rr=1 }
  if(r<0){ rr=-1 }
  invisible(rr)
}


eprunif <- function(ep,len,filt){ 
  r = runif(min=-ep,max=+ep,n=len)
  invisible(r*filt+1.0)
}

weightfilt <- function(x,Hw0,Hw1,WeightIn,WeightOUT){
  y=seq(0,0,length=length(x))
  for(i in 1:length(x)){
   if(x[i]>Hw0 & x[i]<Hw1){ y[i] = WeightIn }else{ y[i] = WeightOUT}
  }
  invisible(y)
}
weightdiff <- function(vdiff,x,Hw0,Hw1,WeightIn,WeightOUT){
  for(i in 1:length(vdiff)){
   if(x[i]>Hw0 & x[i]<Hw1){ vdiff[i] = vdiff[i]*WeightIn }else{vdiff[i] = vdiff[i]*WeightOUT}
  }
  invisible(vdiff)
}
weightdiffZero <- function(vdiff,x,Hw0,Hw1){
  for(i in 1:length(vdiff)){
   if(x[i]>Hw0 & x[i]<Hw1){ vdiff[i] = vdiff[i] }else{vdiff[i]=0}
  }
  invisible(vdiff)
}

weightIgnoreNearZero <- function(x,y,xv,yv,yigzero){
  w=rep(1,length(x))
 if(yigzero>0){
  ap = approx(xv,yv,xout=x,yleft=0,yright=0)
  for(i in 1:length(x)){
   if(abs(ap$y[i])<yigzero){ w[i] = 0 }
  }
 }#if
  invisible(w)
}

weightfiltlower <- function(y,masklow){
  w=rep(1,length(y))
  for(i in 1:length(y)){
   if(y[i]<masklow){ w[i]=0 }
  }
  invisible(w)
}

calcslope <- function(x,y){
  l=length(x)
  ta=data.frame(x,y)
  sa=subset(ta,ta$x>=0.490 & ta$x<=0.494)
  sb=subset(ta,ta$x>=0.504 & ta$x<=0.506)
   y0 = max(sa$y)
   x0 = sa$x[which(sa$y==y0)]
   y1 = max(sb$y)
   x1 = sb$x[which(sb$y==y1)]
   slopeMAX = (y1-y0)/(x1-x0)
   y0 = min(sa$y)
   x0 = sa$x[which(sa$y==y0)]
   y1 = min(sb$y)
   x1 = sb$x[which(sb$y==y1)]
   slopeMIN = (y1-y0)/(x1-x0)
  invisible((slopeMAX+slopeMIN)/2.0)
}

calccutslope <- function(x,y,gmn1,gmx1,gmn2,gmx2){
  l=length(x)
  ta=data.frame(x,y)
  sa=subset(ta,ta$x>=gmn1 & ta$x<=gmx1)
  sb=subset(ta,ta$x>=gmn2 & ta$x<=gmx2)
  x1 = mean(sa$x)
  x2 = mean(sb$x)
  y1 = mean(sa$y)
  y2 = mean(sb$y)
  slope = (y2-y1)/(x2-x1)
  cut = y2-slope*x2
  invisible(list(slope=slope,cut=cut))
}

#smooth <- function(x,y,bin){
#  if(bin>=1){
#    xs = c()
#    ys = c()
#    ii=0
#   for(i in (bin+1):(length(x)-bin)){
#     ii=ii+1
#     xs[ii] = mean(x[(i-bin):(i+bin)])
#     ys[ii] = mean(y[(i-bin):(i+bin)])
#   }#for
#   return(list(x=xs,y=ys))
#  }else{
#   return(list(x=x,y=y))
#  }
#}

smooth <- function(x, y, bin) {
  if (bin >= 1) {
    # Use the filter function from base R to compute the rolling mean
    w <- rep(1, 2 * bin + 1) / (2 * bin + 1)
    xs <- stats::filter(x, w, sides = 2)
    ys <- stats::filter(y, w, sides = 2)
    
    # Remove NA values from the filtered result
    valid_idx <- !is.na(xs) & !is.na(ys)
    return(list(x = xs[valid_idx], y = ys[valid_idx]))
  } else {
    return(list(x = x, y = y))
  }
}

#########


stdin <- commandArgs(trailingOnly = T)
lenstdin <- length(stdin)
fn <- c()
rawfn <- c()
rawfn = stdin
fn = stdin
fn[1] = paste(dir,stdin[1],sep="")

spcfile=fn[1]
parfile=fn[2]
parfileOUT=fn[3]
sinkout=paste("sink_",rawfn[1],".log",sep="")

da <- read.csv(spcfile,header=F,skip=1,col.names=c("gc","amp"))
df <- data.frame(ginv=1/da$gc,amp=da$amp)
df <- df[order(df$ginv),]
pdfname = paste(rawfn[1],".pdf",sep="")
filename=rawfn[1]

par <- read.table(parfile,header=F,skip=1,sep=":")

filtA=c(1,1,1,1,1,1,1,1)
filtW=c(1,1,1,1,1,1,1,1)
filtG=c(1,1,1,1,1,1,1,1)
filtSubNatA=c(1,1,1,1,1,1)
filtSubNatW=c(1,1,1,1,1,1)
filtSubNatG=c(1,1,1,1,1,1)
asubNat=c(0,0,0,0,0,0)
wsubNat=c(0.01,0.01,0.01,0.01,0.01,0.01)
gsubNat=c(0.5,0.5,0.5,0.5,0.5,0.5)
smoothingbin=0

for(i in 1:length(par$V1)){
 if(par$V1[i]!="NOTE"){ 
   if(par$V2[i]=="check initial guess"){ initial_guess=par$V3[i] }
   if(par$V2[i]=="perturb initial guess"){ initial_perturb=par$V3[i] }
   if(par$V2[i]=="ep perturb"){ ep_perturb=par$V3[i] }
   if(par$V2[i]=="smoothing(D=0)"){ smoothingbin=par$V3[i] }

   if(par$V2[i]=="number of hot batch"){ nbatch_hot=par$V3[i] }
   if(par$V2[i]=="number of cooling batch"){ nbatch_cooling=par$V3[i] }
   if(par$V2[i]=="number of cold batch"){ nbatch_cold=par$V3[i] }
   if(par$V2[i]=="Temperature hot"){ Thot=par$V3[i] }
   if(par$V2[i]=="Temperature cold"){ Tcold=par$V3[i] }
   if(par$V2[i]=="number of iteration"){ nite=par$V3[i] }
  
   if(par$V2[i]=="g1 carbonate"){ g1=par$V3[i] }
   if(par$V2[i]=="g2 carbonate"){ g2=par$V3[i] }
   if(par$V2[i]=="g3 carbonate"){ g3=par$V3[i] }
   if(par$V2[i]=="width carbonate"){ wCO2neg=par$V3[i]; filtW[3]=0 }
   if(par$V2[i]=="amplitude carbonate"){ aCO2neg=par$V3[i]; if(par$V1[i]=="LOCK"){filtA[3]=0} }
  
   if(par$V2[i]=="g1 O2-"){ g1O=par$V3[i] }
   if(par$V2[i]=="g2 O2-"){ g2O=par$V3[i] }
   if(par$V2[i]=="g3 O2-"){ g3O=par$V3[i] }
   if(par$V2[i]=="width O2-"){ wO2neg=par$V3[i]; filtW[5]=0 }
   if(par$V2[i]=="amplitude O2-"){ aO2neg=par$V3[i]; if(par$V1[i]=="LOCK"){filtA[5]=0} }
  
   if(par$V2[i]=="g1 Hole center"){ g1Hole=par$V3[i] }
   if(par$V2[i]=="g2 Hole center"){ g2Hole=par$V3[i] }
   if(par$V2[i]=="g3 Hole center"){ g3Hole=par$V3[i] }
   if(par$V2[i]=="width Hole center"){ wHole=par$V3[i]; filtW[6]=0 }
   if(par$V2[i]=="amplitude Hole center"){ aHole=par$V3[i]; if(par$V1[i]=="LOCK"){filtA[6]=0} }
  
   if(par$V2[i]=="g1 add1"){ g1add1=par$V3[i] }
   if(par$V2[i]=="g2 add1"){ g2add1=par$V3[i] }
   if(par$V2[i]=="g3 add1"){ g3add1=par$V3[i] }
   if(par$V2[i]=="width add1"){ wadd1=par$V3[i]; filtW[7]=0 }
   if(par$V2[i]=="amplitude add1"){ aadd1=par$V3[i]; if(par$V1[i]=="LOCK"){filtA[7]=0} }
  
#  if(par$V2[i]=="slope(linear BG)"){ slope=par$V3[i]; if(par$V1[i]=="LOCK"){filtA[1]=0} }
   if(par$V2[i]=="cut(linear BG)"){ cut0=par$V3[i];cutfilt=1; if(par$V1[i]=="LOCK"){cutfilt=0} }
   if(par$V2[i]=="cut auto 1=auto 0=tune"){ cutauto=par$V3[i];cutfilt=1; if(cutauto==1){cutfilt=0} }
   if(par$V2[i]=="LinearBG region min1"){ ginvbgmn1=par$V3[i] }
   if(par$V2[i]=="LinearBG region max1"){ ginvbgmx1=par$V3[i] }
   if(par$V2[i]=="LinearBG region min2"){ ginvbgmn2=par$V3[i] }
   if(par$V2[i]=="LinearBG region max2"){ ginvbgmx2=par$V3[i] }
#   if(par$V2[i]=="cut sign ibatch"){ ibatchcut=par$V3[i]; }

   if(par$V2[i]=="amplitude Mn"){ aMn=par$V3[i]; if(par$V1[i]=="LOCK"){filtA[1]=0} }
   if(par$V2[i]=="width scale Mn"){ wwMn=par$V3[i]; if(par$V1[i]=="LOCK"){filtW[1]=0} }

   if(par$V2[i]=="amplitude Native"){ aNat=par$V3[i]; if(par$V1[i]=="LOCK"){filtA[2]=0} }
   if(par$V2[i]=="broadening Native"){ wscaleNative=par$V3[i]; if(par$V1[i]=="LOCK"){filtW[2]=0} }

   if(par$V2[i]=="subamp1 Native"){ asubNat[1]=par$V3[i]; if(par$V1[i]=="LOCK"){filtSubNatA[1]=0} }
   if(par$V2[i]=="subamp2 Native"){ asubNat[2]=par$V3[i]; if(par$V1[i]=="LOCK"){filtSubNatA[2]=0} }
   if(par$V2[i]=="subamp3 Native"){ asubNat[3]=par$V3[i]; if(par$V1[i]=="LOCK"){filtSubNatA[3]=0} }
   if(par$V2[i]=="subamp4 Native"){ asubNat[4]=par$V3[i]; if(par$V1[i]=="LOCK"){filtSubNatA[4]=0} }
   if(par$V2[i]=="subamp5 Native"){ asubNat[5]=par$V3[i]; if(par$V1[i]=="LOCK"){filtSubNatA[5]=0} }
   if(par$V2[i]=="subamp6 Native"){ asubNat[6]=par$V3[i]; if(par$V1[i]=="LOCK"){filtSubNatA[6]=0} }
   if(par$V2[i]=="subwid1 Native"){ wsubNat[1]=par$V3[i]; if(par$V1[i]=="LOCK"){filtSubNatW[1]=0}  }
   if(par$V2[i]=="subwid2 Native"){ wsubNat[2]=par$V3[i]; if(par$V1[i]=="LOCK"){filtSubNatW[2]=0}  }
   if(par$V2[i]=="subwid3 Native"){ wsubNat[3]=par$V3[i]; if(par$V1[i]=="LOCK"){filtSubNatW[3]=0}  }
   if(par$V2[i]=="subwid4 Native"){ wsubNat[4]=par$V3[i]; if(par$V1[i]=="LOCK"){filtSubNatW[4]=0}  }
   if(par$V2[i]=="subwid5 Native"){ wsubNat[5]=par$V3[i]; if(par$V1[i]=="LOCK"){filtSubNatW[5]=0}  }
   if(par$V2[i]=="subwid6 Native"){ wsubNat[6]=par$V3[i]; if(par$V1[i]=="LOCK"){filtSubNatW[6]=0}  }
   if(par$V2[i]=="subgiv1 Native"){ gsubNat[1]=par$V3[i]; if(par$V1[i]=="LOCK"){filtSubNatG[1]=0}  }
   if(par$V2[i]=="subgiv2 Native"){ gsubNat[2]=par$V3[i]; if(par$V1[i]=="LOCK"){filtSubNatG[2]=0}  }
   if(par$V2[i]=="subgiv3 Native"){ gsubNat[3]=par$V3[i]; if(par$V1[i]=="LOCK"){filtSubNatG[3]=0}  }
   if(par$V2[i]=="subgiv4 Native"){ gsubNat[4]=par$V3[i]; if(par$V1[i]=="LOCK"){filtSubNatG[4]=0}  }
   if(par$V2[i]=="subgiv5 Native"){ gsubNat[5]=par$V3[i]; if(par$V1[i]=="LOCK"){filtSubNatG[5]=0}  }
   if(par$V2[i]=="subgiv6 Native"){ gsubNat[6]=par$V3[i]; if(par$V1[i]=="LOCK"){filtSubNatG[6]=0}  }
  
   if(par$V2[i]=="amplitude BG1"){ aBG1=par$V3[i]; if(par$V1[i]=="LOCK"){filtA[4]=0} }
   if(par$V2[i]=="g BG1"){ gBG1=par$V3[i]; if(par$V1[i]=="LOCK"){filtG[4]=0} }
   if(par$V2[i]=="width BG1"){ wBG1=par$V3[i]; if(par$V1[i]=="LOCK"){filtW[4]=0} }
  
   if(par$V2[i]=="amplitude BG2"){ aBG2=par$V3[i]; if(par$V1[i]=="LOCK"){filtA[8]=0} }
   if(par$V2[i]=="g BG2"){ gBG2=par$V3[i]; if(par$V1[i]=="LOCK"){filtG[8]=0} }
   if(par$V2[i]=="width BG2"){ wBG2=par$V3[i]; if(par$V1[i]=="LOCK"){filtW[8]=0} }
  
   if(par$V2[i]=="refine ite"){ ite_refine=par$V3[i] }
   if(par$V2[i]=="initial ep"){ ep0=par$V3[i] }
  
   if(par$V2[i]=="ymax(gross)"){ ymax1=par$V3[i] }
   if(par$V2[i]=="ymax(close-up)"){ ymax2=par$V3[i] }

   if(par$V2[i]=="Fitting Weight region min"){ Hw0=par$V3[i] }
   if(par$V2[i]=="Fitting weight region max"){ Hw1=par$V3[i] }
   if(par$V2[i]=="Weight inside fitting region"){ WeightIn=par$V3[i] }
   if(par$V2[i]=="Weight outside fitting region"){ WeightOUT=par$V3[i] }
   if(par$V2[i]=="2nd deriv weight region min"){ H2w0=par$V3[i] }
   if(par$V2[i]=="2nd deriv weight region max"){ H2w1=par$V3[i] }
   if(par$V2[i]=="2nd deriv weight"){ Weight2=par$V3[i] }
   if(par$V2[i]=="IgnoreNearZero 2nd deriv"){ yigzero=par$V3[i] }
   if(par$V2[i]=="mask lower amp"){ masklow=par$V3[i] }

   if(par$V2[i]=="SINK"){ isink=par$V3[i] }
 }#NOTE
}


if(cutauto==0){ 
  slopeMn=calcslope(df$ginv,df$amp)
  shiftBGlin = 0.490
}else if(cutauto==1){ 
  shiftBGlin = 0.0
  cutslope = calccutslope(df$ginv,df$amp,ginvbgmn1,ginvbgmx1,ginvbgmn2,ginvbgmx2)
  cut0 = cutslope$cut
  slopeMn = cutslope$slope
  cutfilt=0
}

 df0=df
 smoothspec = smooth(df0$ginv,df0$amp,smoothingbin)
 df = data.frame(ginv=smoothspec$x,amp=smoothspec$y)
 print(head(df))


xraw = df0$ginv
xsmt = df$ginv

h1 = 1/g1
h2 = 1/g2
h3 = 1/g3
hplt = seq(min(xraw),max(xraw),length= 500)
spc = spcpowderderiv(hplt,h1,h2,h3,wCO2neg)
apf_CO2neg <- approxfun(spc$x,spc$y*1e-4,yleft=0.0,yright=0.0)
apf_CO2negx <- apf_CO2neg(xraw)

spc = spcpowderderiv(hplt,1/g1O,1/g2O,1/g3O,wO2neg)
apf_O2neg <- approxfun(spc$x,spc$y*1e-4,yleft=0.0,yright=0.0)
apf_O2negx <- apf_O2neg(xraw)

spc = spcpowderderiv(hplt,1/g1Hole,1/g2Hole,1/g3Hole,wHole)
apf_Hole <- approxfun(spc$x,spc$y*1e-4,yleft=0.0,yright=0.0)
apf_Holex <- apf_Hole(xraw)

spc = spcpowderderiv(hplt,1/g1add1,1/g2add1,1/g3add1,wadd1)
apf_add1 <- approxfun(spc$x,spc$y*1e-4,yleft=0.0,yright=0.0)
apf_add1x <- apf_add1(xraw)


xdv = derivA(df$ginv,df$amp)$x
ydv = derivA(df$ginv,df$amp)$y
ddv = data.frame(xdv=xdv,ydv=ydv)
dvmax = max(abs(ydv))
weightfilt1 <- weightfilt(xsmt,Hw0,Hw1,WeightIn,WeightOUT)
weightfilt1low <- weightfiltlower(df$amp,masklow)
weightfilt1dvzero <- weightIgnoreNearZero(df$ginv,df$amp,xdv,ydv/dvmax,yigzero)
weightfilt2 <- weightfilt(xdv,H2w0,H2w1,Weight2,0.0)

 qdiff=1e+30
 ep=ep0

 fitAMP=c(aMn,aNat,aCO2neg,aBG1,aO2neg,aHole,aadd1,aBG2)
 fitWID=c(wwMn,wscaleNative,wCO2neg,wBG1,wO2neg,wHole,wadd1,wBG2)
 fitGIV=c(0,0,0,1/gBG1,0,0,0,1/gBG2)


icnv=0
cnvite = c()
cnvsig = c()
cnvdif = c()

if(isink==1){ sink(sinkout) }

fitAMP1 = fitAMP
fitWID1 = fitWID
fitGIV1 = fitGIV
asubNat1 = asubNat
wsubNat1 = wsubNat
gsubNat1 = gsubNat
cut1=cut0
if(initial_perturb==1){
  fitAMP1 = fitAMP * eprunif(ep_perturb,length(fitAMP),filtA)
  fitWID1 = fitWID * eprunif(ep_perturb,length(fitWID),filtW)
  fitGIV1 = fitGIV * eprunif(ep_perturb/100,length(fitGIV),filtG)
  asubNat1 = asubNat * eprunif(ep_perturb,length(asubNat),filtSubNatA)
  wsubNat1 = wsubNat * eprunif(ep_perturb,length(wsubNat),filtSubNatW)
  gsubNat1 = gsubNat * eprunif(ep_perturb/100,length(gsubNat),filtSubNatG)
  cut1 = cut0*eprunif(ep_perturb,1,c(cutfilt))
}

nbatch = nbatch_hot + nbatch_cooling + nbatch_cold
nb12 = nbatch_hot + nbatch_cooling
nb1 = nbatch_hot
Temp = numeric(nbatch)
Temp[1] = Tcold
Temp[2:nb1] = Thot
Temp[(nb12+1):nbatch] = Tcold
Trate = (Tcold/Thot)**(1.0/(nbatch_cooling-1))
for(i in 1:nbatch_cooling){ Temp[nb1+i]=Thot*Trate**(i-1) }
#Temp[(nb1+1):nb12] = seq(Thot,Tcold,length=nbatch_cooling)
print(Temp)

#npm<-seq(1,1,length=nbatch)
#for(i in 1:ibatchcut){npm[i]=1}
#for(i in (ibatchcut+1):nbatch){npm[i]=2}

if(initial_guess==0){

filtAsave=filtA
filtA=filtA*0

ifail=0
for(ibatch in 1:nbatch){
 Temp0 = Temp[ibatch]
 if(ibatch>nbatch_hot){filtA=filtAsave}
for(i in 1:nite){

model = BGlinear(xraw,slope=slopeMn,cut=cut1,shift=shiftBGlin) + MnMarker(xraw,fitWID1[1])*fitAMP1[1] + 
        NativeFit(xraw,asubNat1,wsubNat1,gsubNat1) + apf_CO2negx * fitAMP1[3] + 
        Gau6(fitGIV1[4],fitWID1[4],xraw) * fitAMP1[4] + apf_O2negx * fitAMP1[5] + apf_Holex * fitAMP1[6] +
        Gau6(fitGIV1[8],fitWID1[8],xraw) * fitAMP1[8] + apf_add1x * fitAMP1[7] 
modelsmth = smooth(xraw,model,smoothingbin)
vdiff1 = ( df$amp - modelsmth$y )**2
vdiff2 = ( ydv - derivA(modelsmth$x,modelsmth$y)$y )**2
diff1 = sum( vdiff1*weightfilt1*weightfilt1dvzero*weightfilt1low ) 
diff2 = sum( vdiff2*weightfilt2 ) /dvmax
diff = diff1+diff2

judge=exp(-diff/Temp0)
judgeRmd=runif(1,min=0,max=1)

if( diff<qdiff ){
 qdiff=diff
 fitAMP2 = fitAMP1; fitWID2 = fitWID1; fitGIV2 = fitGIV1; asubNat2=asubNat1;wsubNat2=wsubNat1;gsubNat2=gsubNat1; cut2=cut1; 
 fitAMP  = fitAMP1; fitWID  = fitWID1; fitGIV  = fitGIV1; asubNat =asubNat1;wsubNat =wsubNat1;gsubNat =gsubNat1; cut0=cut1;
 print(sprintf("ib %d ite %d ep %e T %e fitAMP[3] %f qdiff %.8f *****",ibatch,i,ep,Temp0,fitAMP[3],qdiff),quote=F)
 icnv=icnv+1;cnvsig[icnv]=fitAMP2[3];cnvdif[icnv]=diff;
 ifail=0
}else if(judge>judgeRmd){
 fitAMP = fitAMP1; fitWID = fitWID1; fitGIV = fitGIV1; asubNat=asubNat1;wsubNat=wsubNat1;gsubNat=gsubNat1; cut0=cut1; 
 print(sprintf("ib %d ite %d ep %e T %e fitAMP[3] %f  diff %.8f",ibatch,i,ep,Temp0,fitAMP[3],diff),quote=F)
#icnv=icnv+1;cnvsig[icnv]=fitAMP[3];cnvdif[icnv]=diff;
#ifail=0
}else{
 ifail=ifail+1
 if(ifail==ite_refine){
  if(ep>1e-9){ep=ep/2.0}
  ifail=0
  fitAMP  = fitAMP2; fitWID  = fitWID2; fitGIV  = fitGIV2; asubNat =asubNat2;wsubNat =wsubNat2;gsubNat =gsubNat2; cut0=cut2;
 }
}#diff if

# next guess
  fitAMP1 = fitAMP * eprunif(ep,length(fitAMP),filtA)
  fitWID1 = fitWID * eprunif(ep,length(fitWID),filtW)
  fitGIV1 = fitGIV * eprunif(ep/100,length(fitGIV),filtG)
  asubNat1 = asubNat * eprunif(ep,length(asubNat),filtSubNatA)
  wsubNat1 = wsubNat * eprunif(ep,length(wsubNat),filtSubNatW)
  gsubNat1 = gsubNat * eprunif(ep/100,length(gsubNat),filtSubNatG)
  cut1 = cut0*eprunif(ep/10,1,c(cutfilt))


}#ite
 ep=ep0/10
 fitAMP  = fitAMP2; fitWID  = fitWID2; fitGIV  = fitGIV2; asubNat =asubNat2;wsubNat =wsubNat2;gsubNat =gsubNat2; cut0=cut2;
}#ibatch


fitAMP = fitAMP2
fitWID = fitWID2
fitGIV = fitGIV2
asubNat = asubNat2
wsubNat = wsubNat2
gsubNat = gsubNat2
cut0=cut2

}else{cnvsig=c(1);cnvdif=c(1)}#initial_guess

modelBGL = BGlinear(xraw,slope=slopeMn,cut=cut0,shift=shiftBGlin) 
modelMn = MnMarker(xraw,fitWID[1])*fitAMP[1]
#modelNat = NativeW(xraw,fitWID[4])*fitAMP[4] 
modelNat = NativeFit(xraw,asubNat,wsubNat,gsubNat)
modelCO2neg = apf_CO2neg(xraw) * fitAMP[3] 
modelBG1 = Gau6(fitGIV[4],fitWID[4],xraw) * fitAMP[4] 
modelBG2 = Gau6(fitGIV[8],fitWID[8],xraw) * fitAMP[8] 
modelO2neg = apf_O2neg(xraw) * fitAMP[5] 
modelHole = apf_Hole(xraw) * fitAMP[6] 
modeladd1 = apf_add1(xraw) * fitAMP[7] 

for(i in 1:length(fitAMP)){ print(sprintf("fitAMP[%d] = %.8f",i,fitAMP[i])) }
for(i in 1:length(fitWID)){ print(sprintf("fitWID[%d] = %.8f",i,fitWID[i])) }
for(i in 1:length(fitGIV)){ print(sprintf("fitGIV[%d] = %.8f",i,fitGIV[i])) }
for(i in 1:length(asubNat)){ print(sprintf("asubNat[%d] = %.8f",i,asubNat[i])) }
for(i in 1:length(wsubNat)){ print(sprintf("wsubNat[%d] = %.8f",i,wsubNat[i])) }
for(i in 1:length(gsubNat)){ print(sprintf("gsubNat[%d] = %.8f",i,gsubNat[i])) }

# print(sprintf("slope = %.8f",fitAMP[1]),quote=F)
# print(sprintf("cut =%.8f",fitAMP[2]),quote=F)
# print(sprintf("aMn =%.8f",fitAMP[3]),quote=F)
# print(sprintf("aNat =%.8f",fitAMP[4]),quote=F)
# print(sprintf("aCO2neg =%.8f",fitAMP[5]),quote=F)
# print(sprintf("aBG1 =%.8f",fitAMP[6]),quote=F)


wpar=par
for(i in 1:length(wpar$V1)){
 if(wpar$V2[i]=="width carbonate"){ wpar$V3[i]=fitWID[3] }
 if(wpar$V2[i]=="amplitude carbonate"){ wpar$V3[i]=fitAMP[3] }
 if(wpar$V2[i]=="slope(linear BG)"){ wpar$V3[i]=slopeMn }
 if(wpar$V2[i]=="cut(linear BG)"){ wpar$V3[i]=cut0 }
 if(wpar$V2[i]=="amplitude Mn"){ wpar$V3[i]=fitAMP[1] }
 if(wpar$V2[i]=="width scale Mn"){ wpar$V3[i]=fitWID[1] }
 if(wpar$V2[i]=="amplitude Native"){ wpar$V3[i]=fitAMP[2] }
 if(wpar$V2[i]=="broadening Native"){ wpar$V3[i]=fitWID[2] }
 if(wpar$V2[i]=="amplitude BG1"){ wpar$V3[i]=fitAMP[4] }
 if(wpar$V2[i]=="g BG1"){ wpar$V3[i]=1/fitGIV[4] }
 if(wpar$V2[i]=="width BG1"){ wpar$V3[i]=fitWID[4] }
 if(wpar$V2[i]=="width O2-"){ wpar$V3[i]=fitWID[5] }
 if(wpar$V2[i]=="amplitude O2-"){ wpar$V3[i]=fitAMP[5] }
 if(wpar$V2[i]=="width Hole"){ wpar$V3[i]=fitWID[6] }
 if(wpar$V2[i]=="amplitude Hole"){ wpar$V3[i]=fitAMP[6] }
 if(wpar$V2[i]=="subamp1 Native"){ wpar$V3[i]=asubNat[1] }
 if(wpar$V2[i]=="subamp2 Native"){ wpar$V3[i]=asubNat[2] }
 if(wpar$V2[i]=="subamp3 Native"){ wpar$V3[i]=asubNat[3] }
 if(wpar$V2[i]=="subamp4 Native"){ wpar$V3[i]=asubNat[4] }
 if(wpar$V2[i]=="subamp5 Native"){ wpar$V3[i]=asubNat[5] }
 if(wpar$V2[i]=="subamp6 Native"){ wpar$V3[i]=asubNat[6] }
 if(wpar$V2[i]=="subwid1 Native"){ wpar$V3[i]=wsubNat[1] }
 if(wpar$V2[i]=="subwid2 Native"){ wpar$V3[i]=wsubNat[2] }
 if(wpar$V2[i]=="subwid3 Native"){ wpar$V3[i]=wsubNat[3] }
 if(wpar$V2[i]=="subwid4 Native"){ wpar$V3[i]=wsubNat[4] }
 if(wpar$V2[i]=="subwid5 Native"){ wpar$V3[i]=wsubNat[5] }
 if(wpar$V2[i]=="subwid6 Native"){ wpar$V3[i]=wsubNat[6] }
 if(wpar$V2[i]=="subgiv1 Native"){ wpar$V3[i]=gsubNat[1] }
 if(wpar$V2[i]=="subgiv2 Native"){ wpar$V3[i]=gsubNat[2] }
 if(wpar$V2[i]=="subgiv3 Native"){ wpar$V3[i]=gsubNat[3] }
 if(wpar$V2[i]=="subgiv4 Native"){ wpar$V3[i]=gsubNat[4] }
 if(wpar$V2[i]=="subgiv5 Native"){ wpar$V3[i]=gsubNat[5] }
 if(wpar$V2[i]=="subgiv6 Native"){ wpar$V3[i]=gsubNat[6] }
}
write.table(wpar,parfileOUT,sep=":",quote=F,row.names=F)

dfcnv = data.frame(NumUpdate=c(1:length(cnvsig)),cnvsig=cnvsig,cnvdif=cnvdif)
write.csv(dfcnv,paste(parfileOUT,".cnv.csv",sep=""),quote=F,row.names=F)


##############

pdf(pdfname, width = 28/2.54, height = 22/2.54)
xmax <- 0.507
xmin <- 0.489

marc <- c(3, 4, 1, 1) 
mgpc <- c(4, 0.5, 0) 
omac <- c(1, 1, 1, 1) 

split.screen(figs = c(2, 3)) 
par(oma=omac)
tck =  0.03
cexaxis = 0.8 

for(iscr in 1:5){
if(iscr==1){ ymax <-  ymax1; ymin <- -ymax; xmin <- 0.489; xmax <- 0.507; ylog=0 }
if(iscr==2){ ymax <-  ymax2; ymin <- -ymax; xmin <- 0.489; xmax <- 0.507; ylog=0 }
if(iscr==3){ ymax <-  5e-2; ymin <- -ymax; xmin <- 0.489; xmax <- 0.507; ylog=0 }
#if(iscr==4){ ymax <-  ymax1; ymin <- -ymax; xmin <- 1.975; xmax <- 2.04 }
#if(iscr==5){ ymax <-  ymax2; ymin <- -ymax; xmin <- 1.975; xmax <- 2.04 }
if(iscr==4){ ymax <-  max(cnvsig); ymin <- min(cnvsig); xmin <- 0; xmax <- length(cnvsig); ylog=0 }
#if(iscr==5){ ymax <-  max(cnvdif); ymin <- min(cnvdif); xmin <- 0; xmax <- length(cnvdif); ylog=1 }
if(iscr==5){ ymax <-  1e+1; ymin <- 1e-4; xmin <- 0; xmax <- length(cnvdif); ylog=1 }
screen(iscr) 
par(mar=marc)
par(mgp=mgpc)
if(ylog==0){plot(1,1,type='n',ylim=c(ymin,ymax),xlim=c(xmin,xmax),xaxt="n",yaxt="n",xaxs="r",yaxs="r",xlab=NA,ylab=NA,bty="n") }
if(ylog==1){plot(1,1,type='n',ylim=c(ymin,ymax),xlim=c(xmin,xmax),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab=NA,ylab=NA,bty="n",log="y") }
box()

if(iscr==2){
abline(v=Hw0,lty=2,lwd=0.5)
abline(v=Hw1,lty=2,lwd=0.5)
abline(v=H2w0,lty=3,lwd=0.5)
abline(v=H2w1,lty=3,lwd=0.5)
abline(h=masklow,lty=2,lwd=0.5)
}
if(iscr==3){
abline(v=Hw0,lty=2,lwd=0.5)
abline(v=Hw1,lty=2,lwd=0.5)
abline(v=H2w0,lty=3,lwd=0.5)
abline(v=H2w1,lty=3,lwd=0.5)
abline(h=masklow,lty=2,lwd=0.5)
}

if(iscr<=3){ mtext(TeX("$\\beta H/h\\nu$"),side=1,line=1.5,cex=0.8) }
if(iscr>=4){ mtext("Number of updates",side=1,line=1.5,cex=0.8) }
if(iscr<=3){ mtext("Intensity (arb.u.)",side=2,line=3,cex=0.8) }
if(iscr==4){ mtext("Amplitude of radical (arb.u.)",side=2,line=3,cex=0.8) }
if(iscr==5){ mtext("Difference (arb.u.)",side=2,line=3,cex=0.8) }

xv <- pretty(c(xmin,xmax))
if(iscr<=3){ xvs <- sprintf("%.3f",xv) }else{ xvs <- sprintf("%.0f",xv) }
axis(side=1,las=F,tck=tck,cex.axis=cexaxis,at=xv,labels=xvs,family="sans")
axis(side=3,las=F,tck=tck,at=xv,labels=NA)

#yv <- pretty(c(ymin,ymax))
#if(iscr<=5){ yvs <- sprintf("%.2f",yv) }else{ yvs <- sprintf("%.3f",yv) }
#axis(side=2,las=T,tck=tck,cex.axis=cexaxis,at=yv,labels=yvs,family="sans")
if(ylog==0){
axis(side=2,las=T,tck=tck,cex.axis=cexaxis,family="sans")
axis(side=4,las=T,tck=tck,labels=NA)
}else{
ticks <- seq( 10, -10, by=-1);tickValues <- 10**ticks;t19 <- seq(1,9,by=1)
tickStrings <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(side=2,las=T,tck=tck,cex.axis=1.0,family="sans",at=tickValues,labels=tickStrings)
for(i in -10:10){ axis(side=2,las=T,tck=tck/2,at=t19*10**i,labels=NA) }
}

colPal2 <- colorRampPalette(c("#8C36FF","#2CA7FF","#11D4AC","#BAE627","#FF8E19"))
cp <- colPal2(9) 

#### plot
if(iscr==1){
  ap <- approx(df$ginv,df$amp,xout=xsmt)
  ap0 <- approx(df0$ginv,df0$amp,xout=xraw)
  
  modelTOT=modelBGL+modelMn+modelNat+modelCO2neg+modelBG1+modelBG2+modelO2neg+modelHole+modeladd1
  modelTOTs = smooth(xraw,modelTOT,smoothingbin)$y
  vdiff1 = ( df$amp - modelTOTs )**2
  vdiff2 = ( ydv - derivA(xsmt,modelTOTs)$y )**2
  diff1 = sum( vdiff1*weightfilt1 ) 
  diff2 = sum( vdiff2*weightfilt2 ) /dvmax
  diff = diff1+diff2
  dfdiff = data.frame(diff1=diff1,diff2=diff2,diff=diff)
  write.csv(dfdiff,paste(parfileOUT,".diff.csv",sep=""),quote=F,row.names=F)
  
  lg=c(
  "BGL",
  "Mn",
  "Nat",
  "CO2-",
  "BG1",
  "BG2",
  "O2-",
  "Hole",
  "add1",
  "TOT"
  )
  
# NatResidual=ap$y - modelBGL - modelMn - modelCO2neg - modelBG1 - modelBG2 - modelO2neg - modelHole - modeladd1

  wdf<-data.frame(ginv=xsmt,g=1/xsmt,measure=ap$y, 
                  modelBGL=smooth(xraw,modelBGL,smoothingbin)$y, 
                  modelMn=smooth(xraw,modelMn,smoothingbin)$y,
                  modelNat=smooth(xraw,modelNat,smoothingbin)$y,
                  modelCO2neg=smooth(xraw,modelCO2neg,smoothingbin)$y,
                  modelBG1=smooth(xraw,modelBG1,smoothingbin)$y,
                  modelBG2=smooth(xraw,modelBG2,smoothingbin)$y,
                  modelO2neg=smooth(xraw,modelO2neg,smoothingbin)$y,
                  modelHole=smooth(xraw,modelHole,smoothingbin)$y,
                  modeladd1=smooth(xraw,modeladd1,smoothingbin)$y,
                  modelTOT=smooth(xraw,modelTOT,smoothingbin)$y 
                 )
  spdfname=gsub("\\.pdf$", "", pdfname)
  write.csv(wdf,file=paste(spdfname,"_opt.csv",sep=""),quote=F,row.names=F)
}#iscr==1, write

if(iscr<=2){

  points(df$ginv,df$amp,cex=0.2,col="#66666600",pch=21,bg="#666666")
#  points(df$ginv[weightfilt1dvzero==0],df$amp[weightfilt1dvzero==0],cex=0.1,pch=21,col="#FFFFFF00",bg="#0000FF")
#  points(df$ginv[weightfilt1dvzero==1],df$amp[weightfilt1dvzero==1],cex=0.1,pch=21,col="#FFFFFF00",bg="#000000")
# if(sum(abs(modelBGL))>0.0){    lines(xraw,modelBGL,col=cp[1]) }
# if(sum(abs(modelMn))>0.0){     lines(xraw,modelMn,col=cp[2]) }
# if(sum(abs(modelNat))>0.0){    lines(xraw,modelNat,col=cp[3]) }
# if(sum(abs(modelCO2neg))>0.0){ lines(xraw,modelCO2neg,col=cp[4]) }
# if(sum(abs(modelBG1))>0.0){    lines(xraw,modelBG1,col=cp[5]) }
# if(sum(abs(modelBG2))>0.0){    lines(xraw,modelBG2,col=cp[6]) }
# if(sum(abs(modelO2neg))>0.0){  lines(xraw,modelO2neg,col=cp[7]) }
# if(sum(abs(modelHole))>0.0){   lines(xraw,modelHole,col=cp[8]) }
# if(sum(abs(modeladd1))>0.0){   lines(xraw,modeladd1,col=cp[9]) }
# lines(xraw,modelTOT,col="red")
  if(sum(abs(modelBGL))>0.0){    lines(wdf$ginv,wdf$modelBGL,col=cp[1]) }
  if(sum(abs(modelMn))>0.0){     lines(wdf$ginv,wdf$modelMn,col=cp[2]) }
  if(sum(abs(modelNat))>0.0){    lines(wdf$ginv,wdf$modelNat,col=cp[3]) }
  if(sum(abs(modelCO2neg))>0.0){ lines(wdf$ginv,wdf$modelCO2neg,col=cp[4]) }
  if(sum(abs(modelBG1))>0.0){    lines(wdf$ginv,wdf$modelBG1,col=cp[5]) }
  if(sum(abs(modelBG2))>0.0){    lines(wdf$ginv,wdf$modelBG2,col=cp[6]) }
  if(sum(abs(modelO2neg))>0.0){  lines(wdf$ginv,wdf$modelO2neg,col=cp[7]) }
  if(sum(abs(modelHole))>0.0){   lines(wdf$ginv,wdf$modelHole,col=cp[8]) }
  if(sum(abs(modeladd1))>0.0){   lines(wdf$ginv,wdf$modeladd1,col=cp[9]) }
  lines(wdf$ginv,wdf$modelTOT,col="red")
  legend((xmin+xmax)/2,ymax,col=c(cp,"red"),legend=lg,cex=0.5,ncol=3,bg="white",lwd=1,xjust=0.5)
  #legend("topright",col=c(cp,"red"),legend=lg,cex=0.5,ncol=3,bg="white",lwd=1)
  
}else if(iscr==3){
  dv = derivA(df$ginv, df$amp)
  points(dv$x,dv$y/dvmax,cex=0.2,col="#66666600",pch=21,bg="#666666")
# points(dv$x[abs(dv$y/dvmax)<yigzero],dv$y[abs(dv$y/dvmax)<yigzero]/dvmax,cex=0.1,col="#0000FF")
  dv = derivA(xraw, modelTOT)
  lines(dv$x,dv$y/dvmax,cex=0.1,col="#ff0000")
}else if(iscr==4){
  points(1:length(cnvsig),cnvsig,cex=0.3,col="blue")
  lines(1:length(cnvsig),cnvsig,col="blue")
}else if(iscr==5){
  points(1:length(cnvdif),cnvdif,cex=0.3,col="blue")
  lines(1:length(cnvdif),cnvdif,col="blue")
}

}#iscr

if(isink==1){ 
  warnings()
  sink() 
}






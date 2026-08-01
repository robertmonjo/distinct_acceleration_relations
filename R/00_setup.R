# Monjo & Banik (2025), ApJ 992, 35 — reproducible pipeline
# Sourced in order by run_all.R (shared session state).
# Stage 0: packages, constants, models, data loading, galaxy/cluster preprocessing.

cex0 <- 1.5


#library("rgdal")
library("jpeg")
library("png")
library("fields")

Flecha = function(xmin,xmax,ymin,ymax, mode, dx=0.1,dy=0.1,px=0.5,py=0.5, color="lightblue",border="blue", new=TRUE, ...)
{
  Dx = xmax - xmin  
  Dy = ymax - ymin
  
  if(new)plot.new()
  
  if(mode=="right")
    polygon(c(xmin,xmin+(px-dx)*Dx,xmin+px*Dx,xmax, xmin+px*Dx, xmin+(px-dx)*Dx, xmin),
            c(ymin+dy*Dy,ymin+dy*Dy, ymin, (ymax+ymin)/2,  ymax, ymax-dy*Dy, ymax-dy*Dy), col=color, border=border,...)
  
  if(mode=="left")
    polygon(c(xmax,     xmax-(px-dx)*Dx,  xmax-px*Dx,   xmin,     xmax-px*Dx, xmax-(px-dx)*Dx, xmax),
            c(ymax-dy*Dy,ymax-dy*Dy, ymax,     (ymax+ymin)/2,  ymin,    ymin+dy*Dy, ymin+dy*Dy), col=color, border=border,...)
  
  if(mode=="up")
    polygon(c(xmin+dx*Dx,xmin+dx*Dx,     xmin,       (xmax+xmin)/2,  xmax,      xmax-dx*Dx, xmax-dx*Dx),
            c(ymin,     ymin+(py-dy)*Dy, ymin+py*Dy,     ymax,      ymin+py*Dy,  ymin+(py-dy)*Dy, ymin ), col=color, border=border,...)
  
  if(mode=="down")
    polygon(c(xmin+dx*Dx,xmin+dx*Dx,     xmin,       (xmax+xmin)/2,  xmax,      xmax-dx*Dx, xmax-dx*Dx),
            c(ymax,      ymax-(py-dy)*Dy, ymax-py*Dy,     ymin,      ymax-py*Dy,  ymax-(py-dy)*Dy, ymax ), col=color, border=border,...)
  
}


acos. = function(pi_gg)
{
  ggg = (1/pi_gg)*(pi/2)
  gamma1 = atan(ggg)
  gamma1[ggg>0 & !is.na(ggg)] = ggg[ggg>0& !is.na(ggg)]^0.68/1.34
  gamma2 = gamma1
  gamma3 = gamma1
  gamma4 = gamma1
  gamma2[ggg>1 & !is.na(ggg)] = acos((gamma1/ggg)[ggg>1 & !is.na(ggg)])
  gamma3[ggg>1 & !is.na(ggg)] = acos((gamma2/ggg)[ggg>1 & !is.na(ggg)])
  gamma4[ggg>1 & !is.na(ggg)] = acos(cos(gamma3[ggg>1 & !is.na(ggg)])/0.91)
  return(gamma4)
}

polygon_ = function(x,v,y1,y2,GRID=F, cols="gray", dens=NULL, add=T)
{
  years = x; v1 = v;
  xI  = sort(x, decreasing=T); 
  y2I = y2; for(i in 1:length(y2)) y2I[i] = y2[length(y2)-i+1]; 
  
  lg  = is.na(y1)==F; nye = length(x[lg])
  lgI = is.na(y2I)==F;
  
  xx = c(x[lg],  xI[lgI])
  yy = c(y1[lg], y2I[lgI])
  
  if(add==F)
    plot(xx,yy*NA,ylim=c(min(yy,na.rm=T),max(yy,na.rm=T)))
  polygon(xx, yy, col=cols, border=cols, anomaly=dens)
  if(GRID)grid(lwd=cex0*c(1,1,2.7,2.7)[f])
}


#http://astroweb.cwru.edu/SPARC/
#data/MassModels_Lelli2016c.mrt.txt
#load("data/figure_Milgrom4.RData")

#### Unidades ####
c0 = 3*10^8
Msol = 	 1.9891e30
ua = 1.496e+11
kpc = 3261.8116478174 * 365*24*3600*c0
pc = kpc/1000
kms = 1000
T0 = 13.7*10^9*365*24*3600 
GN =  6.674e-11
g00 = 2*pi
UZ = 4222.735 #Megaparsec
Mpc = 0.000237 
lp = 1.616199e-35 #metros
rk = 200*kpc
#vk = sqrt(ak0*rk)

E0 = 2.5e4*Msol^0.5
D0 = 6.25e3*Msol^0.5*kpc^-1
20*GN*D0^2 #GM/r^2

xx =1:25
plot(xx,(((GN*Msol*8*10^10/(xx*kpc))^2 + 2*(GN*Msol*8*10^10)*c0/T0/6)^0.25-(GN*8*Msol*10^10/(xx*kpc))^0.5)/kms)


gc = 0.5*pi
gu = pi/3
nn = 10^{1:20}
mm = 10^{1:20}*Msol
r0 = 5*kpc

mm/(4*pi*r0^3)
d0 = Msol/(0.1*pc)^3
gm = asin(sqrt(sin(gu)^2 + (sin(gc)^2-sin(gu)^2)*1/(1 + 1/(4/3*pi*2*GN*T0^2*d0))))
a0 = 2*c0/T0*cos(gm)/gm

gm = asin(sqrt(sin(gu)^2 + (sin(gc)^2-sin(gu)^2)*2*GN*mm/r0/(r0^2/T0^2 + 2*GN*mm/r0)))
gm = asin(sqrt(sin(gu)^2 + (sin(gc)^2-sin(gu)^2)*1/(1 + 1/(4/3*pi*2*GN*T0^2*mm/(4/3*pi*r0^3)))))
a0 = 2*c0/T0*cos(gm)/gm
a_0 = 1.2e-10


data = read.table("data/RAR_data_Indranil.txt")
log10_gN_values = data[,1]
log10_g_values = data[,2]
Dispersion = data[,3]
N_data = data[,4]
sigma_dex_values = Dispersion/sqrt(N_data)


nu_simple = 0.5 + sqrt(0.25 + a_0/10^log10_gN_values)
nu_standard = sqrt(0.5 + sqrt(0.25 + a_0*a_0/100^log10_gN_values))
nu_sharp = sqrt(a_0/10^log10_gN_values)

for(i in 1:length(nu_sharp))
  if(nu_sharp[i] < 1.0)
  nu_sharp[i] = 1.0

nu_sum_rule_values = 1 + sqrt(a_0/10^log10_gN_values)
nu_sum_rule_half_a0_values = 1 + sqrt(0.5*a_0/10^log10_gN_values)
nu_sum_rule_quarter_a0_values = 1 + sqrt(0.25*a_0/10^log10_gN_values)
log10_nu_sum_rule_values = log10(nu_sum_rule_values)
log10_nu_sum_rule_half_a0_values = log10(nu_sum_rule_half_a0_values)
log10_nu_sum_rule_quarter_a0_values = log10(nu_sum_rule_quarter_a0_values)
#plot(log10_gN_values, log10_gN_values + log10_nu_sum_rule_values)
nu_MLS_values = 1./(1 - exp(-sqrt((10^log10_gN_values)/a_0)))
log10_nu_MLS_values = log10(nu_MLS_values)
log10_nu_simple_values = log10(nu_simple)
log10_nu_standard_values = log10(nu_standard)
log10_nu_sharp_values = log10(nu_sharp)
log10_nu_obs_values = log10_g_values - log10_gN_values

#plot(log10_gN_values, log10_gN_values + log10_nu_MLS_values)
#plt.xlim(xmin = log10_gN_values.min() - 0.01, xmax = log10_gN_values.max() + 0.01)
#plt.xlabel(r'$\textrm{log}_{10}$ $g_N$ [m/s$^2$]')
plot(log10_gN_values, log10_nu_MLS_values - log10_nu_obs_values, col='blue', ylab= 'MLS')
lines(log10_gN_values, log10_nu_simple_values - log10_nu_obs_values, yerr = sigma_dex_values, col='black')
segments(lwd=4,x0=log10_gN_values, y0=log10_nu_simple_values - log10_nu_obs_values-sigma_dex_values, log10_gN_values,  log10_nu_simple_values - log10_nu_obs_values + sigma_dex_values, col='black')## 'Simple')
lines(log10_gN_values, log10_nu_standard_values - log10_nu_obs_values, yerr = sigma_dex_values, col='red')
segments(lwd=4,x0=log10_gN_values, y0=log10_nu_standard_values - log10_nu_obs_values-sigma_dex_values, log10_gN_values,  log10_nu_standard_values - log10_nu_obs_values + sigma_dex_values, col='red')## 'Simple')
lines(log10_gN_values, log10_nu_sharp_values - log10_nu_obs_values, yerr = sigma_dex_values, col='green')
segments(lwd=4,x0=log10_gN_values, y0=log10_nu_sharp_values - log10_nu_obs_values-sigma_dex_values, log10_gN_values,  log10_nu_standard_values - log10_nu_obs_values + sigma_dex_values, col='green')## 'Simple')



###### 0. SPARC: Mass Discrepancy-Acceleration Relation ######


#http://astroweb.cwru.edu/SPARC/
#http://astroweb.cwru.edu/SPARC/SPARC_Lelli2016c.mrt
sparc = read.table("data/MassModels_Lelli2016c.mrt.txt",skip=25)
sparc.h = read.fwf("data/MassModels_Lelli2016c.mrt.txt",skip=9, nrows=10, widths = nchar(c("54- 59 ","F6.2","   ","km/s","         ","Vbul   ","Bulge velocity contribution (3)")))
colnames(sparc) = sparc.h$V6
spa = aggregate(sparc[,-1],by=list(sparc$`  ID   `),mean)

ids = unique(sparc$`  ID   `)





# smcGa07.Vk/smcGa07.V es defineixen a la seccio seguent; mk00 es recupera com a alias de mk000.

###  1. McGaugh et al. 2007 #####
######  1.1.0. Mass models used in McGaugh et al. 2007, ApJ, 659, 149 ###### 
#http://astroweb.case.edu/ssm/data/
#https://iopscience.iop.org/article/10.1086/511807
mcGa07 = read.table("data/S_McGaugh2007_mass_discrepancy.txt",sep="\t",header=T)
smcGa07 = read.table("data/McGaugh2007.txt",sep="\t",header=T)
mcGa07 = mcGa07[!is.na(smcGa07$R),]
namess = smcGa07$Name
namesm = mcGa07$Name
namesg = table(namess)
library("stringr")
namesG = namesg
names(namesG) = str_to_upper(names(namesg))
names(namesG)[c(1,2,3)] = c("F563-1","F563-V2","F568-V1" )

smcGa07.Vk = sqrt(smcGa07$Vst^2 + smcGa07$Vgas^2)
smcGa07.V = (smcGa07.Vk^2*kms^2 + 2*smcGa07.Vk^2*kms^2*smcGa07$R*kpc*c0/(9*T0))^0.25/kms
ak000 = smcGa07.Vk^2/(smcGa07$R*kpc)*kms^2
aobs0 = smcGa07$Vobs^2/(smcGa07$R*kpc)*kms^2
di000 = smcGa07$Vobs^2/smcGa07.Vk^2
mk000 = smcGa07.Vk^2*(smcGa07$R*kpc)*kms^2
mk00 = mk000
# seccio 1.1.1 (despres de Fig02); Fig02 l'usa, per tant cal calcular-lo abans.
ak00 = smcGa07.Vk^2/(smcGa07$R*kpc)*kms^2
di00 = smcGa07$Vobs^2/smcGa07.Vk^2
gg = 1/ak00*(2*c0/T0)/(di00^2 - 1 )
gg[gg>30] = NA
gg[gg<1] = NA
ak00[(ak00/((smcGa07$R*kpc)/T0^2))^0.25 >30] = NA
gg = round(gg,2)
gamma1 = acos.((pi/2)/gg)
gamma1[gamma1<1.4] = NA
gg000 = 1/ak000*(2*c0/T0)/(di000^2 - 1 )
gal_vevH = (sqrt(2)*smcGa07.Vk*kms*T0)/(smcGa07$R*kpc)

smcGa07.R = aggregate(smcGa07$R,by=list(namess),max,na.rm=T)
smcGa07.VN = aggregate(smcGa07.Vk,by=list(namess),max,na.rm=T)
smcGa07.VO = aggregate(smcGa07.V,by=list(namess),mean,na.rm=T)
smcGa07.mass = smcGa07.VN$x^2*smcGa07.R$x*kpc*kms^2/GN


plot(gal_vevH, 1/gg000, log="xy",ylim=c(0.01,0.5),xlim=c(10,1000),pch=20)
### plot(c(rar.escape_newton^2/rar.radius, ak000)/(1.2*10^-10),c(rar.escape_newton/rar.escape_hubble, gal_vevH),log="xy")
lmak = lm(gal_vevH~0+c(ak000/(1.2e-10)))
     
##### Por cada galaxia por separado #####

seq_gammagc = seq(0.44,0.5,0.001)*pi
#seq_gammagc = seq(0.2,0.33,0.01)*pi
seq_galeps = unique(c(1,seq(1,10,0.1),seq(10,200,1)))

gal_xi2_0.47 = array(NA, dim=c(length(names(namesg)), length(seq_galeps)),
                dimnames = list(names(namesg), seq_galeps))
gal_xi2_0.48 = array(NA, dim=c(length(names(namesg)), length(seq_galeps)),
                     dimnames = list(names(namesg), seq_galeps))
gal_xi2_0.5 = array(NA, dim=c(length(names(namesg)), length(seq_galeps)),
                    dimnames = list(names(namesg), seq_galeps))
gal_xi22 = array(NA, dim=c(length(names(namesg)), length(seq_galeps), length(seq_gammagc)),
                dimnames = list(names(namesg), seq_galeps, seq_gammagc))

#### Fig000 calculamos las curvas de las galaxias ######

smcGa07.Rg = rep(NA,length(smcGa07$R))
for(i in 1:length(namesg))
{
  cat(names(namesg)[i],"\n")
  lg_gal = namess==names(namesg)[i]
  smcGa07.Rg[lg_gal] = max(smcGa07$R[lg_gal])
  
  smcGa07.Vk = sqrt(smcGa07$Vst^2 + smcGa07$Vgas^2)
  smcGa07.V = (smcGa07.Vk^2*kms^2 + 2*smcGa07.Vk^2*kms^2*smcGa07$R*kpc*c0/(9*T0))^0.25/kms
  
  ## plot(smcGa07.Vk[lg_gal]^2/(smcGa07$R[lg_gal]*kpc)*kms^2, smcGa07$Vobs[lg_gal]^2/smcGa07.Vk[lg_gal]^2,log="xy")
  ## Acceleration of Keppler
  ak000 = smcGa07.Vk[lg_gal]^2/(smcGa07$R[lg_gal]*kpc)*kms^2
  aobs0 = smcGa07$Vobs[lg_gal]^2/(smcGa07$R[lg_gal]*kpc)*kms^2
  di000 = smcGa07$Vobs[lg_gal]^2/smcGa07.Vk[lg_gal]^2
  mk000 = smcGa07.Vk[lg_gal]^2*(smcGa07$R[lg_gal]*kpc)*kms^2
  gal_vevH = (sqrt(2)*smcGa07.Vk[lg_gal]*kms*T0)/(smcGa07$R[lg_gal]*kpc)
  gammaU = pi/3
  
  
  gg000_btfr = mean(1/ak000*(2*c0/T0)/(di000^2),na.rm=T)
  ## ggg00 = 1/(aobs0-ak000)*(c0/T0)
  gg000 = 1/ak000*(2*c0/T0)/(di000^2 - 1 )
  
  #gg000[gg000>30] = NA
  #gg000[gg000<1] = NA
  #ak000[(ak000/((smcGa07$R[lg_gal]*kpc)/T0^2))^0.25 >30] = NA
  
  gg000 = round(gg000,2)
  g_sys1 = acos.((pi/2)/gg000)
  g_sys1[g_sys1< 1] = NA
  

  
  for(k in 1:length(seq_galeps))
  {
    gammagc = 0.47*pi
    gaempty = pi/3
    gal_eps = seq_galeps[k]
    #quotient_vevh = (1)/(gal_eps^2*gal_vevH^(-2) + 1)
    #quotient_vevh = (gal_vevH-gal_eps)^2/(gal_eps^2 + gal_vevH^2)
    quotient_vevh = abs(gal_vevH^2-gal_eps^2)/(gal_eps^2 + gal_vevH^2)
    g_sys0_pred = asin((sin(gaempty)^2 +  (sin(gammagc)^2-sin(gaempty)^2)*quotient_vevh)^0.5)
    gg000_pred = g_sys0_pred/cos(g_sys0_pred)
    gal_xi2_0.48[i,k] = mean(abs(log(gg000_pred)-log(gg000)),na.rm=T)
    
    gammagc = 0.48*pi
    gaempty = pi/3
    gal_eps = seq_galeps[k]
    #quotient_vevh = (1)/(gal_eps^2*gal_vevH^(-2) + 1)
    #quotient_vevh = (gal_vevH-gal_eps)^2/(gal_eps^2 + gal_vevH^2)
    quotient_vevh = abs(gal_vevH^2-gal_eps^2)/(gal_eps^2 + gal_vevH^2)
    g_sys0_pred = asin((sin(gaempty)^2 +  (sin(gammagc)^2-sin(gaempty)^2)*quotient_vevh)^0.5)
    gg000_pred = g_sys0_pred/cos(g_sys0_pred)
    gal_xi2_0.48[i,k] = mean(abs(log(gg000_pred)-log(gg000)),na.rm=T)
    
    gammagc = 0.50*pi
    gaempty = pi/3
    gal_eps = seq_galeps[k]
    #quotient_vevh = (1)/(gal_eps^2*gal_vevH^(-2) + 1)
    #quotient_vevh = (gal_vevH-gal_eps)^2/(gal_eps^2 + gal_vevH^2)
    quotient_vevh = abs(gal_vevH^2-gal_eps^2)/(gal_eps^2 + gal_vevH^2)
    g_sys0_pred = asin((sin(gaempty)^2 +  (sin(gammagc)^2-sin(gaempty)^2)*quotient_vevh)^0.5)
    gg000_pred = g_sys0_pred/cos(g_sys0_pred)
    gal_xi2_0.5[i,k] = mean(abs(log(gg000_pred)-log(gg000)),na.rm=T)
  }
  
  for(k in 1:length(seq_galeps))
    for(j in 1:length(seq_gammagc))
    {
    gammagc = seq_gammagc[j]  
    gaempty = pi/3 ###seq_gammagc[j] ##pi/3
    gal_eps = seq_galeps[k]
    quotient_vevh = abs(gal_vevH^2-gal_eps^2)/(gal_eps^2 + gal_vevH^2)
    g_sys1_pred = asin((sin(gaempty)^2 +  (sin(gammagc)^2-sin(gaempty)^2)*quotient_vevh)^0.5)
    gg001_pred = g_sys1_pred/cos(g_sys1_pred)
    gal_xi22[i,k,j] = mean(abs(log(gg001_pred)-log(gg000)),na.rm=T)
  }
}  
round(pchisq(apply(gal_xi22,1,min,na.rm=T)[namesg>1],namesg[namesg>1]-1),2)
round(pchisq(apply( gal_xi2_0.48,1,min,na.rm=T)[namesg>1],namesg[namesg>1]-1),2)
round(pchisq(apply( gal_xi2_0.5,1,min,na.rm=T)[namesg>1],namesg[namesg>1]-1),2)

galeps_0.47 = rep(NA,length(namesg))
galeps_0.48 = rep(NA,length(namesg))
galeps_0.5 = rep(NA,length(namesg))
galeps_k = rep(NA,length(namesg))
gammagc_k = rep(NA,length(namesg))
 
#### Fig00A_galaxies.pdf (anomalies) ####

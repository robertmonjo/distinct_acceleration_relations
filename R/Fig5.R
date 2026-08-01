# Monjo & Banik (2025), ApJ 992, 35 — reproducible pipeline
# Sourced in order by run_all.R (shared session state).
# Figure C3 + cluster tables. -> outputs/Fig5.pdf, outputs/table_cluster_RAR_*.txt



for(iclu in 1:10)
{
  
  eps_clusA[iclu,1] = seq_eps[order(clus_xi2A[iclu,])[1]]
  eps_clusA[iclu,2:3] = range(seq_eps[clus_xi2A[iclu,] < p95[iclu]])
  
  eps_clusB[iclu,1] = seq_eps[order(clus_xi2B[iclu,])[1]]
  eps_clusB[iclu,2:3] = range(seq_eps[clus_xi2B[iclu,] < p95[iclu]])
  
  eps_clusb[iclu,1] = seq_eps[order(clus_xi2b[iclu,])[1]]
  eps_clusb[iclu,2:3] = range(seq_eps[clus_xi2b[iclu,] < p95[iclu]])
  
  mejor_galC = order(apply(clus_xi2C[iclu,,],2,min))[1]
  mejor_epsC = order(apply(clus_xi2C[iclu,,],1,min))[1]
  
  eps_clusC[iclu,1] = seq_eps[order(clus_xi2C[iclu,,mejor_galC])[1]]
  eps_clusC[iclu,2:3] = range(seq_eps[clus_xi2C[iclu,,mejor_galC] < p67[iclu]])
  gal_clusC[iclu,1] = seq_gal[order(clus_xi2C[iclu,mejor_epsC,])[1]]
  gal_clusC[iclu,2:3] = range(seq_gal[clus_xi2C[iclu,mejor_epsC,] < p67[iclu]])
  
  
  clus_z[iclu] = max(rar$z[rar$Name==clusters[iclu]])
  clus_mass[iclu] = max(rar.mass[rar$Name==clusters[iclu]])
  clus_massmin[iclu] = min(rar.mass[rar$Name==clusters[iclu]])
  clus_rad[iclu] = max(rar.radius[rar$Name==clusters[iclu]])
  clus_radmin[iclu] = min(rar.radius[rar$Name==clusters[iclu]])
  clus_ve[iclu,1] = mean(rar.escape_newton[rar$Name==clusters[iclu]])
  clus_ve[iclu,2:3] = range(rar.escape_newton[rar$Name==clusters[iclu]])
  clus_vH[iclu,1] = mean(rar.escape_hubble[rar$Name==clusters[iclu]])
  clus_vH[iclu,2:3] = range(rar.escape_hubble[rar$Name==clusters[iclu]])
  clus_vevH2[iclu,1] = mean((rar.escape_newton[rar$Name==clusters[iclu]]/rar.escape_hubble[rar$Name==clusters[iclu]])^2)
  clus_vevH2[iclu,2:3] = range((rar.escape_newton[rar$Name==clusters[iclu]]/rar.escape_hubble[rar$Name==clusters[iclu]])^2)
  
}

### RAR Table 01 ####

table_clusB = cbind(eps_clusB,
                    round(pchisq(apply(clus_xi2B[, ],1,min),rar.num-1),2))
colnames(table_clusB) = c("eps_med","eps_low","eps_upp","p_value")
write.table(table_clusB, "outputs/Table1_general_model.txt", quote = FALSE)


table_clusb = cbind(eps_clusb,
                    round(pchisq(apply(clus_xi2b[, ],1,min),rar.num-1),2))
colnames(table_clusb) = c("eps_med","eps_low","eps_upp","p_value")
write.table(table_clusb, "outputs/Suppl_table_general_model_variant.txt", quote = FALSE)


#### 2-parameter general model for clusters ####
table_clusC = cbind(eps_clusC, round(gal_clusC,3),
                    round(pchisq(apply(clus_xi2C[, , ],1,min),rar.num-1),2))
colnames(table_clusC) = c("eps_med","eps_low","eps_upp","ggal_med","ggal_low","ggal_upp","p_value")
write.table(table_clusC, "outputs/Table1_cluster_model.txt", quote = FALSE)
par(mfrow=c(1,1),oma=c(3,3,3,3))
plot(clus_mass/clus_rad,1/eps_clusB[,1])
     

#### IMPORTANT: Relation between epsilon and mass density ####


cor(sqrt(smcGa07.mass/(4/3*pi*(smcGa07.R$x*kpc*4)^3)/(3/(8*pi*GN*T0^2))),galeps_0.48)
lm(galeps_0.48 ~ (sqrt(smcGa07.mass/(4/3*pi*(smcGa07.R$x*kpc*4)^3)/(3/(8*pi*GN*T0^2)))))
lm(c(galeps_0.48^2-1/6) ~ 0+c((smcGa07.mass/(4/3*pi*(smcGa07.R$x*kpc*4)^3)/(3/(8*pi*GN*T0^2)))))
lmlog = lm(c(log(galeps_0.48^2-1/6)) ~ c(log(smcGa07.mass/(4/3*pi*(smcGa07.R$x*kpc*4)^3)/(3/(8*pi*GN*T0^2)))))
summary(lmlog)

par(mfrow=c(1,1),oma=c(3,3,3,3))
plot(sqrt(smcGa07.mass/(4/3*pi*(smcGa07.R$x*kpc*3)^3)/(3/(8*pi*GN*T0^2))),galeps_0.5,
     xlab="Square root of mass density between 50 and 200 kpc",xlim=c(0,160),ylim=c(0,160))
points(sqrt(clus_massmin/(4/3*pi*((clus_radmin)^3))/(3/(8*pi*GN*T0^2))),eps_clusB[,1],pch=20) ##### ojuuuuu #####
abline(0,1)

par(mfrow=c(1,1),oma=c(3,3,3,3))
plot(sqrt(smcGa07.mass/(4/3*pi*(smcGa07.R$x*kpc*4)^3)/(3/(8*pi*GN*T0^2))),galeps_k,
     xlab="Square root of mass density between 50 and 200 kpc",log="xy",xlim=c(1,160),ylim=c(1,160))
points(sqrt(clus_massmin/(4/3*pi*((clus_radmin)^3))/(3/(8*pi*GN*T0^2))),eps_clusB[,1],pch=20) ##### ojuuuuu #####
abline(0,1)

#### RAR_FigureC3: Relation between epsilon and mass density ####
## summary(c(smcGa07.R$x*4,clus_radmin/kpc))

{
  
  pdf(paste0("outputs/Fig5.pdf"), width = 7, height = 4.5)
  {
par(fig=c(0,0.48,0.1,1),oma=c(5.9,3,0.5,3),mar=c(0,0,0,0))

    plot.new()
    rect(-10,-10,10,10,col=rgb(1,0.9,0.9,0.1))
    
par(fig=c(0,0.48,0.1,1),oma=c(5.9,3,0.5,3),mar=c(0,0,0,0),new=TRUE)
    
#### 1-par model for galaxies and clusters #####
smcGa07.VO = aggregate(smcGa07.V,by=list(namess),mean,na.rm=T)
galeps_k_err = apply(cbind(galeps_0.47,galeps_0.48,galeps_k),1,sd)/sqrt(2)
gammagc_k_err = apply(cbind(0.47*pi,0.48*pi,gammagc_k),1,sd)/2

xlab = expression(paste(sqrt(rho[typ]/rho[vac])))
ylab = expression(paste(epsilon[H]))

x0 = sqrt(smcGa07.mass/(4/3*pi*(smcGa07.R$x*kpc*4)^3)/(3/(8*pi*GN*T0^2)))
#x1 = sqrt(smcGa07.mass/(4/3*pi*(smcGa07.R$x*kpc*1.1)^3)/(3/(8*pi*GN*T0^2)))
#cor(x1^0.92, galeps_0.48); cor(x0, galeps_0.48); cor(x0-5/6, galeps_0.48);
plot(x0, galeps_0.48-5/6, lwd=1.5,
     xlab=xlab,ylab=ylab,xlim=c(3,180),ylim=c(3,180),log="xy",pch=0,axes=FALSE,
     col=rainbow(length(namesg))[order(smcGa07.VO$x,decreasing=TRUE)[-55]])

segments(x0,galeps_0.48-galeps_k_err,x0,galeps_0.48+galeps_k_err,
         col=rainbow(length(namesg))[order(smcGa07.VO$x,decreasing=TRUE)[-55]])

x0 = sqrt(clus_massmin/(4/3*pi*((clus_radmin)^3))/(3/(8*pi*GN*(1+clus_z)^(-1)*T0^2))) 
points(x0, eps_clusB[,1],
       pch=20,cex=1.5,col=rainbow(length(clusters))) ##### ojuuuuu #####
segments(x0,eps_clusB[,2],x0,eps_clusB[,3],
         col=rainbow(length(clusters)))

box()
abline(0,1)
axis(side=1,c(5,10,20,50,100),padj =-0.8)
axis(side=2,c(5,10,20,50,100),padj =0.8)
abline(v=c(5,10,20,50,100),h=c(5,10,20,50,100),lty=2,col="gray95")
mtext(side=1,xlab,line=1.99)
mtext(side=2,ylab,line=1.6)
#typ=4r \sim 50-200*kpc
#gamma[center]==0.48*pi

#### 2-par model for galaxies and clusters #####

par(fig=c(0.52,1,0.1,1),oma=c(5.9,3,0.5,3),mar=c(0,0,0,0),new=TRUE)

plot.new()
rect(-10,-10,10,10,col=rgb(1,0.9,0.9,0.1))

par(fig=c(0.52,1,0.1,1),oma=c(5.9,3,0.5,3),mar=c(0,0,0,0),new=TRUE)

xlab = expression(paste(log(1/epsilon[H])))
ylab = expression(paste(cos(gamma[cen])))

x0 = -log(galeps_k)[namesg>1]-abs(seq(0,0.05,length.out= sum(namesg>1)))
xer1 = -log(galeps_k*0.3)[namesg>1]
xer1[is.na(xer1)] = x0
plot(x0, cos(gammagc_k)[namesg>1],pch=0,xlab=xlab, ylab=ylab,axes=FALSE,ylim=c(0,0.16),
     col=rainbow(length(namesg))[order(smcGa07.VO$x,decreasing=TRUE)[-55]], lwd=1.5)
segments(xer1,cos(gammagc_k)[namesg>1],-log(galeps_k+galeps_k_err/2)[namesg>1],cos(gammagc_k)[namesg>1],
         col=rainbow(length(namesg))[order(smcGa07.VO$x,decreasing=TRUE)[-55]])
segments(x0,cos(gammagc_k-gammagc_k_err)[namesg>1],x0,cos(gammagc_k+gammagc_k_err)[namesg>1],
         col=rainbow(length(namesg))[order(smcGa07.VO$x,decreasing=TRUE)[-55]])
box()
abline(a=cos(0.46*pi),b=0.02)

points(-log(eps_clusC[-3,1]), cos(gal_clusC[-3,1]*pi),pch=20,cex=1.5,col=rainbow(length(clusters)))
segments(-log(eps_clusC[-3,1]), cos(gal_clusC[-3,2]*pi),-log(eps_clusC[-3,1]), cos(gal_clusC[-3,3]*pi),col=rainbow(length(clusters)))
segments(-log(eps_clusC[-3,2]), cos(gal_clusC[-3,1]*pi),-log(eps_clusC[-3,3]), cos(gal_clusC[-3,1]*pi),col=rainbow(length(clusters)))

axis(side=1,seq(-5,0,1),padj =-0.8)
axis(side=4,seq(0,0.14,0.02),padj =-0.8)
mtext(side=1,xlab,line=1.8)
mtext(side=4,ylab,line=1.6)
abline(v=seq(-5,0,1),h=seq(0,0.14,0.02),lty=2,col="gray95")

cor(-log(galeps_k)[namesg>1],cos(gammagc_k)[namesg>1])
lmg = lm(cos(gammagc_k[namesg>1]) ~ c(-log(galeps_k[namesg>1])))
#0.122(4)           0.0185(2)
summary(lmg )
acos(0.12219)/pi; acos(0.12219+0.004201)/pi; acos(0.12219-0.004201)/pi;
exp(0.12219/0.02)

par(fig=c(0,1,0,0.181),mar=c(0,0,0,0),oma=c(0.01,0.1,1,0.1),new=TRUE) 
plot.new()
legend("bottomleft",legend=clusters,pch=20,col=clus_col,ncol=2,pt.cex=cex0*0.61, cex=cex0*0.45, bty="n")
legend("bottomright",legend=names(namesG)[order(smcGa07.VO$x,decreasing=TRUE)[-55]],pch=0,col=rainbow(length(namesg))[-55],ncol=10,pt.cex=cex0*0.48, cex=cex0*0.365, bty="n")

par(fig=c(0,1,0,0.225),mar=c(0,0,0,0), oma=c(0,0.1,0,0.45), new=T)
plot.new()
legend("topleft",cex=0.801,legend=expression(bold("Clusters:")),pch=20,col="white", bty="n")
legend("topright",cex=0.801,legend=expression(bold("Galaxies:")),pch=20,col="white", bty="n", pt.cex=0.001)



##cos(gammagc_k) = cos(0.46*pi) - 0.02*log(galeps_k) ### epsilon entre 1 y 450
}
dev.off()
}

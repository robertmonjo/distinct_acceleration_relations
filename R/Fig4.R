# Monjo & Banik (2025), ApJ 992, 35 — reproducible pipeline
# Sourced in order by run_all.R (shared session state).
# Figure C2: fitting of the radial acceleration relation. -> outputs/Fig4.pdf




gal_xi2 = gal_xi2_0.48
##gal_xi2[gal_xi2>0.2 & !is.na(gal_xi2)] = 0.2
nonan = !is.na(gal_xi2[,1])
plot(seq_galeps, gal_xi2[nonan,][1,],ylim=c(0,0.2))
galeps_ok = rep(NA,length(namesg))
galeps_ko = rep(NA,length(namesg))

par(mfrow=c(10,7),mar=c(0,0,0,0),oma=c(3,3,3,3))
for(i in 1:sum(nonan))
{
  plot(seq_galeps, gal_xi2[nonan,][i,],ylim=c(0,1), xlim=c(1,120),axes=FALSE,typ="l", log="x")
  text(15, 0.45,names(namesg)[nonan][i],cex=0.8)
  box()
  
  galeps_ko[nonan][i] = seq_galeps[order(gal_xi2[nonan,][i,])[1]]
  galeps_ok[nonan][i] = galeps_ko[nonan][i] 
  
  abline(v=galeps_ok[nonan][i],col="red", lty=2)
  
  
  if((i-1)%%6==0)
    axis(side=2,seq(0.0,0.2,0.02), padj =0.8, las=2)
  if((i)%%6==0)
    axis(side=4,seq(0.0,0.2,0.02),padj =-0.8, las=2)
  if(i<=6)
    axis(side=3,c(1,2,10,20,50,100,200,500), padj =0.8 )
  if(i>=55)
    axis(side=1,c(1,2,10,20,50,100,200,500), padj =-0.8 )
  
}


#### 2. RAR ####

rar = read.table("data/clusterRAR.dat")
rar.names = c("Name", "z", "radius", "log(gbar)", "log(gtot)","err_log(gtot)_low","err_log(gtot)_up")
rar.units = c("","","(kpc)","[m/s^2]","[m/s^2]","[m/s^2]","[m/s^2]")
colnames(rar) = rar.names
rar = as.data.frame(rar)
clusters = as.character(unique(rar$Name))

rar.acc_total_up = 10^(rar$`log(gtot)`+ rar$`err_log(gtot)_up`)
rar.acc_total_do = 10^(rar$`log(gtot)`- rar$`err_log(gtot)_low`)
rar.acc_total_err = (rar.acc_total_up-rar.acc_total_do)/2
rar.acc_total = 10^rar$`log(gtot)`
rar.acc_newton = 10^rar$`log(gbar)`
rar.acc_newton_err = rar.acc_total_err/rar.acc_total*rar.acc_newton
rar.radius = rar$radius*kpc
rar.mass = rar.acc_newton*rar.radius^2/GN
rar.escape_newton = sqrt(2*rar.acc_newton*rar.radius)
rar.escape_hubble = rar.radius/T0*(1+rar$z)

diff_accel = T0/(c0*(1+rar$z))*(rar.acc_total-rar.acc_newton)
diff_accel_up = T0/(c0*(1+rar$z))*(rar.acc_total_up-rar.acc_newton)
diff_accel_do = T0/(c0*(1+rar$z))*(rar.acc_total_do-rar.acc_newton)
diff_accel_err = (diff_accel_up-diff_accel_do)/2

quot_accel = (rar.acc_total/rar.acc_newton)
quot_accel_up = (rar.acc_total_up/rar.acc_newton)
quot_accel_do = (rar.acc_total_do/rar.acc_newton)
quot_accel_err = (quot_accel_up-quot_accel_do)/2

gamma_U = pi/3
gamma_gc = pi/2
vE2vH2 = (rar.escape_newton/rar.escape_hubble)^2
sin2gamma_Mgal = sin(gamma_U)^2 + (sin(gamma_gc)^2 - sin(gamma_U)^2)*1/(1+rar.escape_newton^2/rar.escape_hubble^2)
gamma_Mgal = asin(sqrt(sin2gamma_Mgal))
gamma_0gal = gamma_Mgal/cos(gamma_Mgal)
gamma_MM = acos.(T0/((1+rar$z)*c0)*(rar.acc_total-rar.acc_newton))
plot(vE2vH2, 1/gamma_MM,log="x")

gamma_M2 = acos.(T0/((1+rar$z)*c0)*(rar.acc_newton*((rar.acc_total/rar.acc_newton)^2-1)))
plot(vE2vH2, gamma_M2,log="x")

#save.image("outputs/figure_Milgrom5.RData")
#save.image("outputs/figure_Milgrom6.RData")
#load("data/figure_Milgrom4.RData")
#load("data/figure_Milgrom5.RData")


plot(((rar.escape_newton-rar.escape_hubble)/rar.escape_hubble), (pi/2-gamma_MM)/(pi/2-pi/3),log="xy")
plot((rar.escape_newton/rar.escape_hubble), (pi/2-gamma_MM)/(pi/2-pi/3),log="xy")


#####  Entonces ######
par(mfrow=c(1,1))
plot(((rar.escape_newton/rar.escape_hubble)^2)/(55^2+(rar.escape_newton/rar.escape_hubble)^2), (sin(pi/2)^2-sin(gamma_MM)^2)/(sin(pi/2)^2-sin(pi/3)^2),log="xy")
abline(0,1)

eps_H = sqrt(((rar.escape_newton/rar.escape_hubble)^2)/((sin(pi/2)^2-sin(gamma_MM)^2)/(sin(pi/2)^2-sin(pi/3)^2))-(rar.escape_newton/rar.escape_hubble)^2)
quantile(eps_H,na.rm=T,c(0.10,0.5,0.9)) #45 \pm 14

seq_gal = c(seq(0.4,0.498,length.out=98),0.499, 0.5)
seq_eps = unique(c(seq(0.01,1,0.01),seq(1,10,0.1),seq(10,50,0.5),seq(50,100,1),seq(100,200,2),seq(200,1000,20),seq(1000,10000,200)))

clus_xi2g0 = array(NA, dim=c(length(clusters), length(seq_eps)),
                 dimnames = list(clusters, seq_eps))

clus_xi2A = array(NA, dim=c(length(clusters), length(seq_eps)),
                         dimnames = list(clusters, seq_eps))
clus_xi2B = array(NA, dim=c(length(clusters), length(seq_eps)),
                  dimnames = list(clusters, seq_eps))
clus_xi2b = array(NA, dim=c(length(clusters), length(seq_eps)),
                  dimnames = list(clusters, seq_eps))
clus_xi2C = array(NA, dim=c(length(clusters), length(seq_eps),length(seq_gal)),
                  dimnames = list(clusters, seq_eps, seq_gal))


gempty_space = pi/3 ## seq_emp[iemp]*pi   #0.739 #pi/5
g_galaxy =  0.49*pi  #### seq_gal[igal]*pi
gblack_hole = 0.5*pi

for(ieps in 1:length(seq_eps))
{

  eps_H0 = seq_eps[ieps]
  eps_H0B = seq_eps[ieps]
   
  
  vEvH = (rar.escape_newton/rar.escape_hubble)
  vE2vH2 = vEvH^2
  
  dens_div_dens = abs(6/5*vEvH^2-eps_H0B^2)/(eps_H0B^2+vEvH^2)
  gamma_M = asin(sqrt(6/5*sin(gempty_space)^2-dens_div_dens*(6/5*sin(gempty_space)^2-sin(gblack_hole)^2)))
  gamma_0A= gamma_M/cos(gamma_M)
  pred_accelA = ((1+rar$z)*c0)/(T0*gamma_0A*rar.acc_newton)

  dens_div_dens = abs(vEvH^2-eps_H0B^2)/(eps_H0B^2+vEvH^2)
  gamma_M = asin(sqrt(sin(gempty_space)^2-dens_div_dens*(sin(gempty_space)^2-sin(gblack_hole)^2)))
  gamma_0b= gamma_M/cos(gamma_M)
  pred_accelb = ((1+rar$z)*c0)/(T0*gamma_0b*rar.acc_newton)
  
  dens_div_dens = abs(vEvH^2-eps_H0B^2)/(eps_H0B^2+vEvH^2)
  gamma_M = asin(sqrt(sin(gempty_space)^2-dens_div_dens*(sin(gempty_space)^2-sin(g_galaxy)^2)))
  gamma_0B= gamma_M/cos(gamma_M)
  pred_accelB = ((1+rar$z)*c0)/(T0*gamma_0B*rar.acc_newton)
  
  for(iclu in 1:length(clusters))
  {
    lgclus = rar$Name == clusters[iclu]
    clus_xi2g0[iclu, ieps] = sum((1/gamma_0B[lgclus]-diff_accel[lgclus])^2/diff_accel_err[lgclus]^2)
    clus_xi2A[iclu, ieps] = sum((pred_accelA[lgclus]-quot_accel[lgclus])^2/quot_accel_err[lgclus]^2)
    clus_xi2B[iclu, ieps] = sum((pred_accelB[lgclus]-quot_accel[lgclus])^2/quot_accel_err[lgclus]^2)
    clus_xi2b[iclu, ieps] = sum((pred_accelb[lgclus]-quot_accel[lgclus])^2/quot_accel_err[lgclus]^2)
  }
  
  for(igal in 1:length(seq_gal))
  {
    g_galaxyC = seq_gal[igal]*pi ####pi/2
    
    dens_div_dens = vE2vH2/(eps_H0B^2+vE2vH2)
    gamma_m2 = asin(sqrt(sin(g_galaxyC)^2-dens_div_dens*(sin(g_galaxyC)^2-sin(gempty_space)^2)))
    gamma_0C = gamma_m2/cos(gamma_m2)
    pred_accelC = ((1+rar$z)*c0)/(T0*gamma_0C*rar.acc_newton)
    
    for(iclu in 1:length(clusters))
    {
      lgclus = rar$Name == clusters[iclu]
      clus_xi2C[iclu, ieps, igal] = sum((pred_accelC[lgclus]-quot_accel[lgclus])^2/quot_accel_err[lgclus]^2)
    }
  }
  
  cat(ieps,"/", length(seq_eps)," ",igal,"/",length(seq_gal),"\n")
   cat(range(clus_xi2B[, ieps]),"\n")
}
apply(clus_xi2g0[,  ],1,min,na.rm=T)
apply(clus_xi2C[, , ],1,min,na.rm=T)
apply(clus_xi2B[,  ],1,min,na.rm=T)
apply(clus_xi2b[,  ],1,min,na.rm=T)
apply(clus_xi2A[,  ],1,min,na.rm=T)

clus_xi2C[clus_xi2C>500] = 500
rar.num = rep(NA,length(clusters))
names(rar.num) = clusters
for(iclu in 1:length(clusters))
  rar.num[iclu] = sum(rar$Name == clusters[iclu])
  
p67 = qchisq(0.67,rar.num)
p95 = qchisq(0.95,rar.num)
p99 = qchisq(0.99,rar.num)
p999 = qchisq(0.999,rar.num)
p67S = qchisq(0.67,sum(rar.num))
p95S = qchisq(0.95,sum(rar.num))

pchisq(apply(clus_xi2B,1,min, na.rm=T),rar.num-1)
pchisq(apply(clus_xi2b,1,min, na.rm=T),rar.num-1)
pchisq(apply(clus_xi2C,1,min, na.rm=T),rar.num-1)
pchisq(apply(clus_xi2A,1,min, na.rm=T),rar.num-1)



### RAR_FigureC2 ##### 
######  Fitting epsilon and gamma ####
###### Ojuuu: pasar a Xi2 frente a epsilon ####

eps_clusA = array(NA,dim=c(length(clusters),3), dimnames=list(clusters,c("med","low","upp")))
eps_clusB = array(NA,dim=c(length(clusters),3), dimnames=list(clusters,c("med","low","upp")))
eps_clusb = array(NA,dim=c(length(clusters),3), dimnames=list(clusters,c("med","low","upp")))
eps_clusC = array(NA,dim=c(length(clusters),3), dimnames=list(clusters,c("med","low","upp")))
gal_clusC = array(NA,dim=c(length(clusters),3), dimnames=list(clusters,c("med","low","upp")))
clus_mass  = array(NA,dim=length(clusters),dimnames = list(clusters))
clus_z  = array(NA,dim=length(clusters),dimnames = list(clusters))
clus_massmin  = array(NA,dim=length(clusters),dimnames = list(clusters))
clus_rad = array(NA,dim=length(clusters),dimnames = list(clusters))
clus_radmin = array(NA,dim=length(clusters),dimnames = list(clusters))
clus_ve = array(NA,dim=c(length(clusters),3), dimnames=list(clusters,c("med","low","upp")))
clus_vH = array(NA,dim=c(length(clusters),3), dimnames=list(clusters,c("med","low","upp")))
clus_vevH2 = array(NA,dim=c(length(clusters),3), dimnames=list(clusters,c("med","low","upp")))
{
  pdf("outputs/Fig4.pdf", width = 9, height = 5)
  {
    {
      par(mfrow=c(2,5), mar=c(0,0,0,0), oma=c(3.5,3.5,2,30.3))
      for(iclu in 1:10)
      {
        rv = order(1/seq_eps)
        u_seq_eps = rev(1/seq_eps)
        image(u_seq_eps,seq_gal,clus_xi2C[iclu,rv,]/p95[iclu],col=hcl.colors(100, "terrain"),ylim=c(0.407,0.501),xlim=c(0.0,0.062),zlim=c(0,1),axes=FALSE)
        contour(u_seq_eps,seq_gal,clus_xi2C[iclu,rv,],levels = c(p67[iclu],p95[iclu]),labels = c("67%","95%"), add=T)
        contour(u_seq_eps,seq_gal,clus_xi2C[iclu,rv,],levels = c(p99[iclu]),labels = c(""),lty=2, add=T)
        if(iclu == 3)
          contour(u_seq_eps,seq_gal,clus_xi2C[iclu,rv,],levels = c(200),labels = c(""),lty=2, add=T)
        points(u_seq_eps[order(apply(clus_xi2C[iclu,rv,],1,min))[1]],
               seq_gal[order(apply(clus_xi2C[iclu,rv,],2,min))[1]], pch=20,cex=2,col="red")
        
        text(1/63,0.412,clusters[iclu])
        ## text(1/63,0.412,paste0(clusters[iclu], " (",rar.num[iclu],")"))
        abline(h=seq(0.4,0.5,0.02),lty=2,col="gray92")
        abline(v=seq(0,0.1,0.02),lty=2,col="gray92")
        
        box()
        if(iclu %in% c(1,6))
          axis(2, seq(0.4,0.5,0.02),padj =0.8 )
        
        #if(iclu %in% c(5,10))
        #  axis(4, seq(0.4,0.5,0.02),padj =-0.8 )
        
        if(iclu %in% c(1:5))
        {
          axis(3,c(0, 0.01,0.02,0.03,0.04,0.05), paste(c(0, 0.01,0.02,0.03,0.04,0.05)),padj =0.8 ) # axis(3, seq(0,200,50),padj =0.8 )
          axis(3,0.1,padj =0.8 )
        }

        if(iclu %in% c(6:10))
        {
          axis(1,c(0, 0.01,0.02,0.03,0.04,0.05), paste(c(0, 0.01,0.02,0.03,0.04,0.05)),padj =-0.8 ) # axis(1, seq(0,200,50),padj =-0.8 )
          axis(1,0.1,padj =-0.8 )
        }
        
        
      }
      par(fig=c(0,0.65,0,1),oma=c(3.5,3.5,2,2),mar=c(0,0,0,0),new=TRUE)  
      plot.new()
      mtext(side=1,expression(paste("Square root of the relative density (", 1/epsilon[H],")")),line=2,cex=0.95)
      mtext(side=2,expression(paste("Dominant galaxy projective angle ( ",gamma[cen]/pi," )")),line=1.6,cex=0.95)
      
    }
    clus_col = rainbow(length(clusters))
    par(fig=c(0.61,0.95,0.45,1), mar=c(0,0,0,0), oma=c(3,1,2,3.4),new=TRUE)
    {
      
      clus_xi2 = clus_xi2b #g_center = pi/2
      
      plot(u_seq_eps,pchisq(clus_xi2[1,rv],rar.num[1]-1),ylim=c(0,1),xlim=c(0.005,0.2),log="x",axes=FALSE,typ="l",lwd=2,col="white")  #col=hcl.colors(100, "terrain"),
      for(iclu in 1:10)
      {
        xi2_norm = pchisq(clus_xi2[iclu,rv],rar.num[iclu]-1)
        #xi2_norm = clus_xi2B[iclu,]/max(clus_xi2B[iclu,])
        lines(u_seq_eps,xi2_norm,lwd=0.8,lty=2,col=clus_col[iclu])  #col=hcl.colors(100, "terrain"),
        points(u_seq_eps[order(xi2_norm)[1]], min(xi2_norm,na.rm=T), pch=20,cex=3,col=clus_col[iclu])
        if(iclu==3)
          points(u_seq_eps[order(clus_xi2[iclu,])[1]], 1, pch=20,cex=3,col=clus_col[iclu])
      }
      abline(v=1/c(1,2,5,10,20,50,100,200),lty=2,col="gray92")
      abline(h=seq(0,1,0.1),lty=2,col="gray92")
      box()
      #axis(2, seq(0,1,0.1),padj =0.8 )
      axis(4, seq(0,1,0.1),padj =-0.8 )
      axis(3, 1/c(1,2,5,10,20,50,100,200,500,1000),padj =0.8 )
      axis(1, 1/c(1,2,5,10,20,50,100,200,500,1000),padj =-0.8 )
      mtext(side=4,expression(paste(chi^2," p-value")),line=2.0,cex=0.85)
      mtext(side=1,expression(paste("Square root of the relative density (", 1/epsilon[H],")")),line=2,cex=0.85)
      
    }
    par(fig=c(0.1,0.615,0,0.35), mar=c(0,0,0,0), oma=c(3.5,3.5,2,2),new=TRUE)
    image.plot(seq_eps,seq_gal,clus_xi2C[iclu,,]/p95[iclu],col=hcl.colors(100, "terrain"),legend.only=TRUE, ylim=c(0.407,0.501),xlim=c(5,200),zlim=c(0,1),axes=FALSE)
    par(fig=c(0.42,0.575,0,0.30), mar=c(0,0,0,0), oma=c(3.5,3.5,2,2),new=TRUE)
    plot.new()
    mtext(side=4,expression(paste(chi^2," p-value")),line=1.6)
    
    par(fig=c(0.71,0.9,0,0.44),mar=c(0,0,0,0), new=T,oma=c(0,1,0,1), new=T)
    plot.new()
    #legend("topright",legend=c(expression(paste(gamma[galaxy]==pi/2,", ",gamma[neig]==pi/3,", ",epsilon[H]==57))),col=c("gray80"),lty=c(1),lwd=c(3), bty="n")
    legend("left",legend=clusters,pch=20,col=clus_col,ncol=2,pt.cex=1.5, cex=0.87, bty="n")
    
    par(fig=c(0.71,0.9,0,0.38),mar=c(0,0,0,0), new=T,oma=c(0,1,0,1), new=T)
    plot.new()
    legend("topleft",legend=expression(bold("Clusters:")),pch=20,col="white", bty="n")
    
    
    # par(fig=c(0,1,0.93,1),mar=c(0,0,0,0),oma=c(0,0,0,0), new=T)
    # plot.new()
    # text(-0.02,0.55,substitute(paste(bold("a) "))), cex=1)
    # text(0.62,0.55,substitute(paste(bold("b) "))),cex=1)
    
  }
  
  dev.off()
 }

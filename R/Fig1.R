# Monjo & Banik (2025), ApJ 992, 35 — reproducible pipeline
# Sourced in order by run_all.R (shared session state).
# Figure 1: cluster radial acceleration relation. -> outputs/Fig01_cluster_RAR.pdf




    if(FALSE)
    {
    plot(seq_eps,pchisq(clus_xi2[1,],rar.num[1]-1),col="white",ylim=c(0,1),xlim=c(5,1000),log="x",axes=FALSE,typ="l",lwd=2)  #col=hcl.colors(100, "terrain"),
    #plot(seq_galeps, gal_xi2[nonan,][1,],col="white",ylim=c(0,0.2), xlim=c(1,380),axes=FALSE,typ="l", log="x")
    box()
    for(i in 1:sum(nonan))
    {
      lines(seq_galeps, gal_xi2[nonan,][i,])
      text(15, 0.45,names(namesg)[nonan][i],cex=0.8)
      
      
      galeps_ko[nonan][i] = seq_galeps[order(gal_xi2[nonan,][i,])[1]]
      galeps_ok[nonan][i] = galeps_ko[nonan][i] 
      
      points(galeps_ok[nonan][i],min(gal_xi2[nonan,][i,]),col="red", lty=2)
    }
    axis(2, seq(0,1,0.1),padj =0.8 )
    axis(4, seq(0,1,0.1),padj =-0.8 )
    axis(3, c(1,2,5,10,20,50,100,200,500,1000),padj =0.8 )
    axis(1, c(1,2,5,10,20,50,100,200,500,1000),padj =-0.8 )
    
    par(fig=c(0,1,0,1),oma=c(3.5,3.5,2,2),mar=c(0,0,0,0),new=TRUE)  
    plot.new()
    mtext(side=1,expression(paste("Square root relative neighborhood density (", epsilon[H],")")),line=2)
    mtext(side=2,expression(paste(chi^2," p-value")),line=1.6)

    }
    


g_galaxy = 0.5*pi
gblack_hole = 0.5*pi
gempty_space = pi/3
gempty_space = pi/3

##### 0.333*pi model ####
eps_H0 = median(eps_H,na.rm=T)/2
vE2vH2 = (rar.escape_newton/rar.escape_hubble)^2
dens_div_dens = (rar.escape_newton/rar.escape_hubble)^2/(eps_H0^2+(rar.escape_newton/rar.escape_hubble)^2)
gamma_m1 = asin(sqrt(sin(g_galaxy)^2-dens_div_dens*(sin(gblack_hole)^2-sin(gempty_space)^2)))
gamma_01 = gamma_m1/cos(gamma_m1)
##  plot(c0/(T0*gamma_01*rar.acc_newton),quot_accel); abline(0,1)

vE2vH2_sort = seq(min(vE2vH2 ), max(vE2vH2), length.out = 300)
vEvH_sort = sqrt(vE2vH2_sort)

vE2vH2_sort2 = 10^seq(0, 8, length.out = 500)
vEvH_sort2 = sqrt(vE2vH2_sort2)


dens_div_dens = vE2vH2_sort/(mean(eps_clusC[-3,1])^2+vE2vH2_sort)
gamma_Mmed = asin(sqrt(sin(g_galaxy)^2-dens_div_dens*(sin(gblack_hole)^2-sin(gempty_space)^2)))
gamma_0med = gamma_Mmed/cos(gamma_Mmed)
dens_div_dens = vE2vH2_sort/(mean(eps_clusC[-3,3])^2+vE2vH2_sort)
gamma_Mup = asin(sqrt(sin(g_galaxy)^2-dens_div_dens*(sin(gblack_hole)^2-sin(gempty_space)^2)))
gamma_0up = gamma_Mup/cos(gamma_Mup)
dens_div_dens = vE2vH2_sort/(mean(eps_clusC[-3,2])^2+vE2vH2_sort)
gamma_Mdow  = asin(sqrt(sin(g_galaxy)^2-dens_div_dens*(sin(gblack_hole)^2-sin(gempty_space)^2)))
gamma_0dow = gamma_Mdow/cos(gamma_Mdow)


#dens_div_dens = (sqrt(vE2vH2_sort)/25-1)^2/(1+vE2vH2_sort/25^2)
#dens_div_dens = abs(vE2vH2_sort/mean(eps_clusB[,1])^2-1)/(1+vE2vH2_sort/mean(eps_clusB[,1])^2)

eps_min = mean(eps_clusB[-3,2])
eps_max = mean(eps_clusB[-3,3])
eps_seq = seq(eps_min,eps_max,length.out=10)
gamma_0mdd_seq = array(NA,dim=c(length(vEvH_sort),length(eps_seq)))
gamma_0mdd2_seq = array(NA,dim=c(length(vEvH_sort2),length(eps_seq)))
for(iseq in 1:length(eps_seq))
{
  dens_div_dens = abs(vEvH_sort^2-eps_seq[iseq]^2)/(eps_seq[iseq]^2+vEvH_sort^2)
  gamma_M = asin(sqrt(sin(gempty_space)^2-dens_div_dens*(sin(gempty_space)^2-sin(g_galaxy)^2)))
  gamma_0mdd_seq[,iseq] = gamma_M/cos(gamma_M)
  
  dens_div_dens2 = abs(vEvH_sort2^2-eps_seq[iseq]^2)/(eps_seq[iseq]^2+vEvH_sort2^2)
  gamma_M2 = asin(sqrt(sin(gempty_space)^2-dens_div_dens2*(sin(gempty_space)^2-sin(g_galaxy)^2)))
  gamma_0mdd2_seq[,iseq] = gamma_M2/cos(gamma_M2)
}

#lines(sqrt((vE2vH2_sort/4)^2)/eps0,1/gamma_0mdd, lwd=2,col="blue", lty=1)
gamma_0mdd = apply(gamma_0mdd_seq,1,quantile,0.5)
gamma_0mup = apply(gamma_0mdd_seq,1,quantile,0.9)
gamma_0mdo = apply(gamma_0mdd_seq,1,quantile,0.1)

gamma_0mdd2 = apply(gamma_0mdd2_seq,1,quantile,0.5)
gamma_0mup2 = apply(gamma_0mdd2_seq,1,quantile,0.9)
gamma_0mdo2 = apply(gamma_0mdd2_seq,1,quantile,0.1)




#https://milde.users.sourceforge.net/LUCR/Math/unimathsymbols.pdf
#### RAR_Figure01 (clusters) ####
#save.image("outputs/figure_Milgrom5.RData")
#load("data/figure_Milgrom4.RData")
#load("data/figure_Milgrom5.RData")
{

  eps0 = mean(eps_clusC[,1])
  eps0 = quantile(c(eps_clusB[,1]), 0.5, na.rm=TRUE)
  
  pdf(paste0("outputs/Fig01_cluster_RAR.pdf"), width = 10, height = 5.05)
  {
    cex0 = 1.5
    gamma_00 = gamma_01
    
      gempty_space = pi/3
      xmin=2
      xtext = "0.333"
      lty0 = 2
    
      
    
      par(fig=c(0,0.5,0.0941,1), mar=c(3,3,0.5,0.5))#TimesNewRoman"))
      
      plot(((1+rar$z)*c0)/(T0*gamma_01*rar.acc_newton),quot_accel,col="white",xlab="",ylab="",axes=FALSE, xlim=c(xmin,30), ylim=c(xmin,30),log="xy")
      clusters = as.character(unique(rar$Name))
      clus_col = rainbow(length(clusters))
      ##polygon_(seq(0,40,5), seq(0,40,5), seq(0,40,5)/1.2, seq(0,40,5)*1.2,add=TRUE,col=rgb(0.8,0.7,0.9,0.3))
      mtext(side=1,expression(paste("",italic(a)[Tot]/italic(a)[N],"  modelled")),line=1.8) ##, gamma[neigh]==0.333*pi)),line=1.8)
      
      mtext(side=2,expression(paste("",italic(a)[Tot]/italic(a)[N],"  observed")),line=1.6)
      axis(side=1,c(2,3,5,7,10,15,20,30),padj =-0.8 )
      axis(side=2,c(2,3,5,7,10,15,20,30),padj =0.8 )
      abline(v=c(2,3,5,7,10,15,20,30),h=c(2,3,5,7,10,15,20,30),lty=2,col="gray95")
      #abline(0,1, lwd=3/lty0,col=paste0("gray",80/lty0), lty=lty0)
      abline(0,1, lwd=2,col=rgb(0.8,0.4,0.4,1),lty=2)
      box()
      
      for(iclu in 1:length(clusters))
      {
        eps_H0 = eps_clusC[iclu,1]
        eps_H0_do = eps_clusC[iclu,2]
        eps_H0_up = eps_clusC[iclu,3]
        g_galaxy = gal_clusC[iclu,1]*pi
        
        vE2vH2 = (rar.escape_newton/rar.escape_hubble)^2
        dens_div_dens = (rar.escape_newton/rar.escape_hubble)^2/(eps_H0^2+(rar.escape_newton/rar.escape_hubble)^2)
        gamma_m0 = asin(sqrt(sin(g_galaxy)^2-dens_div_dens*(sin(g_galaxy)^2-sin(gempty_space)^2)))
        gamma_00 = gamma_m0/cos(gamma_m0)
        
        eps_H0 = eps_H0_do
        vE2vH2 = (rar.escape_newton/rar.escape_hubble)^2
        dens_div_dens = (rar.escape_newton/rar.escape_hubble)^2/(eps_H0^2+(rar.escape_newton/rar.escape_hubble)^2)
        gamma_m2 = asin(sqrt(sin(g_galaxy)^2-dens_div_dens*(sin(g_galaxy)^2-sin(gempty_space)^2)))
        gamma_do = gamma_m2/cos(gamma_m2)
        
        eps_H0 = eps_H0_up
        vE2vH2 = (rar.escape_newton/rar.escape_hubble)^2
        dens_div_dens = (rar.escape_newton/rar.escape_hubble)^2/(eps_H0^2+(rar.escape_newton/rar.escape_hubble)^2)
        gamma_m2 = asin(sqrt(sin(g_galaxy)^2-dens_div_dens*(sin(g_galaxy)^2-sin(gempty_space)^2)))
        gamma_up = gamma_m2/cos(gamma_m2)
        
        lgclus = rar$Name == clusters[iclu]
        points((((1+rar$z)*c0)/(T0*gamma_00*rar.acc_newton))[lgclus], quot_accel[lgclus], pch=20, cex=cex0, col= clus_col[iclu])
        segments((((1+rar$z)*c0)/(T0*gamma_00*rar.acc_newton))[lgclus], quot_accel_do[lgclus], (((1+rar$z)*c0)/(T0*gamma_00*rar.acc_newton))[lgclus], quot_accel_up[lgclus], col= clus_col[iclu])
        segments((((1+rar$z)*c0)/(T0*gamma_do*rar.acc_newton))[lgclus], quot_accel[lgclus], (((1+rar$z)*c0)/(T0*gamma_up*rar.acc_newton))[lgclus], quot_accel[lgclus], col= clus_col[iclu])
        
      }
    
      
    eps00 = eps0/sqrt(2)
 
    par(fig=c(0.5,1,0.0941,1), new=T, mar=c(3,3,0.5,0.5))
    plot(vEvH/eps0*eps00, diff_accel,col="white",xlab="",ylab="",axes=FALSE, log="xy",ylim=c(0.004,1.5))
    box()
    clusters = as.character(unique(rar$Name))
    clus_col = rainbow(length(clusters))

    uno = rep(1,length(vE2vH2_sort))
    #polygon_(sqrt(vE2vH2_sort)/eps0, 1/gamma_0med,  1/gamma_0dow,1/gamma_0up,add=TRUE,col=rgb(0.8,0.7,0.9,0.3))
    
    #polygon_(c(40-4*2,40+8*2), c(1,1),   c(0.001,0.001),  c(2,2), add=TRUE,col=rgb(0.65,0.99,0.12,0.18))
    #polygon_(c(56-12*2,56+22*2), c(1,1),   c(0.001,0.001),  c(2,2), add=TRUE,col=rgb(0.65,0.99,0.12,0.18))
    
    polygon_(sqrt(vE2vH2_sort)/eps0*eps00, uno*cos(0.466*pi)/(0.466*pi),  uno*cos(0.477*pi)/(0.477*pi), uno*cos(0.456*pi)/(0.456*pi), add=TRUE,col=rgb(0.12,0.65,0.99,0.08))
    polygon_(sqrt(vE2vH2_sort)/eps0*eps00, 1/gamma_0mdd,  1/gamma_0mdo,1/gamma_0mup,add=TRUE,col=rgb(0.75,0.75,0.75,0.3))
    polygon_(sqrt(vE2vH2_sort)/eps0*eps00, 1/gamma_0med,  1/gamma_0dow,1/gamma_0up,add=TRUE,col=rgb(0.8,0.4,0.4,0.3))
    #mtext(side=1,expression(paste("Newton-Hubble speed rate [", sqrt(2)*v[N]/v[H],"]")),line=1.8)
    mtext(side=1,expression(v[N]/v[H]),line=1.8)
    mtext(side=2,expression(paste("[ ",italic(a)[Tot]-italic(a)[N]," ]/[ c/t ]")),line=1.6)
    axis(side=1,c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200),
         paste(c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200)),padj =-0.8 )
    abline(v=c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100),h=c(1,2,5,10,20,50,100)/100,lty=2,col="gray95")
    axis(side=2,c(1,5)/100,padj =0.8 )
    axis(side=2,2/100,padj =0.8 )
    axis(side=2,c(0.1,0.2,0.5),padj=0.8)
    axis(side=2,1,padj=0.8)
    
    for(iclu in 1:length(clusters))
    {
      lgclus = rar$Name == clusters[iclu]
      eps_H0 = eps_clusB[iclu,1]
      eps_H0_do = eps_clusB[iclu,2]
      eps_H0_up = eps_clusB[iclu,3]
      eps00 = eps_H0/sqrt(2)
      points(vEvH[lgclus]/eps_H0*eps00 , diff_accel[lgclus], pch=20, cex=cex0, col= clus_col[iclu])
      segments(vEvH[lgclus]/eps_H0*eps00 , diff_accel_do[lgclus], vEvH[lgclus]/eps_H0*eps00 , diff_accel_up[lgclus], col = clus_col[iclu])
    }
    #legend("topright",legend=clusters,pch=20,col=clus_col, bty="n")
    #points(rar.escape_newton/rar.escape_hubble, diff_accel,pch=21)
    
    legend("bottomright",ncol=2,legend=c("MOND-like approach","General model", "Cluster approach"),cex=0.90, col=c(rgb(0.12,0.65,0.99,0.08),rgb(0.75,0.75,0.75,0.3),rgb(0.8,0.4,0.4,0.3)),lty=c(1,1,1),lwd=c(6, 6.5,6), bty="n")
    legend("bottomright",ncol=2,legend=c("MOND-like approach","General model", "Cluster approach"),cex=0.90,col=c(rgb(0.12,0.65,0.99,0.2),"gray70",rgb(0.8,0.4,0.4,1)),lty=c(1,1,2),lwd=c(1.8,1.8,1.8), bty="n")
    
  
    abline(h=c(0.5,1),col="gray65",lwd=1,lty=2)
    text(2.15,0.55,"Empty-space limit",cex=0.6,col="gray75")
    text(2.05,1.10,"Causality limit",cex=0.6,col="gray75")
    ## text(0.052,0.09,expression(paste("Milgrom's value for ",gamma[0]^{-1})),cex=0.6,col="gray75")
    ##text(0.048,0.092,"Milgrom's",cex=0.55,col="gray75")
    ##text(0.050,0.072,expression(paste("value for ",gamma[0]^{-1})),cex=0.55,col="gray75")
    
    eps00 = eps0/sqrt(2)
    ##lines(sqrt(vE2vH2_sort)/80,1/gamma_0med, lwd=3,col="gray70")
    lines(sqrt(vE2vH2_sort)/eps0*eps00,  uno*cos(0.466*pi)/(0.466*pi), lwd=2, col=rgb(0.12,0.65,0.99,0.2),lty=1)
    lines(sqrt(vE2vH2_sort)/eps0*eps00,1/gamma_0mdd, lwd=2,col="gray70", lty=1)
    lines(sqrt(vE2vH2_sort)/eps0*eps00,1/gamma_0med, lwd=2,col=rgb(0.8,0.4,0.4,1), lty=2)
    #lines(sqrt(vE2vH2_sort)/eps0,1/gamma_0mup, lwd=1,col="gray40", lty=2)
    #lines(sqrt(vE2vH2_sort)/eps0,1/gamma_0mdo, lwd=1,col="gray40", lty=2)
    ###lines(sqrt(vE2vH2_sort),1/gamma_0up)
    ###lines(sqrt(vE2vH2_sort),1/gamma_0dow)
    
    par(fig=c(0,1,0,0.0991),mar=c(0,0,0,0), new=T,oma=c(0,1,0,1), new=T)
    plot.new()
    legend("right",legend=clusters,pch=20,col=clus_col,ncol=10,pt.cex=cex0, bty="n")
    #legend("right",ncol=2,legend=c("MOND-like approach","General model", "Cluster approach"),cex=0.90, col=c(rgb(0.12,0.65,0.99,0.08),rgb(0.75,0.75,0.75,0.3),rgb(0.8,0.4,0.4,0.3)),lty=c(1,1,1),lwd=c(6, 6.5,6), bty="n")
    #legend("right",ncol=2,legend=c("MOND-like approach","General model", "Cluster approach"),cex=0.90,col=c(rgb(0.12,0.65,0.99,0.2),"gray70",rgb(0.8,0.4,0.4,1)),lty=c(1,1,2),lwd=c(1.8,1.8,1.8), bty="n")
    
    #par(fig=c(0,1,0,0.1641),mar=c(0,0,0,0), new=T,oma=c(0,1,0,1), new=T)
    #plot.new()
    legend("left",legend=expression(bold("Clusters:")),pch=20,col="white", bty="n")
    #legend("topright",legend=expression(bold("Hyperconical modified gravity (HMG):")),pch=20,col="white", bty="n")
    
    # par(fig=c(0,1,0.93,1),mar=c(0,0,0,0),oma=c(0,0,0,0), new=T)
    # plot.new()
    # text(-0.02,0.55,"a) ",cex=1.5)
    # text(0.51,0.55,"b)",cex=1.5)
    
  }
  dev.off()
}

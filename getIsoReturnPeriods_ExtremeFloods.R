###=====================================================================================================
### Guillaume Evin
### 16/03/2020, Grenoble
###  IRSTEA
### guillaume.evin@inrae.fr
### 
### Evin & Piton (Mars 2020): Codes d?velopp?s dans le cadre de l'action SNRH CrueBiv - INRAE, UR ETNA
### "Analyse bivari?e des liens entre magnitude et dur?e des crues en zones Alpine et Pyr?n?enne"
###=====================================================================================================

###______________________________________
### getIsoReturnPeriods_ExtremeFloods.R
###______________________________________
### Code extraction de contours equi-probables pour diff?rentes p?riodes de retour

# REFERENCES

# Serinaldi, Francesco. 2015. "Dismissing Return Periods!" Stochastic Environmental Research and Risk Assessment
# 29 (4): 1179-89. https://doi.org/10.1007/s00477-014-0916-1.

#Salvadori, G., G. R. Tomasicchio, and F. D'Alessandro. 2014. "Practical Guidelines for Multivariate Analysis and 
#Design in Coastal and Off-Shore Engineering." Coastal Engineering 88 (June): 1-14. 
#https://doi.org/10.1016/j.coastaleng.2014.01.011.


#======================================================================
# load outputs of previous treatments
#======================================================================

# load outputs from the bivariate fit:
# - listCop,listQcdf,listDcdf,listFitQ,listFitD,PkBar,PkBar,seqt
load(paste0(path.datacrues,"fitBivDist.RData"))

# load stations Metadata
load(paste0(path.datacrues,"S1/Metadata.stations.S1.RData"))
nS = nrow(Metadata.stations.S1)

# contour in [0,1]^2 space
nUV = 1000
sequ = seqv = seq(from=0,to=1,length.out=nUV)

# return periods
vecRP = c(2,5,10,20)
nRP = length(vecRP)

# Eq 11 in Salvadori, G., C. De Michele, and F. Durante. 2011. "On the Return Period and Design 
#in a Multivariate Framework." Hydrol. Earth Syst. Sci. 15 (11): 3293-3305. https://doi.org/10.5194/hess-15-3293-2011.
#
# pk = 1-Kc(t) = mu (mean elapsed time between two events) divided by T (return period)
# here 2 events per year -> 0.5 year separating two events in average => muByT = 0.5/T

# prepare array
matIsoRP_U_K = matIsoRP_V_K = matIsoRP_U_Kbar = matIsoRP_V_Kbar = 
  matIsoRP_Q_K = matIsoRP_D_K = matIsoRP_Q_Kbar = matIsoRP_D_Kbar = array(dim=c(nS,nRP,nUV))

for(iSt in 1:nS){
  cat(paste0(iSt,'-'))
  
  # retrieve copula
  cop = listCop[[iSt]]
  
  #___________________________________________________________#
  #                     Using K(t)
  #___________________________________________________________#
  
  for(iRP in 1:nRP){
    # mean elapsed time between two events (0.5 year) divided by T years (return period)
    muByT = 0.5/vecRP[iRP]
    
    # find t corresponding to this pk(t)
    t = seqt[which.min(abs(muByT-Pk[iSt,]))]
    
    # for different values of u, find v such that C(u,v)=t
    seqvcond = vector(length=nUV)
    
    # loop over u
    for(i in 1:nUV){
      if(all(BiCopCDF(rep(sequ[i],nUV), seqv, cop)-t<0)){
        seqvcond[i]=NA
      }else{
        # for this u, find v such that C(u,v)=t
        iMin = which.min(abs(BiCopCDF(rep(sequ[i],nUV), seqv, cop)-t))
        seqvcond[i]=seqv[iMin]
      }
    }
    
    matIsoRP_U_K[iSt,iRP,] = sequ
    matIsoRP_V_K[iSt,iRP,] = seqvcond
  }
  
  
  #___________________________________________________________#
  #                     Using Kbar(t)
  #___________________________________________________________#
  
  for(iRP in 1:nRP){
    # mean elapsed time between two events (0.5 year) divided by T years (return period)
    muByT = 0.5/vecRP[iRP]
    # find t corresponding to this pk(t)
    t = seqt[which.min(abs(muByT-PkBar[iSt,]))]
    
    # for different values of u, find v such that C(u,v)=t
    seqvcond = vector(length=nUV)
    
    # loop over u
    for(i in 1:nUV){
      if(all(BiCopCDF(rep(sequ[i],nUV), seqv, cop)-t<0)){
        seqvcond[i]=NA
      }else{
        # copula family
        if(cop$family==5){
          famSurv = 5
        }else{
          famSurv = cop$family+10
        }
        # survival cdf: Eq. 17 in Salvadori et al. 2014
        Surv.cdf = BiCopCDF(u1=rep(1-sequ[i],nUV),u2=1-seqv,family=famSurv, par=cop$par, par2=cop$par2)
        # for this u, find v such that Cbar(1-u,1-v)=t
        iMin = which.min(abs(Surv.cdf-t))
        seqvcond[i]=seqv[iMin]
      }
    }
    
    matIsoRP_U_Kbar[iSt,iRP,] = sequ
    matIsoRP_V_Kbar[iSt,iRP,] = seqvcond
  }
  
  
  #======================================================================
  # BACK TO THE DIMENSIONS OF THE VARIABLES
  #======================================================================
  
  # marginal distributions
  fitQ = listFitQ[[iSt]]
  fitD = listFitD[[iSt]]
  
  for(iRP in 1:nRP){
    ######################### K ###################################
    
    #_________________ Qmax ____________________
    seq.tmp =  matIsoRP_U_K[iSt,iRP,]
    notNA = !is.na(seq.tmp)&(!seq.tmp%in%c(0,1))
    matIsoRP_Q_K[iSt,iRP,notNA] = qgpd(seq.tmp[notNA], loc=fitQ$threshold, scale=fitQ$mle[1], shape=fitQ$mle[2])
    
    #_________________ Dequiv ____________________
    seq.tmp =  matIsoRP_V_K[iSt,iRP,]
    notNA = !is.na(seq.tmp)&(!seq.tmp%in%c(0,1))
    matIsoRP_D_K[iSt,iRP,notNA] = qgamma(seq.tmp[notNA],shape = fitD$estimate[1], rate = fitD$estimate[2])
    
    
    
    ######################### Kbar ###################################
    
    #_________________ Qmax ____________________
    matIsoRP_Q_Kbar[iSt,iRP,] = matIsoRP_Q_K[iSt,iRP,]
    
    #_________________ Dequiv ____________________
    seq.tmp =  matIsoRP_V_Kbar[iSt,iRP,]
    notNA = !is.na(seq.tmp)&(!seq.tmp%in%c(0,1))
    matIsoRP_D_Kbar[iSt,iRP,notNA] = qgamma(seq.tmp[notNA],shape = fitD$estimate[1], rate = fitD$estimate[2])
    
  }
}

# save outputs
save(matIsoRP_U_K,matIsoRP_V_K,matIsoRP_Q_K,matIsoRP_D_K,
     matIsoRP_U_Kbar,matIsoRP_V_Kbar,matIsoRP_Q_Kbar,matIsoRP_D_Kbar,
     file=paste0(path.datacrues,"getIsoReturnPeriods.RData"))

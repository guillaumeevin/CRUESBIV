###=====================================================================================================
### Guillaume Evin
### 16/03/2020, Grenoble
###  IRSTEA
### guillaume.evin@inrae.fr
### 
### Evin & Piton (Mars 2020): Codes développés dans le cadre de l'action SNRH CrueBiv - INRAE, UR ETNA
### "Analyse bivariée des liens entre magnitude et durée des crues en zones Alpine et Pyrénéenne"
###=====================================================================================================

###__________________________
### FitBivDist_ExtremeFloods.R
###__________________________
### Code ajustement des distributions bivariéees sur l'ensemble de crues S1 
### (128 stations présentant des chroniques >25 ans et au moins 2 crues par an en moyenne)


# REFERENCES

# Serinaldi, Francesco. 2015. "Dismissing Return Periods!" Stochastic Environmental Research and Risk Assessment
# 29 (4): 1179-89. https://doi.org/10.1007/s00477-014-0916-1.

#Salvadori, G., G. R. Tomasicchio, and F. D'Alessandro. 2014. "Practical Guidelines for Multivariate Analysis and 
#Design in Coastal and Off-Shore Engineering." Coastal Engineering 88 (June): 1-14. 
#https://doi.org/10.1016/j.coastaleng.2014.01.011.


# load libraries
library(fitdistrplus)
library(VineCopula)
library(lmomco)
library(ismev)
library(evd)


#======================================================================
# Load metadata for the set S_1
load(paste0(path.datacrues,"S1/Metadata.stations.S1.RData"))

# number of stations: 128
nS = nrow(Metadata.stations.S1)

# station codes
code.station = Metadata.stations.S1$Code


#======================================================================
# initialise objects
#======================================================================

# sequence of probabilities
seqt = seq(0.001, 1, 0.001)

# isoproba associée à un niveau de retour
Pk = PkBar = matrix(nrow=nS,ncol=length(seqt))

# list of data used for the fitting
listQ = listD = list()

# list of fitted copulas for each station
listCop = list()

# list of fitted distributions for peak flows and durations, for each station
listFitQ = listFitD = listQcdf = listDcdf = list()

vecQ10 = vector(length=nS)

#======================================================================
# for each station, we fit a bivariate distribution model (copula + margins)
#======================================================================

for(iSt in 1:nS){
  # load large floods (2 floods/year in average)
  load(paste0(path.datacrues,'S1/',code.station[iSt],'_DFcruesS1.RData'))
  
  # retrieve peak and equivalent duration
  Q.all = DFcruesS1$peak
  D.all = DFcruesS1$Deqt
  
  # remove nas
  nz = !is.na(Q.all)
  Q.raw = Q.all[nz]
  D = D.all[nz]
  listD[[iSt]] = D
  
  # first fit a gpd distribution on peaks : we have selected 2 floods per year (npy=2)
  # the peak associated to a return period = 10 years is used to adimensionalized the peaks
  #fitQall = gpd.fit(xdat=Q.raw, threshold = min(Q.raw),show=FALSE,npy=2)
  #Qadim = qgpd(p=0.5/10, loc=fitQall$threshold, scale=fitQall$mle[1], shape=fitQall$mle[2])
  Q = Q.raw #/Qadim
  listQ[[iSt]] = Q
  
  # fit gumbel distribution on peaks : we have selected 2 floods per year (npy=2)
  fitQ = gpd.fit(xdat=Q, threshold = min(Q),show=FALSE,npy=2)
  listFitQ[[iSt]] = fitQ
  vecQ10[iSt] = qgpd(p=1-0.5/10, loc=fitQ$threshold, scale=fitQ$mle[1], shape=fitQ$mle[2])
  
  # fit gamma distribution on durations
  fitD = fitdist(D, "gamma")
  listFitD[[iSt]] = fitD
  
  # cdf values
  Q.cdf = pgpd(Q, loc=fitQ$threshold, scale=fitQ$mle[1], shape=fitQ$mle[2])
  D.cdf = pgamma(D,shape = fitD$estimate[1], rate = fitD$estimate[2])
  listQcdf[[iSt]] = Q.cdf
  listDcdf[[iSt]] = D.cdf
  
  # Copula selection among archimedean copulas (for ease of computation)
  # The method of selection is based on the AIC criteria
  cop <- BiCopSelect(Q.cdf, D.cdf, familyset = 3:10,rotations = FALSE)
  listCop[[iSt]] = cop
  
  # scatterplot avec simulation
  U = BiCopSim(N=1000,obj = cop)
  X = U
  X[,1] = qgpd(p = U[,1], loc=fitQ$threshold, scale=fitQ$mle[1], shape=fitQ$mle[2])
  X[,2] = qgamma(p = U[,2],shape = fitD$estimate[1], rate = fitD$estimate[2])
  
  # Serinaldi (2015) Eq 15: Pk = 1-Kc(t), avec Kc(t) = t-lambda(t), see also ?BiCopLambda
  # ->Pk(t) = 1-t+lambda(t), lambda est donné par la fonction BiCopLambda
  Pk[iSt,] = 1 - seqt + BiCopLambda(cop, PLOT = FALSE)$theoLambda
  
  # +10 to get the corresponding survival copula, except for the Frank copula for which the survival copula is the same
  # due to its symmetry (as for elliptical copulas)
  if(cop$family==5){
    famSurv = 5
  }else{
    famSurv = cop$family+10
  }
  Ubar = 1-BiCopSim(N=10000,family=famSurv, par=cop$par, par2=cop$par2) 
  V = BiCopCDF(u1=Ubar[,1],u2=Ubar[,2],family=famSurv, par=cop$par, par2=cop$par2)
  PkBar[iSt,] = sapply(seqt, function(x) mean(V<x))
}

# save all objects in a RData
save(listQ, listD, listCop, listQcdf, listDcdf, listFitQ, listFitD, Pk, PkBar, seqt, vecQ10,
     file=paste0(path.datacrues,"fitBivDist.RData"))
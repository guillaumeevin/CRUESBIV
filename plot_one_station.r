###=====================================================================================================
### Guillaume Evin
### 16/03/2020, Grenoble
###  IRSTEA
### guillaume.evin@inrae.fr
### 
### Evin & Piton (Mars 2020): Codes d?velopp?s dans le cadre de l'action SNRH CrueBiv - INRAE, UR ETNA
### "Analyse bivari?e des liens entre magnitude et dur?e des crues en zones Alpine et Pyr?n?enne"
###=====================================================================================================


###__________________
### plot_one_station.r
###__________________


# produce various plots for one station


#======================================================================
# LOAD RESULTS
#======================================================================
# - Metadata.stations.S1: data.frame avec caracteristiques stations
load(paste0(path.datacrues,"S1/Metadata.stations.S1.RData"))

# - listCop,listQcdf,listDcdf,listFitQ,listFitD,PkBar,PkBar,seqt
load(paste0(path.datacrues,"fitBivDist.RData"))

# - listQ, listD, listCop, listQcdf, listDcdf, listFitQ, listFitD, Pk, PkBar, seqt, vecQ10: liste contenant les ajustements
# bivaries et les distributions marginales pour chaque station
# - matIsoRP_U_K,matIsoRP_V_K,matIsoRP_Q_K,matIsoRP_D_K,
# matIsoRP_U_Kbar,matIsoRP_V_Kbar,matIsoRP_Q_Kbar,matIsoRP_D_Kbar: matrices nStation x nPeriodeRetour x nPoints des contours pour
# une meme periode de retour. U et V dans l'espace des proba, Q et D dans l'espace des pics/durees, K et Kbar pour les deux definitions
load(paste0(path.datacrues,"getIsoReturnPeriods.RData"))

# metadata stations
nS = nrow(Metadata.stations.S1)
code.station = Metadata.stations.S1$Code

# return periods
vecRP = c(2,5,10)
nRP = 3

# index station
iSt = which(code.station=="O0384010")

# header
tagStation = paste0('Station ',code.station[iSt],'\nSuperficie = ',Metadata.stations.S1$Superficie[iSt],' km2')

#======================================================================
# FIGURE 1: SCATTERPLOT Q/D
#======================================================================
# Q and D observes
Q = listQ[[iSt]]
D = listD[[iSt]]

plot(D,Q,pch=20,cex=1.5,main=tagStation,
     xlab="Duree equivalente (heures)",ylab="Pics de crue (m3/s)",cex.lab=1.3)



#======================================================================
# FIGURE 2: AJUSTEMENT GPD SUR LES DEBITS MAX
#======================================================================
# fit gumbel distribution on peaks : we have selected 2 floods per year (npy=2)
fitQ = listFitQ[[iSt]]

# note gpd.diag fails when the shape parameter < 0, for the frequency plot
if(fitQ$mle[2]>0) gpd.diag(fitQ)


#======================================================================
# FIGURE 2: AJUSTEMENT GAMMA SUR LES DUREES
#======================================================================
fitD = listFitD[[iSt]]
plot(fitD)


#======================================================================
# FIGURE 3: AJUSTEMENT DE LA COPULE
#======================================================================
cop = listCop[[iSt]]
print(plot(cop))


#======================================================================
# FIGURE 4: SIMULATIONS BIVARIEES
#======================================================================
# Q and D observes
Q = listQ[[iSt]]
D = listD[[iSt]]

# ajustements pour cette station
cop = listCop[[iSt]]
fitQ = listFitQ[[iSt]]
fitD = listFitD[[iSt]]

# scatterplot avec simulation
U = BiCopSim(N=1000,obj = cop)
X = U
X[,1] = qgpd(p = U[,1], loc=fitQ$threshold, scale=fitQ$mle[1], shape=fitQ$mle[2])
X[,2] = qgamma(p = U[,2],shape = fitD$estimate[1], rate = fitD$estimate[2])

# scatterplot avec simulation
plot(X[,2],X[,1],pch=20,cex=1.5,col="gray",main=tagStation,
     xlab="Duree equivalente (heures)",ylab="Pics de crue (m3/s)",cex.lab=1.3)
points(D,Q,pch=20,cex=1.5,col="black")

#======================================================================
# FIGURE 5: ISOCOURBES DE MEME PERIODES DE RETOUR
#======================================================================

# type de lignes
vec.lty = c(3,2,1)

# Q and D
Q = listQ[[iSt]]
D = listD[[iSt]]

plot(D,Q,pch=20,cex=1.5,xlim=range(D),ylim=range(Q),frame.plot = FALSE,main=tagStation,
     xlab="Duree equivalente (heures)",ylab="Pics de crue (m3/s)",cex.lab=1.3)
for(iRP in 1:nRP){
  # add line
  x = matIsoRP_D_Kbar[iSt,iRP,]
  y = matIsoRP_Q_Kbar[iSt,iRP,]
  
  # add line
  nz = (x!=0 & y!=0) & (!is.na(x) & !is.na(y)) 
  x = x[nz]
  y = y[nz]
  lines(x,y,lty=vec.lty[iRP])
}
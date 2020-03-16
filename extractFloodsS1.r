###=====================================================================================================
### Guillaume Evin
### 16/03/2020, Grenoble
###  IRSTEA
### guillaume.evin@inrae.fr
### 
### Evin & Piton (Mars 2020): Codes développés dans le cadre de l'action SNRH CrueBiv - INRAE, UR ETNA
### "Analyse bivariée des liens entre magnitude et durée des crues en zones Alpine et Pyrénéenne"
###=====================================================================================================

###__________________
### extractFloodsS1.r
###__________________
### extract a second set of floods and stations

#================================================================================
# load Metadata for all stations
load(paste0(path.QtVarData,"Metadata.stations.RData"))

# number of stations
nS = nrow(Metadata.stations)

# stations code
code.stations = Metadata.stations$Code


#======================================================================
# number of floods and mean number of floods per year
vecLength = vecPeriodCovered = vector(length=nS)
for(iSt in 1:nS){
  # long RData
  load(paste0('./DGPR_CruesBivariees/DATA/DATACRUES/S0/',code.station[iSt],'/DFcrues.RData'))
  # number of floods: one row per flood in df.crues
  vecLength[iSt] = nrow(df.crues)
  vecPeriodCovered[iSt] = as.numeric(diff(range(as.Date(df.crues$td))))/365.25
}
NbEventperYear = vecLength/vecPeriodCovered


#======================================================================
# Selection of stations in S_1: at least 2 events per year and 25 years of data
isInS1 = which(NbEventperYear>2&vecPeriodCovered>25)
nS1 = length(isInS1)

for(iS1 in 1:nS1){
  # load infos et floods S0
  iSt = isInS1[iS1]
  superficie = Metadata.stations$Superficie[iSt]
  load(paste0(path.datacrues,'/S0/',code.station[iSt],'/DFcrues.RData'))
    
  # retrieve peak and equivalent duration
  Qraw = df.crues$peak
  Draw = df.crues$Deqt
    
  # selection of floods with the highest peaks iot retain 2 events / year
  nbFloods2retain = floor(vecPeriodCovered[iSt]*2)
  nbFloodsRaw = length(Qraw)
  rank_Q = rank(Qraw)
  sel = rank_Q > (nbFloodsRaw-nbFloods2retain)
  
  # filter duration > 30 days (should not happen)
  isDnottooLong = which(df.crues$D<(30*24))
  iSel1 = isDnottooLong[sel]

  # save selection in a data.frame  
  DFcruesS1 = df.crues[iSel1,]
  save(DFcruesS1,file=paste0(path.datacrues,'/S1/',code.station[iSt],'_DFcruesS1.RData'))
}

# save Metadata for S_1
Metadata.stations.S1 = Metadata.stations[isInS1,]
save(Metadata.stations.S1,file=paste0(path.datacrues,'S1/Metadata.stations.S1.RData'))
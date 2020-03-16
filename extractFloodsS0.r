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
### extractFloodsS0.r
###__________________
### CODE de découpage des séries temporelles Qtvar (extraction banque Hydro - logiciel HYDRO2), basé sur une
### première version développée par G. Piton (guillaume.piton@inrae.fr)
### Dans le repertoire .../S0, crée des sous-repertoires pour chaque station où il
### range les hydrogrammes découpés selon un critère de dépassement de seuil.

# load libraries
library(zoo)
library(xts)
library(date)

#================================================================================
# Function to extract characteristics from flood:
# 1. Interpolate starting and end time
# 2. Extract flood peak and corresponding time
# 3. Compute a volume (simple trapezoidal integration)
# 4. Compute an equivalent duration (volume divided by the peak)
#
# INPUTS
# - dfOneFlood: data.frame with columns "time" (POSIXct) and "q" (flood values, numeric), corresponding to a 
# period above the threshold (with bounds below the threshold)
# - Qseuil: threshold defining flood periods
#
# OUTPUT
# a list with:
# - times: 3 times (POSIXct): td (start), tf (end), tpeak (peak)
# - carac: 3 caracteristics: peak (m3/s), d (raw duration, in hours) ,v (volume, in hm3) ,deqv (equivalent duration, in hours)
extract.flood = function(dfOneFlood,Qseuil){
  # number of recordss
  nC = nrow(dfOneFlood)
  
  ####### starting time: interpolation using the distance from the second element
  # retrieve flow and time for elements 1 and 2
  q1 = dfOneFlood$q[1]
  q2 = dfOneFlood$q[2]
  t1 = dfOneFlood$time[1]
  t2 = dfOneFlood$time[2]
  # simple linear interpolation for time, assume that flow at this time corr to the thresold
  pct = (Qseuil-q1)/(q2-q1)
  dfOneFlood$time[1] = t1 + (t2-t1)*pct
  dfOneFlood$q[1] = Qseuil
  
  ####### end time: interpolation using the distance from the n-1 element
  # retrieve flow and time for elements n-1 and n
  qBefLast = dfOneFlood$q[nC-1]
  qLast = dfOneFlood$q[nC]
  tBefLast = dfOneFlood$time[nC-1]
  tLast = dfOneFlood$time[nC]
  # simple linear interpolation for time, assume that flow at this time corr to the thresold
  pct = (qBefLast-Qseuil)/(qBefLast-qLast)
  dfOneFlood$time[nC] = tBefLast + (tBefLast-tBefLast)*pct
  dfOneFlood$q[nC] = Qseuil
  
  ####### flood peak and time of peak: simply the maximum flow recorded and corr. time
  peak = max(dfOneFlood$q)
  i.peak = max(which.max(dfOneFlood$q))
  tpeak = dfOneFlood$time[i.peak]
  
  ####### duration in hours
  td = dfOneFlood$time[1]
  tf = dfOneFlood$time[nC]
  d = as.numeric(difftime(tf,td,units="hours"))
  if(is.na(d)){
    print(dfOneFlood)
    stop('extract.flood: cannot compute the difference of time')
  }
  
  ####### volume: trapezoidal integration
  q.moy = (dfOneFlood$q[2:nC] + dfOneFlood$q[1:(nC-1)])/2
  
  # diff in seconds
  t.diff = as.numeric(difftime(dfOneFlood$time[2:nC],dfOneFlood$time[1:(nC-1)],units="secs"))
  
  # volume in hm3
  v = sum(q.moy*t.diff)/100^3
  
  # equivalent duration in hours ( convert volume into m3, and convert secs to hours)
  deqv = 2*(v*100^3)/peak/3600
  
  return(list(time = c(td,tf,tpeak), carac = c(peak,d,v,deqv)))
}

#================================================================================
# load Metadata for all stations
load(paste0(path.QtVarData,"Metadata.stations.RData"))

# number of stations
nS = nrow(Metadata.stations)

# stations code
code.stations = Metadata.stations$Code

# number of values first kept around the peak (must be large)
# This number is relative to the temporal resolution (can be minutes or hours)
Nint<-50

# Initialize matrix of codes concerning floods which are kept or discarded
# - col. 1: nb of floods conserved
# - col. 2: nb of floods discarded with more than 12 h between successive measurements
# - col. 3: nb of floods discarded with an issue with the dates (negative equivalent duration)
# - col. 4: nb of floods discarded because too many points equal to Qmax
code.crues = matrix(0,nrow=nS,ncol=4)

# loop over all stations
for(iS in 1:nS){
  # initialise vectors
  # vec.index.xx indicate the position in the vector of raw data (index)
  # vec.xx is a vector of string indicating times (YYYY-MM-DD-HH-MM)
  # td indicate the start, tf the end and tpeak is the peaking time
  # peak (m3/s), d (raw duration, in hours) ,v (volume, in hm3) ,deqv (equivalent duration, in hours)
  vec.index.td = vec.td = vec.index.tf = vec.tf = vec.index.peak = vec.tpeak = vec.peak = 
    vec.D = vec.V = vec.Deqt = rep(NA,10000)
  vec.merge = rep(F,10000)
  
  # station code
  code = code.stations[iS]
  
  # create folder
  path.station = paste0(path.datacrues,"/S0/",code)
  if(dir.exists(path=path.station)==FALSE)
  {
    dir.create(path=path.station,showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  #_____________ Dstar:  caracteristic duration
  # Cipriani, Thomas, Tristan Toilliez, et Eric Sauquet. 2012. "Estimation régionale des débits décennaux et 
  # durées caractéristiques de crue en France". La Houille Blanche, n 4-5 (octobre): 5-13.
  # https://doi.org/10.1051/lhb/2012024.
  superficie = Metadata.stations$Superficie[iS]
  # Eq 3 in Cipriani et al. (2012): caracteristic duration in hours
  Dstar = 7.45*superficie^0.27
  
  #_____________ Read data  
  x<-scan(paste0(path.QtVarData,code,".txt"),what=character(0),sep=";",skip=4)       # Read data
  plage<-which(x=="920")                                                             # Floods records
  vec.Q.all <- as.numeric(x[plage+5])                                                # varying temporal resolution (QtVar)
  # QtVar times
  vec.time.all <- as.POSIXct(paste0(substr(x[plage+3],1,4),"-",substr(x[plage+3],5,6),'-',
                                    substr(x[plage+3],7,8),' ',x[plage+4]),tz="GMT",
                             format="%Y-%m-%d %H:%M")
  
  # lenght (time) of each interval
  vec.length.per = as.numeric(diff(vec.time.all), units = "hours")
  
  # remove very high or negative values for the computation of the module
  vec.length.per[vec.length.per>100] = NA
  vec.length.per[vec.length.per<0] = NA
  n.per = length(vec.length.per)
  
  # Module = mean flow
  module = sum(vec.Q.all[1:n.per]*vec.length.per,na.rm=T)/sum(vec.length.per,na.rm=T)
  
  # Threshold = 4 X module
  Qseuil = 4*module
  
  ###### find periods exceeding a threshold, use "rle" function, efficient way to retrieve this type of sequences
  # Exceed threshold
  rle.Q = rle(vec.Q.all>Qseuil)
  
  # Corresponding indices
  rle.index = c(1,cumsum(rle.Q$lengths))
  does.exceed = which(rle.Q$values)
  
  # index of starting and ending periods
  vec.index.td.all = rle.index[does.exceed]
  vec.index.tf.all = rle.index[does.exceed+1]+1
  
  # number of floods
  n.crue.all = length(does.exceed)
  
  # index of the end of the last flood
  vec.index.tf.all[n.crue.all] = min(vec.index.tf.all[n.crue.all],length(vec.Q.all))
  
  # count floods
  cnt.crue = 0
  
  ###### Retrieve properties of each flood, merge successive floods if floods oare close
  # (= time between end preceeding floods and start next floods < caracteristic duration)
  for(i.crue in 1:n.crue.all){
    # index start and end of flood
    index.td = vec.index.td.all[i.crue]
    index.tf = vec.index.tf.all[i.crue]
    dd = index.td:index.tf
    dd.ext = (index.td-Nint):(index.tf+Nint)
    
    # times between successive measurements
    if(min(dd)>2){
      time.elapsed.measurements = difftime(vec.time.all[dd], vec.time.all[dd-1], units = "days")
    }else{
      time.elapsed.measurements = difftime(vec.time.all[dd+1], vec.time.all[dd], units = "days")
    }
    
    # if there is more than 0.5 days between two successive measurements, we do not store the flood characteristics
    if((index.tf-index.td>1 & index.td>1) & max(time.elapsed.measurements)<0.5){
      Q.crue = vec.Q.all[dd]
      time.crue = vec.time.all[dd]
      equal.max = Q.crue==max(Q.crue)
      # add condition iot avoid many successive measurements with equal values
      if(length(unique(time.crue[equal.max]))<5){
        # merge if difference between peak time of the previous flood and starting time of the
        # current flood are less than Dstar (caracteristic time)
        if(cnt.crue>0){
          if(difftime(vec.time.all[index.td],vec.time.all[vec.index.peak[cnt.crue]],units = 'hours')<Dstar){
            do.merge= T
          }else{
            do.merge = F
          }
        }else{
          do.merge = F
        }
        
        if(do.merge){
          # if merge, just update the end of the flood and re-extract flood carac with the longest period
          vec.index.tf[cnt.crue] = index.tf
          
          # df pour la crue seulement
          dd.merge = vec.index.td[cnt.crue]:index.tf
          vec.index.peak[cnt.crue] = dd.merge[which.max(vec.Q.all[dd.merge])[1]]
          df.c = data.frame(time=vec.time.all[dd.merge],q=vec.Q.all[dd.merge]) 
          flood.out = extract.flood(df.c,Qseuil)
          
          vec.merge[cnt.crue] = T
        }else{
          # increment counter
          cnt.crue = cnt.crue + 1
          
          # record index
          vec.index.td[cnt.crue] = index.td
          vec.index.tf[cnt.crue] = index.tf
          vec.index.peak[cnt.crue] = dd[which.max(vec.Q.all[dd])[1]]
          
          # df pour la crue seulement
          df.c = data.frame(time=vec.time.all[dd],q=vec.Q.all[dd]) 
          flood.out = extract.flood(df.c,Qseuil) 
        }
        
        # fill vectors
        vec.td[cnt.crue] = format(flood.out$time[1],"%Y-%m-%d %H%-%M-%S")
        vec.tf[cnt.crue] = format(flood.out$time[2],"%Y-%m-%d %H%-%M-%S")
        vec.tpeak[cnt.crue] = format(flood.out$time[3],"%Y-%m-%d %H%-%M-%S")
        
        vec.peak[cnt.crue] = flood.out$carac[1]
        vec.D[cnt.crue] = flood.out$carac[2]
        vec.V[cnt.crue] = flood.out$carac[3]
        vec.Deqt[cnt.crue] = flood.out$carac[4]
      }else{
        # code : colonne 4, on incremente le nb de crues ecartees car au moins 3 points avec des dates differentes sont egaux au max de debit
        code.crues[iS,4] = code.crues[iS,4]+1
      }
    }else{
      # code : colonne 2, on incremente le nb de crues ecartees car + de 12h entre deux mesures
      code.crues[iS,2] = code.crues[iS,2]+1
    }
  }
  
  # create data.frame with all floods
  df.crues.all = data.frame(index.td = vec.index.td,
                            index.tf = vec.index.tf,
                            index.peak = vec.index.peak,
                            td = vec.td,
                            tf = vec.tf,
                            tpeak = vec.tpeak,
                            peak = vec.peak,
                            D = vec.D,
                            V = vec.V,
                            Deqt = vec.Deqt)
  
  
  #======================================
  # Filter floods with a date issue
  hasDatePb = which(vec.Deqt[1:cnt.crue]<=0|vec.D[1:cnt.crue]<=0)
  
  # code : colonne 1: nb crues conservees, colonne 3: nb crues avec duree equivalente < 0 (pb dates)
  code.crues[iS,1] = cnt.crue-length(hasDatePb)
  code.crues[iS,3] = length(hasDatePb)
  
  # filter negative duration: problems with dates or at the boundaries
  toKeep = which(vec.Deqt[1:cnt.crue]>0&vec.D[1:cnt.crue]>0)
  
  # trim data.frame to retain these events
  df.crues = df.crues.all[toKeep,]
  vec.merge = vec.merge[toKeep]
  
  # remove duplicated lines (should not happen since we merge close events)
  df.crues = df.crues[!duplicated(df.crues),]
  
  #======================= RData file ===================
  save(df.crues,file = paste0(path.station,"/DFcrues.RData"))
  
  
} #fin de la boucle sur toutes les stations

# compte crues ecartees et conservees
# - colonne 1: nb crues gardees
# - colonne 2: nb crues ecartees avec plus de 12h entre deux points de mesure
# - colonne 3: nb crues ecartees avec probleme dans les dates (duree equivalente negative)
# - colonne 4: nb crues ecartees parce que beaucoup de points egaux a Qmax
code.crues = as.data.frame(code.crues)
colnames(code.crues) = c("Retenues","plus_12h_entre_2_mesures","probleme_dates","au_moins_5_points_egaux_a_Qmax")
save(code.crues,file = paste0(path.datacrues,"code.crues.RData"))

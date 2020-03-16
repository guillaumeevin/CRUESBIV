###=====================================================================================================
### Guillaume Evin
### 16/03/2020, Grenoble
###  IRSTEA
### guillaume.evin@inrae.fr
### 
### Evin & Piton (Mars 2020): Codes développés dans le cadre de l'action SNRH CrueBiv - INRAE, UR ETNA
### "Analyse bivariée des liens entre magnitude et durée des crues en zones Alpine et Pyrénéenne"
###=====================================================================================================

#============================================
#    !!!!!! TO MODIFY !!!!!!!!
# 
#         PATH TO FOLDERS
#============================================

# path to a folder containing QtVar data for all stations, i.e. a list of .txt files named xxx.txt where xxx is the code of the station.
# This folder should also contain "Metadata.stations.RData" which contains a data.frame with the characteristics of the stations:
# 'Code','Nom','QSeuil','zone','Departement','SuperficieReelle','SuperficieTopo','Superficie','ErreurRelSuperficie','Perimetre',
# 'Imin','Imax','Imoy','Zmin','Zmax','Zmoy','Zsd','Lmax','Lequiv','Lrelat','Pente','Q2','Q5','Q10','Q20','Q50','DebitClasse',
# 'X','Y','Altitude','Debut','Fin'
path.QtVarData = "./DGPR_CruesBivariees/DATA/"

# path to a folder which will contain extracted flood characteristics
path.datacrues = "./DGPR_CruesBivariees/DATA/DATACRUES/"

# we create subfolders for each class of floods
dir.create(path=paste0(path.datacrues,"S0/"),showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create(path=paste0(path.datacrues,"S1/"),showWarnings = TRUE, recursive = FALSE, mode = "0777")


#============================================
#         Treatments
#============================================

# Extract the first set of floods Q_0.
source('extractFloodsS0.r')

# Extract the second set of floods Q_1.
source('extractFloodsS1')

# Fit a bivariate distribution Peak/Duration to the larget floods of each station.
source('FitBivDist_ExtremeFloods.R')

# Compute isolignes corresponding to return periods of 2, 5, 10 and 20 years.
source('getIsoReturnPeriods_ExtremeFloods.R')

# Calcule les temps de concentration et durées caractéristiques et le transport solide par évènement et le long du contour décennal.
# isoCurveRP10\_BivUniv.r

# Calcule le transport solide selon les hypothèses du paragraphe \ref{Sect::Qs}.
# ComputeBedloadVolumeTriangularHydrograph.r 
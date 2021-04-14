#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

if (length(args)!=1) {
  stop("One argument must be supplied with name like ${LON}_${LAT}_noVPD.tsv", call.=FALSE)
} 



#making mtclim.mtc43 output from climate data extracted from netcdf using clt 


oldw <- getOption("warn")
options(warn = -1)
require(dplyr,quietly = TRUE, warn.conflicts = FALSE)
require(tidyr,quietly = TRUE, warn.conflicts = FALSE)
require(ggplot2,quietly = TRUE, warn.conflicts = FALSE)
require(stringr,quietly = TRUE, warn.conflicts = FALSE)
require(readr,quietly = TRUE, warn.conflicts = FALSE)
options(warn = oldw)

#get stuff from args
#test
#args<-c('199_25_noVPD.tsv')
#argument should be
#${GLAC}_${CO2}_${O2}_${LON}_${LAT}_noVPD.tsv

glac_co2_o2_lon_lat_year<-str_split_fixed(args[1],"noVPD",2)[1]


glac<-str_split_fixed(glac_co2_o2_lon_lat_year,"_",6)[1]
co2<-str_split_fixed(glac_co2_o2_lon_lat_year,"_",6)[2]
o2<-str_split_fixed(glac_co2_o2_lon_lat_year,"_",6)[3]
lon<-str_split_fixed(glac_co2_o2_lon_lat_year,"_",2)[4]
lat<-str_split_fixed(glac_co2_o2_lon_lat_year,"_",2)[5]
year<-str_split_fixed(glac_co2_o2_lon_lat_year,"_",6)[6]

# print(glac)
# print(co2)
# print(o2)
# print(lon)
# print(lat)
#print(args)

cNames<-c("Yearday","Tmax (deg C)","Tmin (deg C)","Tday (deg C)","Prcp (cm)","Srad (W/m^2)","Q2M","Daylen (Sec)")

noVPD<-read_delim(file = args[1],delim="\t", col_names = cNames, col_types = cols())
#print(args[1])

if(!all(apply(noVPD, 2, is.numeric))){
  print(paste("Missing Data, making dummy mtc43 for",toString(args[1])))
  outDF2<-data.frame(Yearday=c(NA),`Tmax (deg C)`=c(NA),`Tmin (deg C)`=c(NA),
                          `Tday (deg C)`=c(NA),`Prcp (cm)`=c(NA),VPD=c(NA),`Srad (W/m^2)`=c(NA),
                          `Daylen (Sec)`=c(NA))
  
  ###for when i integrate glac,co2,o2 
  # perth<-paste("./",glac,"_",co2,"_",o2,"_",lon,"_",lat,"_",year,".mtc43",sep="")
  perth<-paste("./",glac_co2_o2_lon_lat_year,".mtc43",sep="")
  #perth<-paste("./GLAC_182_28_",lon,"_",lat,".mtc43",sep="")
  #print(perth)
  #control precision of each variable
  write_delim(outDF2, path=perth, delim="\t")
  quit()
}

# #T in K, convert to C
noVPD$`Tmin (deg C)`<-noVPD$`Tmin (deg C)`-273.5
noVPD$`Tmax (deg C)`<-noVPD$`Tmax (deg C)`-273.5
noVPD$`Tday (deg C)`<-noVPD$`Tday (deg C)`-273.5
# #(m/s) to cm/day: *86400 s/day * 100 cm/m
noVPD$`Prcp (cm)`<-noVPD$`Prcp (cm)`*86400*100


# 
# #calculate SVP
# #(pressure in Pa and temperature in degrees Celsius)
# #saturation vapor pressure | eg T= 27 C give P = 3567 Pa (pascals)
# 
sat_vap_pres_W <- function(temp){
   atm_pr = 611.21*exp(((temp*(18.678 - temp/234.5))/(257.14 + temp)))
  return(atm_pr)
}
# 
sat_vap_pres_I <- function(temp){
  atm_pr = 611.15*exp(((temp*(23.036 - temp/333.7))/(279.82 + temp)))
  return(atm_pr)
}



# #conditional assignment/calculation of daily VP based on 0<Tmin?
# # if you don't have a way to calculate vapor pressure directly assume that svp(tave)-svp(tmin) is vpd
noVPD$VP<-NA
noVPD$VP[which(noVPD$`Tmin (deg C)`<0)]<-sat_vap_pres_I(noVPD$`Tmin (deg C)`[which(noVPD$`Tmin (deg C)`<0)])
noVPD$VP[which(noVPD$`Tmin (deg C)`>=0)]<-sat_vap_pres_W(noVPD$`Tmin (deg C)`[which(noVPD$`Tmin (deg C)`>=0)])
# noVPD$VP<-sat_vap_pres_W(noVPD$`Tmin (deg C)`)


noVPD$SVP<-NA
noVPD$SVP[which(noVPD$`Tday (deg C)`<0)]<-sat_vap_pres_I(noVPD$`Tday (deg C)`[which(noVPD$`Tday (deg C)`<0)])
noVPD$SVP[which(noVPD$`Tday (deg C)`>=0)]<-sat_vap_pres_W(noVPD$`Tday (deg C)`[which(noVPD$`Tday (deg C)`>=0)])
# noVPD$SVP<-sat_vap_pres_W(noVPD$`Tday (deg C)`)
# #calc VPD
noVPD$VPD<-noVPD$SVP-noVPD$VP
# 
noVPD$Year<-100
#print(noVPD$)

# #make .mtc43 output
outDF<-noVPD%>%select(Yearday,`Tmax (deg C)`,`Tmin (deg C)`,`Tday (deg C)`,`Prcp (cm)`,VPD,`Srad (W/m^2)`, 'Daylen (Sec)')
#,`Daylen (Sec)`)
# #round everything to 2 decimal pl
outDF2<-round(outDF,2)
# #change yearday to int

#LAT_IGLAC_546_28_87_DAYLEN.txt 
# dlf_fn<-paste("../daylengths/LAT",glac,co2,o2,lat,"DAYLEN.txt",sep="_")
# dlf<-read_delim(file = dlf_fn ,delim="\t", col_names = FALSE, col_types = cols())

#print(dlf)
# outDF2$`Daylen (Sec)`<-dlf$X1

outDF2<-outDF2%>%mutate(Yearday=as.integer(Yearday),`Tmax (deg C)`=sprintf("%0.2f",`Tmax (deg C)`),`Tmin (deg C)`=sprintf("%0.2f",`Tmin (deg C)`),
    `Tday (deg C)`=sprintf("%0.2f",`Tday (deg C)`),`Prcp (cm)`=sprintf("%0.2f",`Prcp (cm)`),VPD=sprintf("%0.2f",VPD),`Srad (W/m^2)`=sprintf("%0.2f",`Srad (W/m^2)`),
    `Daylen (Sec)`=sprintf("%0.2f",as.numeric(`Daylen (Sec)`)))


# # plot(outDF2$`Prcp (cm)`)
# 
###for when i integrate glac,co2,o2 
# perth<-paste("./",glac,"_",co2,"_",o2,"_",lon,"_",lat,"_",year,".mtc43",sep="")
perth<-paste("./",glac_co2_o2_lon_lat_year,".mtc43",sep="")
#perth<-paste("./GLAC_182_28_",lon,"_",lat,".mtc43",sep="")
#print(perth)
#control precision of each variable
write_delim(outDF2, path=perth, delim="\t")
# ####### make precision 2 decimal places!!!!!
# 
# 
# 

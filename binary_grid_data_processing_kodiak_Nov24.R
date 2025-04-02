#!/usr/bin/env Rscript --vanilla
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
###constants, should be the same for all files
#data produced by running *-BGC using iniWriter, or StVar.sh 
header<-c("prcp","tmax","tmin","tavg","tday","tnight","tsoil","vpd","swavgfd","swabs","swtrans","swabs_per_plaisun","swabs_per_plaishade","ppfd_per_plaisun","ppfd_per_plaishade","par","parabs","pa","co2","dayl","epv.proj_lai","summary.daily_et","soil_psi","leaf_psi","gl_s_sun","gl_s_shade","leafc","vegc","leafn","vegn","psn_sun_A","psn_shade_A","gpp","npp","nep","gl_t_wv_sun","gl_t_wv_shade","soilw_outflow","soilw","livestemc","deadstemc","cpool_to_livestemc","cpool_to_deadstemc","m_psi_x")
mets<-c("prcp","tmax","tmin","tavg","tday","tnight","tsoil","vpd","swavgfd","swabs","swtrans","swabs_per_plaisun","swabs_per_plaishade","ppfd_per_plaisun","ppfd_per_plaishade","par","parabs","pa","co2","dayl")
outs<-c("epv.proj_lai","summary.daily_et","soil_psi","leaf_psi","gl_s_sun","gl_s_shade","leafc","vegc","leafn","vegn","psn_sun_A","psn_shade_A","gpp","npp","nep","gl_t_wv_sun","gl_t_wv_shade","soilw_outflow","soilw","livestemc","deadstemc","cpool_to_livestemc","cpool_to_deadstemc","m_psi_x")
varnum <- length(header)
step_length<-"day"

simyears<-50
#
simmonths<-simyears*12
simdays<-simyears*365

#locations with files
location<-"/path/to/your/dayouts/"
setwd(location)

steps<-case_when(step_length == "year" ~ simyears,
                 step_length == "month" ~ simmonths,
                 step_length == "day" ~ simdays)
records<-steps*varnum


allPlots = FALSE

files <- list.files(pattern = ".dayout", full.names = FALSE)
num_files <- length(files)
global <- data.frame(matrix(ncol = length(header), nrow = 0))
colnames(global) <- header

# .Machine$sizeof.long
# alternative 'tidy' method 
my_read_bin <- function(filename, in_header){
  #filename<- global$File[1]
  #in_header <- header
  in_file<-file(filename, "rb")
  float_bytes <- 4
  n_total_elements <- file.info(filename)$size/float_bytes
  
  my_data <- readBin(filename, 
                     what = "numeric", n=9E6, size = float_bytes, endian = "little")
  close(in_file)
  
  df<- as.data.frame(t(matrix(my_data, nrow = varnum)))
  rm(my_data)
  colnames(df)<-in_header

  df$simstep<-seq(1,steps)
  df$yearstep<- case_when(step_length == "day" ~ as.numeric(df$simstep%%365),
                          step_length == "month" ~ as.numeric(df$simstep%%12),
                          step_length == "year" ~ as.numeric(df$simstep))
  
  df_averages<-df %>% 
        summarise(dplyr::across(header,mean))%>%
        rename_with(~str_c(.,"_mean"))
  
  df_mins<-df %>% 
    summarise(dplyr::across(mets,min))%>%
    rename_with(~str_c(.,"_min"))
  
  df_maxs<-df %>% 
    summarise(dplyr::across(mets,max))%>%
    rename_with(~str_c(.,"_max"))

  df_sds<-df %>% 
    summarise(dplyr::across(header,sd))%>%
    rename_with(~str_c(.,"_sd"))
  
  df_meds<- df %>% 
    summarise(dplyr::across(header,median))%>%
    rename_with(~str_c(.,"_median"))
  
  rm(df)
  
  df_stats<- cbind(df_averages, df_mins, df_maxs, df_sds, df_meds)
  
  return(df_stats)
}



filehead<-c("Lon", "Lat", "Genus", "Empirical Model", "Fc")

##examples
#212.5_48.3158_arthropitys_Eud_Fc05_MEDCN_MEDg.annavgout
#2x_187.5_-42.6316_cordaites_Eud_Fc05_MEDCN_MEDg.annavgout 
part_regex<-"([-.0-9]{1,6})_([-.0-9]{1,10})_([a-zS]+)_([A-Za-z]{3})_Fc([0-9]{1,2})_MEDCN_MEDg\\."
this_regex <- case_when(step_length == "day" ~ paste0(part_regex,"dayout"),
                        step_length == "month" ~ paste0(part_regex,"monavgout"),
                        step_length == "year" ~ paste0(part_regex,"annavgout"))

global <- tibble(File = unlist(files)) %>%
  ##reverse order of filename chunks
  #filename chunks 190_10.4211_lepidodendron_SVP_Fc1_MEDCN_MEDg.annavgout
  extract(col = File, into = filehead,
          regex=
            this_regex
          , remove = FALSE)
  
  print("the following tibble should be empty (filename extraction)")
  global[which(is.na(global)),"File"]
  
global <- global %>%
  mutate(Data = lapply(File, my_read_bin, in_header = header)) %>%
  unnest(Data)

global<-global %>% mutate(Lon = as.numeric(Lon), Lat = as.numeric(Lat))

global <- global%>%
  mutate(F_c = case_when(Fc == "05" ~ 0.5,Fc == "1" ~ 1.0, Fc == "10" ~ 1.0))



saveRDS(object = global, file="allDayoutsSummary_CESM_stemBGC_24Sept2024.RDS")


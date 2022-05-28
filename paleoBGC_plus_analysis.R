#script to reproduce analysis and produce graphics for manuscript 
#"Stems Matter: Xylem Physiological Limits Are an Accessible and Critical 
#Improvement to Models of Plant Gas Exchange in Deep Time"
#Data DOI:10.5281/zenodo.6588761
#Author: William J. Matthaeus
#updated 05-27-2021
require(tidyr)
require(ggplot2)
require(dplyr)
require(tibble)
require(stringr)
require(readr)
require(mgcv)
#Daily
#collated with awk:
#stem
# awk '{print FILENAME" "FNR" "$1" "$7" "$8" "$9" "$10" "$11" "$12" "$13" "$20" "$21" "$22" "$25" "$26" "$27" "$28" "$29" "$30" "$31" "$33" "$34" "$35" "$36" "$38" "$41" "$45" "$53" "$54" "$55" "$56" "$57" "$58" "$59" "$23" "$24}' *.dayout.ascii > stem_v6_102621.txt
#base
#awk '{print FILENAME" "FNR" "$1" "$7" "$8" "$9" "$10" "$11" "$12" "$13" "$20" "$21" "$22" "$25" "$26" "$27" "$28" "$29" "$30" "$31" "$33" "$34" "$35" "$36" "$38" "$41" "$45" "$53" "$54" "$55" "$56" "$57" "$58" "$23" "$24}' *.dayout.ascii > base_v6_102621.txt

setwd("/Users/willmatthaeus/Dropbox/Brent_Spiner/stem_bgc_results/dayOuts")
DailyFile=TRUE
MonthlyFile=FALSE

#Monthly
#collated with awk:
#stem
# awk '{print FILENAME" "FNR" "$1" "$7" "$8" "$9" "$10" "$11" "$12" "$13" "$20" "$21" "$22" "$25" "$26" "$27" "$28" "$29" "$30" "$31" "$33" "$34" "$35" "$36" "$38" "$41" "$45" "$53" "$54" "$55" "$56" "$57" "$58" "$59" "$23" "$24}' stem*.monavgout.ascii > stem_monthly_Xbgc_060522.txt
#
# awk '{print FILENAME" "FNR" "$1" "$7" "$8" "$9" "$10" "$11" "$12" "$13" "$20" "$21" "$22" "$25" "$26" "$27" "$28" "$29" "$30" "$31" "$33" "$34" "$35" "$36" "$38" "$41" "$45" "$53" "$54" "$55" "$56" "$57" "$58" "$23" "$24}' paleo*.monavgout.ascii > base_monthly_Xbgc_060522.txt
setwd("/Users/willmatthaeus/Dropbox/Brent_Spiner/CESM2BGC/monthly_Xbgc_060522/")
DailyFile=FALSE
MonthlyFile=TRUE

# save.image(file = "stem_analysis_102721.RData")
# load("stem_analysis_102721.RData")
##
##data import, processing, calculations filtering
if(TRUE){
  #basic file import and df setup
  if(DailyFile){
    base_head<-c("filename","simday",
             "ws.soilw",
             "wf.soilw_outflow",
             "cs.leafc",
             "cs.frootc",
             "cs.livestemc",
             "cs.deadstemc",
             "cs.livecrootc",
             "cs.deadcrootc",
             "epv.proj_lai",
             "epv.psi",
             "epv.psi_leaf",
             "epv.gl_s_sun",
             "epv.gl_s_shade",
             "epv.m_Kl",
             "psn_sun.g",
             "psn_sun.Ci",
             "psn_sun.Ca",
             "psn_sun.A",
             "psn_shade.g",
             "psn_shade.Ci",
             "psn_shade.Ca",
             "psn_shade.A",
             "summary.daily_npp",
             "summary.daily_gpp",
             "summary.vegc",
             "metv.prcp",
             "metv.tmin",
             "psn.Av",
             "psn.Aj",
             "epv.gl_t_sun",
             "epv.gl_t_shade",
             "epv.daily_net_nmin",
             "nf.sminn_leached")

      stem_head<-c("filename","simday",
             "ws.soilw",
             "wf.soilw_outflow",
             "cs.leafc",
             "cs.frootc",
             "cs.livestemc",
             "cs.deadstemc",
             "cs.livecrootc",
             "cs.deadcrootc",
             "epv.proj_lai",
             "epv.psi",
             "epv.psi_leaf",
             "epv.gl_s_sun",
             "epv.gl_s_shade",
             "epv.m_Kl",
             "psn_sun.g",
             "psn_sun.Ci",
             "psn_sun.Ca",
             "psn_sun.A",
             "psn_shade.g",
             "psn_shade.Ci",
             "psn_shade.Ca",
             "psn_shade.A",
             "summary.daily_npp",
             "summary.daily_gpp",
             "summary.vegc",
             "m_psi_x",
             "metv.prcp",
             "metv.tmin",
             "psn.Av",
             "psn.Aj",
             "epv.gl_t_sun",
             "epv.gl_t_shade",
             "epv.daily_net_nmin",
             "nf.sminn_leached")
  } 
  if(MonthlyFile){
      base_head<-c("filename","simmonth",
                   "ws.soilw",
                   "wf.soilw_outflow",
                   "cs.leafc",
                   "cs.frootc",
                   "cs.livestemc",
                   "cs.deadstemc",
                   "cs.livecrootc",
                   "cs.deadcrootc",
                   "epv.proj_lai",
                   "epv.psi",
                   "epv.psi_leaf",
                   "epv.gl_s_sun",
                   "epv.gl_s_shade",
                   "epv.m_Kl",
                   "psn_sun.g",
                   "psn_sun.Ci",
                   "psn_sun.Ca",
                   "psn_sun.A",
                   "psn_shade.g",
                   "psn_shade.Ci",
                   "psn_shade.Ca",
                   "psn_shade.A",
                   "summary.daily_npp",
                   "summary.daily_gpp",
                   "summary.vegc",
                   "metv.prcp",
                   "metv.tmin",
                   "psn.Av",
                   "psn.Aj",
                   "epv.gl_t_sun",
                   "epv.gl_t_shade",
                   "epv.daily_net_nmin",
                   "nf.sminn_leached")
      
      stem_head<-c("filename","simmonth",
                   "ws.soilw",
                   "wf.soilw_outflow",
                   "cs.leafc",
                   "cs.frootc",
                   "cs.livestemc",
                   "cs.deadstemc",
                   "cs.livecrootc",
                   "cs.deadcrootc",
                   "epv.proj_lai",
                   "epv.psi",
                   "epv.psi_leaf",
                   "epv.gl_s_sun",
                   "epv.gl_s_shade",
                   "epv.m_Kl",
                   "psn_sun.g",
                   "psn_sun.Ci",
                   "psn_sun.Ca",
                   "psn_sun.A",
                   "psn_shade.g",
                   "psn_shade.Ci",
                   "psn_shade.Ca",
                   "psn_shade.A",
                   "summary.daily_npp",
                   "summary.daily_gpp",
                   "summary.vegc",
                   "m_psi_x",
                   "metv.prcp",
                   "metv.tmin",
                   "psn.Av",
                   "psn.Aj",
                   "epv.gl_t_sun",
                   "epv.gl_t_shade",
                   "epv.daily_net_nmin",
                   "nf.sminn_leached")
      }
  
if(TRUE){
#read in base model outputs
# Daily<-read.delim("base_v6_102621.txt",sep = " ", header = FALSE, stringsAsFactors = FALSE)
  Daily<-read.delim("base_monthly_Xbgc_060522.txt",sep = " ", header = FALSE, stringsAsFactors = FALSE)
colnames(Daily) <- base_head
drops<-c("simset","m_psi_x")
Daily <- Daily[,!(names(Daily) %in% drops)]
Daily$simset<-"base"
#so i can merge the two
Daily$m_psi_x<-1.0
#check to see there's not missing data
# which(is.na(Daily))
# which(Daily=="")
system("say Hot Coffee!")

stemDaily<-read.delim("stem_monthly_Xbgc_060522.txt",sep = " ", header = FALSE, stringsAsFactors = FALSE)
colnames(stemDaily)<-stem_head
stemDaily$simset<-"stem_leafsoil"
#check to see there's not missing data
# which(is.na(stemDaily))
# which(stemDaily=="")

system("say Hot Coffee!")



Daily<-rbind(Daily, stemDaily)

rm(stemDaily)
# Daily.bak<-Daily
# rm(Daily.bak)
#check to see there's not missing data
# which(is.na(Daily))
# colnames(Daily)

#do this first so you can filter by o2 (time saver)
#or just separate all if montly
if(DailyFile){
Daily$name<-str_split_fixed(Daily$filename,"[.]",3)[,1]
Daily$o2<-str_split_fixed(Daily$name,"[_]",8)[,3]
Daily <- Daily %>% filter(o2 != "21")
}else 
  if(MonthlyFile){
  #paleo_1x_280_120_2.84211_cordaites_MEDCN_MEDg_PRCP.monavgout.ascii
  Daily$name<-str_split_fixed(Daily$filename,".monavgout.",2)[,1]
  #take 'PRCP' out of paleo_iniWriter.sh!!!!!!!! then can remove this next
  Daily$name2<-str_split_fixed(Daily$name,"_MEDCN_",2)[,1]
  #
  temp<-str_split_fixed(Daily$name2,"[_]",6)
  Daily$simtype<-temp[,1]
  Daily$co2<-temp[,3]
  Daily$lon<-temp[,4]
  Daily$lat<-temp[,5]
  Daily$species<-temp[,6]
    #unique(Daily$CO2)
  Daily$sim<-paste(sep = "_",Daily$co2,Daily$species)
  rm(temp)
  }
#check to see there's not missing data
# which(is.na(Daily))
# which(Daily=="")
# system("say Hot Coffee!")


#check to see there's not missing data
# which(is.na(Daily))
# which(Daily=="")
# Daily.bak <- Daily
# Daily<-Daily.bak
if(DailyFile){
    Daily$species<-str_split_fixed(Daily$name,"[_]",8)[,5]
    Daily$simset[grep("Rev", Daily$species)]<-"stemRev"
    Daily$species<-str_split_fixed(Daily$species,"[R]",2)[,1]
    Daily <- Daily %>% filter(simset != "stemRev")
    
    # which(is.na(names))
    # which(names=="")
    # rm(names)
    
    Daily$glac<-str_split_fixed(Daily$name,"[_]",8)[,1]
    Daily$co2<-str_split_fixed(Daily$name,"[_]",8)[,2]
    Daily$location<-str_split_fixed(Daily$name,"[_]",8)[,4]
    
    Daily$CN<-str_split_fixed(Daily$name,"[_]",8)[,6]
    # gs<-str_split_fixed(Daily$name,"[_]",8)[,7]
    # which(gs!="MEDg")
    Daily$g<-str_split_fixed(Daily$name,"[_]",8)[,7]
    Daily$mprcp<-str_split_fixed(Daily$name,"[_]",8)[,8]
    #
    
    
    # which(is.na(Daily))
    # which(Daily=="")
    system("say Hot Coffee!")
    
    Daily <- Daily %>% dplyr::select(!c('filename','CN','g'))
    # colnames(Daily)
    # which(is.na(Daily))
    # which(Daily=="")
    # unique(Daily$mprcp)
    # Daily.bak<-Daily
    # Daily[which(Daily$mprcp=="0pt1PRCP"),]
    Daily$mprcp[which(Daily$mprcp=="0pt1PRCP")]<-"10%"
    Daily$mprcp[which(Daily$mprcp=="0pt2PRCP")]<-"20%"
    Daily$mprcp[which(Daily$mprcp=="0pt3PRCP")]<-"30%"
    Daily$mprcp[which(Daily$mprcp=="0pt4PRCP")]<-"40%"
    Daily$mprcp[which(Daily$mprcp=="0pt5PRCP")]<-"50%"
    Daily$mprcp[which(Daily$mprcp=="0pt6PRCP")]<-"60%"
    Daily$mprcp[which(Daily$mprcp=="0pt7PRCP")]<-"70%"
    Daily$mprcp[which(Daily$mprcp=="0pt8PRCP")]<-"80%"
    Daily$mprcp[which(Daily$mprcp=="0pt9PRCP")]<-"90%"
    Daily$mprcp[which(Daily$mprcp=="1pt0PRCP")]<-"100%"
    # which(is.na(Daily))
    # which(Daily=="")
    
    # system("say Hot Coffee!")
    
    # unique(Daily$mprcp)
    Daily$num.mprcp<-NA
    Daily$num.mprcp[which(Daily$mprcp=="10%")]<-"10"
    Daily$num.mprcp[which(Daily$mprcp=="20%")] <-"20"
    Daily$num.mprcp[which(Daily$mprcp=="30%")] <-"30"
    Daily$num.mprcp[which(Daily$mprcp=="40%")] <-"40"
    Daily$num.mprcp[which(Daily$mprcp=="50%")] <-"50"
    Daily$num.mprcp[which(Daily$mprcp=="60%")] <-"60"
    Daily$num.mprcp[which(Daily$mprcp=="70%")] <-"70"
    Daily$num.mprcp[which(Daily$mprcp=="80%")] <-"80"
    Daily$num.mprcp[which(Daily$mprcp=="90%")] <-"90"
    Daily$num.mprcp[which(Daily$mprcp=="100%")]<-"100"
    Daily$num.mprcp<-as.numeric(Daily$num.mprcp)
    # which(is.na(Daily))

# Daily.bak<-Daily



Daily$sim<-paste(sep = "_",Daily$co2,Daily$o2,Daily$species,Daily$mprcp)
Daily$climate<-paste(sep = "_",Daily$co2,Daily$o2,Daily$mprcp)
Daily$gcm<-paste(sep = "_",Daily$co2,Daily$o2)
}


# which(is.na(Daily))
# which(Daily=="")
}
  #calculate some new variables
  if(TRUE){
##calculate wuei https://doi.org/10.2135/cropsci2002.1220
# hist(Daily$epv.gl_s_shade)
# hist(Daily$epv.gl_s_shade[which(Daily$epv.gl_s_shade>0)])
# hist(Daily$epv.gl_s_sun)

##problems how to calculate where gl_s is close to zer
#have to convert gl_s to mol H20 (won't fix the problem)
# Daily$WUEi.sun<-Daily$psn_sun.A/Daily$epv.gl_s_sun
# Daily$WUEi.sun[which(Daily$epv.gl_s_sun<=0.01)]<-NA
# Daily$WUEi.shade<-Daily$psn_shade.A/Daily$epv.gl_s_shade
# Daily$WUEi.shade[which(Daily$epv.gl_s_shade<=0.01)]<-NA
#weighted average based on 2:1 ratio shade:sun canopy
# Daily$WUEi.t<- (2*Daily$WUEi.shade + Daily$WUEi.sun)/3
# hist(Daily$WUEi.t)  
##calcualte vwc: from bgc source
# /* convert kg/m2 --> m3/m2 --> m3/m3 */
#   vwc = soilw / (1000.0 * sitec->soil_depth);
    Daily$vwc <- Daily$ws.soilw/(1000 * 0.5)
    
    Daily$vegc<-Daily$cs.deadcrootc+Daily$cs.livecrootc+Daily$cs.deadstemc+
      Daily$cs.livestemc+Daily$cs.frootc+Daily$cs.leafc
    
    Daily$live.vegc<-Daily$cs.livecrootc+Daily$cs.livestemc+Daily$cs.frootc+Daily$cs.leafc
    Daily$dead.vegc<-Daily$cs.deadcrootc+Daily$cs.deadstemc
    
    Daily$rootc <- Daily$cs.livecrootc+Daily$cs.frootc+Daily$cs.deadcrootc
    Daily$stemc <- Daily$cs.deadstemc+Daily$cs.livestemc
    
    Daily$lsw_ratio <- (Daily$cs.livestemc/Daily$cs.leafc)/0.22
    hist(Daily$lsw_ratio, breaks = 100)
      
    Daily$WUEi.sun<-Daily$psn_sun.A/Daily$epv.gl_t_sun
    # hist(Daily$WUEp.sun)
    # unique(is.na(Daily$WUEi.sun))
    #looks good?
    Daily$WUEi.shade<-Daily$psn_shade.A/Daily$epv.gl_t_shade
    #weighted average based on 2:1 ratio shade:sun canopy
    Daily$WUEi.t<- (2*Daily$WUEi.shade + Daily$WUEi.sun)/3
    # hist(Daily$WUEi.t)
    #production based
    Daily$WUEe.sun <- Daily$vegc/Daily$epv.gl_t_sun
    # hist(Daily$WUEp.sun)
    # unique(is.na(Daily$WUEi.sun))
    #looks good?
    Daily$WUEe.shade <- Daily$vegc/Daily$epv.gl_t_shade
    #weighted average based on 2:1 ratio shade:sun canopy
    Daily$WUEe.t<- (2*Daily$WUEe.shade + Daily$WUEe.sun)/3
    


##total A
Daily$A.t <- Daily$psn_shade.A+Daily$psn_sun.A
#total gls
Daily$gls.t<- Daily$epv.gl_s_sun + Daily$epv.gl_s_shade
#total gl
Daily$gl.t <- Daily$epv.gl_t_shade + Daily$epv.gl_t_sun
#total gco2
Daily$gco2 <- Daily$psn_shade.g + Daily$psn_sun.g

if(DailyFile){
Daily$simyearday<-(Daily$simday%%365)+1
Daily$simyear<-ceiling(Daily$simday/50)
}else 
  if(MonthlyFile){
    Daily$simyearmonth<-(Daily$simmonth%%12)+1
    Daily$simyear<-ceiling(Daily$simmonth/50)
  }

Daily$A_type_diff<-Daily$psn.Av-Daily$psn.Aj
Daily$A_type<-NA
Daily$A_type[which(Daily$A_type_diff<0)]<-"Av"
Daily$A_type[which(Daily$A_type_diff>=0)]<-"Aj"
which(is.na(Daily$A_type))

Daily$species2<-NA
Daily$species2<-Daily$species
Daily$species2[grep("ebf", Daily$species2)]<-"ebf"
Daily$species2[grep("enf", Daily$species2)]<-"enf"
Daily$species2<-factor(Daily$species2, levels=c("lepidodendron","macroneuropteris","treefern","cordaites","ebf","enf"))

Daily$fossil<-"fossil"
Daily[which(Daily$species %in% c("ebf","enf","ebfArp","ebfMA")),]$fossil<-"modern"


Daily$species<-factor(Daily$species, levels = c("lepidodendron","macroneuropteris",
                                 "treefern","cordaites",
                                 "ebf","enf","ebfArp","ebfMA"))
# cbbPalette <- c("#E69F00", "#009E73", "#0072B2", "#56B4E9","#000000","#000000", "#D55E00", "#CC79A7")

#where moderns are grouped by orignial pft

# Daily$month<-Daily$simyearday*
# which(is.na(Daily))
# rm(Daily.bak)
# seq(1,30,1)%%30
}
  system("say -v Ava  Finished data import!")
}


#do some filtering in prep for graphics
if(TRUE){
  
  #summary data frames
  if(TRUE){
  Annualized <- Daily %>% 
    filter(num.mprcp==100, 
           species2 %in% c('cordaites', 'lepidodendron','macroneuropteris','treefern','ebf'),
           # species %in% c('cordaites', 'lepidodendron','macroneuropteris','treefern',
           #                'ebf','enf','ebfArp','ebfMA'),
           simset %in% c('base','stem_leafsoil')) %>%  
    group_by(simset,gcm,species,simyearday) %>% 
    summarize(species2 = unique(species2),
                fossil = unique(fossil),
                glst.m = mean(gls.t), glst.sd = sd(gls.t),
                At.m = mean(A.t), At.sd = sd(A.t),
                m_psi_x.m = mean(m_psi_x), m_psi_x.sd = sd(m_psi_x),
                lswr.m = mean(lsw_ratio), lswr.sd = sd(lsw_ratio),
                vegc.m = mean(vegc), vegc.sd = sd(vegc),
                nmin.m = mean(epv.daily_net_nmin), nmin.sd = sd(epv.daily_net_nmin),
                nlea.m = mean(nf.sminn_leached), nlea.sd = sd(nf.sminn_leached),
                live.vegc.m = mean(live.vegc), live.vegc.sd = sd(live.vegc),
                dead.vegc.m = mean(dead.vegc), dead.vegc.sd = sd(dead.vegc),
                rootc.m = mean(rootc), rootc.sd = sd(rootc),
                stemc.m = mean(stemc), stemc.sd = sd(stemc),
                leafc.m = mean(cs.leafc), leafc.sd = sd(cs.leafc),
                lai.m = mean(epv.proj_lai), lai.sd = sd(epv.proj_lai),
                swo.m = mean(wf.soilw_outflow), swo.sd = sd(wf.soilw_outflow),
                psi.m = mean(epv.psi), psi.sd = sd(epv.psi),
                psi_leaf.m = mean(epv.psi_leaf), psi_leaf.sd = sd(epv.psi_leaf),
                wueit.m = mean(WUEi.t), wueit.sd = sd(WUEi.t),
                wueet.m = mean(WUEe.t), wueet.sd = sd(WUEe.t))
  
  Annualized
  # unique(Daily$simset)
  Annualized$species2<-factor(as.character(Annualized$species2), levels=c("lepidodendron","macroneuropteris","treefern","cordaites","ebf"))
  
  summary(Annualized$m_psi_x.m)
  # Annualized.100.summary <- Daily %>% filter(num.mprcp==100) %>% group_by(gcm,species,simyearday) %>%
        # summarize(sw.m = mean(ws.soilw), sw.s = sd(ws.soilw), psi.m = mean(epv.psi), psi.s = sd(epv.psi),)
 
  #swo : kg/m2/d 
  #       *365 d/y -> kg/m2/yr
  #       *10^6 km2/m2 -> kg/km2/yr
  #       *4x10^6 km2/congo_km2 -> kg/congo_km2/yr
  #       /9.97×10^11 kg/km3 water -> km3/congo/yr

  
  
  D.stats <- Daily %>% 
    filter(num.mprcp==100, 
           species2 %in% c('cordaites', 'lepidodendron','macroneuropteris','treefern','ebf'),
           # species %in% c('cordaites', 'lepidodendron','macroneuropteris','treefern',
           #                'ebf','enf','ebfArp','ebfMA'),
           simset %in% c('base','stem_leafsoil')) %>%  
    group_by(simset,gcm,species) %>% 
    summarize(species2 = unique(species2),
              min_leaf_psi = min(epv.psi_leaf),
              leaf_psi = mean(epv.psi_leaf),
              psi = mean(epv.psi),
              congo_swo = mean(wf.soilw_outflow*(365)*(10^6)*(4*10^6)/(9.97*10^11))   )
  
  L.p50 <- -2.2
  M.p50 <- -1.5
  TF.p50 <- -3.89
  C.p50 <- -6.2
  
  D.stats$p50<-as.numeric(NA)
  D.stats[which(D.stats$species == "lepidodendron"),]$p50<-L.p50
  D.stats[which(D.stats$species == "macroneuropteris"),]$p50<-M.p50
  D.stats[which(D.stats$species == "treefern"),]$p50<-TF.p50
  D.stats[which(D.stats$species == "cordaites"),]$p50<-C.p50
  
  D.stats$s_m <- D.stats$min_leaf_psi-D.stats$p50
  
  print(D.stats %>% dplyr::select(species,min_leaf_psi,p50,s_m),n=22)
  
    
  }
  
 #subsets and selective removals
  if(TRUE){
  unique(Daily$simset)
  Daily.noRev<-Daily %>% filter(simset %in% c("base","stem_leafsoil"))
  Daily.100<-Daily %>% filter(mprcp %in% c("100%")) 
  Daily.100.noRev<-Daily %>% filter(simset %in% c("base","stem_leafsoil"), mprcp %in% c("100%"))
  
  # Daily.100.G<-Daily.100 %>% filter(co2 %in% c("182")) 
  #
  
  ###selective removals
  #rm(list=setdiff(ls(),c("Daily","base_head","stem_head","Daily.100.noRev.l")))
  #rm(list=setdiff(ls(),c("Daily","base_head","stem_head")))
  #rm(Daily.100)
  # unique(Daily.noRev.cltf$species)
  # 
  # Daily.noRev.cltf<-Daily %>% 
  #   filter(simset %in% c("base","stem")) %>%
  #   filter(species %in% c("cordaites","lepidodendron","treefern"))
  
  #single sims
  if(TRUE){
    vCO2<-"546"
    
    Daily.base.c <- Daily.100%>%
      filter(co2 %in% c(vCO2)) %>% 
      filter(simset %in% c("base"))%>%
      filter(species %in% c("cordaites"))
    
  Daily.stem.c <- Daily.100%>%
    filter(co2 %in% c(vCO2)) %>% 
    filter(simset %in% c("stem_leafsoil"))%>%
    filter(species %in% c("cordaites"))
  
  # Daily.rev.c <-Daily.100%>%
  #   filter(co2 %in% c(vCO2)) %>% 
  #   filter(simset %in% c("stemRev"))%>%
  #   filter(species %in% c("cordaitesRev")) 

  Daily.base.l <- Daily.100%>%
    filter(co2 %in% c(vCO2)) %>% 
    filter(simset %in% c("base"))%>%
    filter(species %in% c("lepidodendron"))
    
  Daily.stem.l <- Daily.100%>%
    filter(co2 %in% c(vCO2)) %>% 
    filter(simset %in% c("stem_leafsoil"))%>%
    filter(species %in% c("lepidodendron"))
  
  # Daily.rev.l <- Daily.100%>%
  #   filter(co2 %in% c(vCO2)) %>% 
  #   filter(simset %in% c("stemRev"))%>%
  #   filter(species %in% c("lepidodendronRev"))

  Daily.base.m <- Daily.100%>%
    filter(co2 %in% c(vCO2)) %>% 
    filter(simset %in% c("base"))%>%
    filter(species %in% c("macroneuropteris"))
  
   Daily.stem.m <- Daily.100%>%
     filter(co2 %in% c(vCO2)) %>% 
    filter(simset %in% c("stem_leafsoil"))%>%
    filter(species %in% c("macroneuropteris"))
  
  # Daily.rev.m <- Daily.100%>%
  #   filter(co2 %in% c(vCO2)) %>% 
  #   filter(simset %in% c("stemRev"))%>%
  #   filter(species %in% c("macroneuropterisRev"))

  Daily.base.p <- Daily.100%>%
    filter(co2 %in% c(vCO2)) %>% 
    filter(simset %in% c("base"))%>%
    filter(species %in% c("treefern"))
    
  Daily.stem.p <- Daily.100%>%
    filter(co2 %in% c(vCO2)) %>% 
    filter(simset %in% c("stem_leafsoil"))%>%
    filter(species %in% c("treefern"))
  
  # Daily.rev.p <- Daily.100%>%
  #   filter(co2 %in% c(vCO2)) %>% 
  #   filter(simset %in% c("stemRev"))%>%
  #   filter(species %in% c("treefernRev"))
  
  Daily.stem.Arp <- Daily.100%>%
    filter(co2 %in% c(vCO2)) %>% 
    filter(simset %in% c("base"))%>%
    filter(species %in% c("ebfArp"))
    
  Daily.stem.Arp <- Daily.100%>%
    filter(co2 %in% c(vCO2)) %>% 
    filter(simset %in% c("stem_leafsoil"))%>%
    filter(species %in% c("ebfArp"))
    
  # Daily.rev.Arp <- Daily.100%>%
  #   filter(co2 %in% c(vCO2)) %>% 
  #   filter(simset %in% c("stemRev"))%>%
  #   filter(species %in% c("ebfArpRev"))

  Daily.base.MA <- Daily.100%>%
    filter(co2 %in% c(vCO2)) %>% 
    filter(simset %in% c("base"))%>%
    filter(species %in% c("ebfMA"))
    
  Daily.stem.MA <- Daily.100%>%
    filter(co2 %in% c(vCO2)) %>% 
    filter(simset %in% c("stem_leafsoil"))%>%
    filter(species %in% c("ebfMA"))
  
  # Daily.rev.MA <- Daily.100%>%
  #   filter(co2 %in% c(vCO2)) %>% 
  #   filter(simset %in% c("stemRev"))%>%
  #   filter(species %in% c("ebfMARev"))
  # 
  }
  #simsets
  if(TRUE){
    Daily.G.c <- Daily.100%>%
      filter(co2 %in% c("546")) %>% 
      filter(species %in% c("cordaites","cordaitesRev"))
    
    Daily.G.l <- Daily.100%>%
      filter(co2 %in% c("546")) %>% 
      filter(species %in% c("lepidodendron", "lepidodendronRev"))
    
   Daily.G.m <- Daily.100%>%
      filter(co2 %in% c("546")) %>% 
      filter(species %in% c("macroneuropteris", "macroneuropterisRev"))
 
    Daily.G.p <- Daily.100%>%
      filter(co2 %in% c("546")) %>% 
      filter(species %in% c("treefern", "treefernRev"))
    
    Daily.G.Arp <- Daily.100%>%
      filter(co2 %in% c("546")) %>% 
      filter(species %in% c("ebfArp","ebfArpRev"))
    
    Daily.G.MA <- Daily.100%>%
      filter(co2 %in% c("546")) %>% 
      filter(species %in% c("ebf","ebfMA"))
  }
  
  ###100 only
  Daily.100.noRev.cltf<-Daily.noRev.cltf %>% filter(mprcp %in% c("100%")) 
  #model
  Daily.100.base.cltf<-Daily.100.noRev.cltf %>% filter(simset %in% c("base"))
  Daily.100.stem.cltf<-Daily.100.noRev.cltf %>%  filter(simset %in% c("stem"))
  #gcm
  Daily.100.noRev.cltf.I <- Daily.100.noRev.cltf%>% filter(gcm %in% "546_28")
  Daily.100.noRev.cltf.G <- Daily.100.noRev.cltf%>% filter(gcm %in% "182_28")
  #species
  Daily.100.noRev.l <- Daily.100.noRev.cltf %>% filter(species %in% c("lepidodendron"))
  
  ###given 100, seperated gcms, break out by species
  # Daily.100.noRev.c.I<-Daily.100.noRev.cltf.I%>%filter(species %in% "cordaites")
  # Daily.100.noRev.tf.I<-Daily.100.noRev.cltf.I%>%filter(species %in% "treefern")
  # Daily.100.noRev.l.I<-Daily.100.noRev.cltf.I%>%filter(species %in% "lepidodendron")
  # Daily.100.noRev.c.G<-Daily.100.noRev.cltf.G%>%filter(species %in% "cordaites")
  # Daily.100.noRev.tf.G<-Daily.100.noRev.cltf.G%>%filter(species %in% "treefern")
  # Daily.100.noRev.l.G<-Daily.100.noRev.cltf.G%>%filter(species %in% "lepidodendron")
  
  # Daily.100.stem.c.I<-  Daily.100.stem.cltf%>% filter(gcm %in% "546_28")%>%
  #   filter(species %in% c("cordaites"))
  # Daily.100.base.c.I<-  Daily.100.base.cltf%>% filter(gcm %in% "546_28")%>%
  #   filter(species %in% c("cordaites"))
  # Daily.100.stem.c.I<-  Daily.100.stem.cltf%>% filter(gcm %in% "546_28")%>%
  #   filter(species %in% c("treefern"))
  # Daily.100.base.c.I<-  Daily.100.base.cltf%>% filter(gcm %in% "546_28")%>%
  #   filter(species %in% c("treefern"))
  
  
  ### 50,40,30 only (chase down lep death)
  # Daily.30to50.noRev.l.IG<-Daily.noRev.cltf %>% filter(mprcp %in% c("30%","40%","50%")) %>%
  #   filter(species %in% c("lepidodendron")) %>% filter(gcm %in% "546_28")
  Daily.10to50.noRev.l.IG<-Daily.noRev.cltf %>% filter(mprcp %in% c("10%","20%","30%","40%","50%")) %>%
    filter(species %in% c("lepidodendron")) %>% filter(gcm %in% "546_28")



  # Daily.40.noRev.l.IG<-Daily.noRev.cltf %>% filter(mprcp %in% c("40%")) %>%
  #   filter(species %in% c("lepidodendron")) %>% filter(gcm %in% "546_28")
  # Daily.50.noRev.l.IG<-Daily.noRev.cltf %>% filter(mprcp %in% c("50%")) %>%
  #   filter(species %in% c("lepidodendron")) %>% filter(gcm %in% "546_28")
  }

}

# pA<-ggplot(data=Daily.100)+
#   geom_point(
#     size=0.1 , color = 'grey', alpha = 0.5, aes(x=simyearday, y=ws.soilw, shape=as.factor(species)))+
#   geom_smooth(aes(x=simyearday, y=ws.soilw, color = as.factor(simset)))+
#   facet_grid(climate~.)+theme_minimal()


# unique(Daily.100.base.c.I$num.mprcp)
# hist(Daily.100.base.c.I$gls.t, breaks=100, xlab = "Daily G.s", main = "Base Model: IGLAC 100%prcp Cordaites")
# hist(Daily.100.stem.c.I$gls.t, breaks=100, xlab = "Daily G.s", main = "Stem Model: IGLAC 100%prcp Cordaites")
# hist(Daily.100.base.c.I$A.t, breaks=100, xlab = "Daily A.rate", main = "Base Model: IGLAC 100%prcp Cordaites")
# hist(Daily.100.stem.c.I$A.t, breaks=100, xlab = "Daily A.rates", main = "Stem Model: IGLAC 100%prcp Cordaites")
# ggplot(data=Daily.100.stem.c.I)+
#   geom_point(aes(epv.psi_leaf, gls.t))+theme_classic()
# ggplot(data=Daily.100.base.c.I)+
#   geom_point(aes(epv.psi_leaf, gls.t))+theme_classic()
# Terry Pitts - Special Investigator - OPM - Defense counter intelligence

#bespoke
if(TRUE){
  
  
  
  pA<- ggplot(data=Daily.stem.c)+
    geom_point(aes(x=epv.psi_leaf,y=1-m_psi_x))+
    theme_void()
  plotName<-"PLC.example.cord.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  load(file = "~/Desktop/tracheid sem/air_seeding.RData")
  hackedat_VA$psi_p_50 <- -hackedat_VA$psi_p_50
  
  temp <- rbind(temp, hackedat_VA)
  
  cbbPalette <- c("#E69F00", "#009E73", "#0072B2", "#56B4E9","#000000", "#D55E00", "#CC79A7")
  # temp$Group<-as.factor(temp$Group)
  # temp$Group <- relevel(temp$Group, "Eudicot")
  pA<- ggplot()+
    geom_point(data=temp[which(temp$Group!="Eudicot"),],
               size = 4,aes(pit_area, psi_p_50, color=Group))+
    geom_errorbar(data=temp[which(temp$Group!="Eudicot"),],
                  width = .5,aes(x = pit_area, ymin=psi_p_50-sd, ymax=psi_p_50+sd))+
    geom_point(data=temp[which(temp$Group %in% c("Eudicot", "Vesselless Angiosperm")),],
              color="black", size = 2, aes(pit_area, psi_p_50))+
    geom_function(fun = f, color = "red")+guides(color=guide_legend(title="Fossil Group"))+
    scale_color_manual(values=cbbPalette)+theme_classic()+theme(text= element_text(size=15))
  plotName<-"eudicot_VA_empirical_regressed_fossil_psi.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  med.simDist<-data.frame(cbind(med.Rdist, med.Ldist, med.Apdist.mm, med.Psidist))
  
  pA<-ggplot(data = med.simDist)+
    geom_histogram(aes(x=med.Rdist,y=stat(density)),bins = 100)+
    theme_bw()+xlab("Radius (um)")+ylab("Count")+
    theme(text= element_text(size=25), axis.title.y = element_blank(),axis.title.x = element_blank())
  plotName<-"med.Rdist.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 5.5, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data = med.simDist)+
    geom_histogram(aes(x=med.Ldist/1000, y=stat(density)),bins = 100)+
    theme_bw()+xlab("Length (mm)")+ylab("Count")+
    theme(text= element_text(size=25), axis.title.y = element_blank(),axis.title.x = element_blank())
  plotName<-"med.Ldist.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 5.5, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data = med.simDist)+
    geom_histogram(aes(x=med.Apdist.mm, y=stat(density)),bins = 100)+
    theme_bw()+xlab("Pit Area per Tracheid (mm2)")+ylab("Count")+
    theme(text= element_text(size=25), axis.title.y = element_blank(),axis.title.x = element_blank())
  plotName<-"med.Apdist.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 5.5, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data = med.simDist)+
    geom_histogram(aes(x=med.Psidist, y=stat(density)),bins = 100)+
    theme_bw()+xlab("Pit Area per Tracheid (mm2)")+ylab("Count")+
    theme(text= element_text(size=25), axis.title.y = element_blank(),axis.title.x = element_blank())
  plotName<-"med.Psidist.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 5.5, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data = med.simDist)+
    geom_histogram(aes(x=med.Psidist, y = stat(density)),bins = 100)+
    stat_function(fun = dnorm,
                  args = list(mean = mean(med.Psidist),
                              sd = sd(med.Psidist)),
                  col = "red",
                  size = 2.5)+
    theme_bw()+xlab("Pit Area per Tracheid (mm2)")+ylab("Count")+
    theme(text= element_text(size=25), axis.title.y = element_blank(),axis.title.x = element_blank())
  plotName<-"med.Psidist.overlay.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 5.5, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  
  pA<- ggplot(data=Annualized)+
    geom_boxplot(aes(x=species,y=glst.m,fill=simset))+
    facet_grid(.~gcm)+
    theme_minimal()
  plotName<-"glst.box.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
 
  
  pA<- ggplot(data=Annualized)+
    geom_point(position = "dodge",aes(x=species2,y=vegc.m, color = as.factor(simset)))+
    geom_errorbar(position = "dodge",width = 0.5, aes(x=species2, ymin=vegc.m-vegc.sd, ymax=vegc.m+vegc.sd, color = as.factor(simset)))+
    facet_grid(gcm~.)+
    theme_minimal()
  plotName<-"dot.vegc.m.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  # & Annualized$simset %in% c("stem_leafsoil")
  # Annualized$species2 %in% c("cordaites","macroneuropteris"))
  cbbPalette <- c("#E69F00", "#009E73", "#0072B2", "#56B4E9","#000000", "#D55E00", "#CC79A7")
  
  
  pA<- ggplot(data=Annualized[which(Annualized$species2 %in% c("cordaites")),])+
    geom_point(size = 2, aes(x=simyearday,y=psi_leaf.m, color = as.factor(simset), shape = as.factor(simset)))+
    geom_errorbar(alpha = 0.33, aes(x=simyearday,ymin=psi_leaf.m-psi_leaf.sd, ymax=psi_leaf.m+psi_leaf.sd, color = as.factor(simset)))+
    geom_hline(yintercept=-6.2, linetype=2, color="black", size=1)+
    scale_color_manual(values = cbbPalette[c(5,6)])+facet_grid(.~gcm)+
    theme_bw()+theme(text= element_text(size=25), strip.background = element_blank())
  plotName<-"psi_leaf_by_yearday.C2.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 12, units = c("in"),
         dpi = 300, limitsize = TRUE)

  pA<- ggplot(data=Annualized[which(Annualized$species2 %in% c("macroneuropteris")),])+
    geom_point(size = 2, aes(x=simyearday,y=psi_leaf.m, color = as.factor(simset), shape = as.factor(simset)))+
    geom_errorbar(alpha = 0.33, aes(x=simyearday,ymin=psi_leaf.m-psi_leaf.sd, ymax=psi_leaf.m+psi_leaf.sd, color = as.factor(simset)))+
    geom_hline(yintercept=-1.5, linetype=2, color="black", size=1)+
    scale_color_manual(values = cbbPalette[c(5,6)])+facet_grid(.~gcm)+
    theme_bw()+theme(text= element_text(size=25), strip.background = element_blank())
  plotName<-"psi_leaf_by_yearday.M2.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 12, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  Annualized.C.wide <- Annualized %>% dplyr::select(gcm, simset, species2,simyearday,rootc.m,stemc.m,leafc.m) %>%
    pivot_longer(cols = c(rootc.m, stemc.m, leafc.m),names_to = "c_type") 
  
  pA<- ggplot(data=Annualized.C.wide)+
    geom_bar(position="stack", stat="identity", aes(x=species2,y=value, fill = c_type))+
    facet_grid(gcm~simset)+#scale_color_manual(values = cbbPalette[c(5,6)])
    theme_bw()+theme(text= element_text(size=25), strip.background = element_blank())
  plotName<-"c.stacks.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 12, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
}

#boxplots
in(TRUE){
  cbbPalette <- c("#E69F00", "#009E73", "#0072B2", "#56B4E9","#000000", "#D55E00", "#CC79A7")
  #
  # swo.m*1000000/(1000*86400) ~ cubic meters per sq km per s
  #     *35.31 (cubic ft per cubic meter)~ to cubic feet per km per s
  #     *0.3861 (square miles per square kilometer)    
  #
  pA <- ggplot()+geom_boxplot(data=Annualized, aes(x=species2, y=swo.m, color = simset))+
    scale_color_manual(values=cbbPalette[c(5,6)])+
    facet_grid(gcm~.,scales = "free_x")+
    theme_bw()+theme(text= element_text(size=25), strip.background = element_blank()) # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  plotName<-"Annualized.swo.box.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 9, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA <- ggplot()+geom_boxplot(data=Annualized, aes(x=species2, y=glst.m, color = simset))+
    scale_color_manual(values=cbbPalette[c(5,6)])+
    facet_grid(gcm~.,scales = "free_x")+
    theme_bw()+theme(text= element_text(size=25), strip.background = element_blank())+ 
    facet_grid(gcm~.,scales = "free_x")
  plotName<-"Annualized.glst.box.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 9, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA <- ggplot()+geom_boxplot(data=Annualized, aes(x=species2, y=At.m, color = simset))+
    scale_color_manual(values=cbbPalette[c(5,6)])+
    facet_grid(gcm~.,scales = "free_x")+
    theme_bw()+theme(text= element_text(size=25), strip.background = element_blank())+ 
    facet_grid(gcm~.,scales = "free_x")
  plotName<-"Annualized.At.box.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  # unique(is.infinite(Annualized$At.m))
  # unique(Annualized$glst.m==0)
  # plot(Annualized$At.m/Annualized$glst.m)
  
  pA <- ggplot()+geom_boxplot(data=Annualized, aes(x=species2, y=wueit.m, color = simset))+
    scale_color_manual(values=cbbPalette[c(5,6)])+
    facet_grid(gcm~.,scales = "free_x")+
    theme_bw()+theme(text= element_text(size=25), strip.background = element_blank())+ 
    facet_grid(gcm~.,scales = "free_x")
  plotName<-"Annualized.wueit.box.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA <- ggplot()+geom_boxplot(data=Annualized, aes(x=species2, y=wueet.m, color = simset))+
    scale_color_manual(values=c("grey75","grey25"))+
    theme_minimal()+
    facet_grid(gcm~fossil,scales = "free_x")
  plotName<-"Annualized.wueet.box.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA <- ggplot()+geom_boxplot(data=Annualized, aes(x=species2, y=live.vegc.m, color = simset))+
    scale_color_manual(values=cbbPalette[c(5,6)])+
    theme_minimal()+
    facet_grid(gcm~.,scales = "free_x")
  plotName<-"Annualized.lvegc.box.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA <- ggplot()+geom_boxplot(data=Annualized, aes(x=species2, y=vegc.m, color = simset))+
    scale_color_manual(values=cbbPalette[c(5,6)])+
    theme_minimal()+
    facet_grid(gcm~.,scales = "free_x")
  plotName<-"Annualized.vegc.box.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA <- ggplot()+geom_boxplot(data=Annualized, aes(x=species2, y=lai.m, color = simset))+
    scale_color_manual(values=cbbPalette[c(5,6)])+
    theme_minimal()+
    facet_grid(gcm~.,scales = "free_x")
  plotName<-"Annualized.laim.box.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  unique(Annualized$gcm)
  pA <- ggplot()+geom_boxplot(data=Annualized, aes(x=species2, y=lswr.m, color = simset))+
    scale_color_manual(values=cbbPalette[c(5,6)])+
    theme_bw()+theme(text = element_text(size = 20))
    # facet_grid(gcm~.,scales = "free_x")
  plotName<-"Annualized.lswr.box.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA <- ggplot()+geom_boxplot(data=Annualized, aes(x=species2, y=m_psi_x.m, color = simset))+
    scale_color_manual(values=cbbPalette[c(5,6)])+
    theme_bw()+theme(text = element_text(size = 20))+
    facet_grid(gcm~.,scales = "free_x")
  plotName<-"Annualized.mx.box.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)

  pA <- ggplot()+geom_boxplot(data=Annualized[which(Annualized$simset=="stem_leafsoil"),], aes(x=species, y=lswr.m, color = simset))+
    scale_color_manual(values=cbbPalette[c(6)])+
    theme_bw()+theme(text = element_text(size = 20))
    # facet_grid(gcm~.,scales = "free_x")
  plotName<-"Annualized.lswr.s.box.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 3.5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
    
  pA <- ggplot()+geom_boxplot(data=Annualized[which(Annualized$simset=="stem_leafsoil"),], aes(x=species, y=m_psi_x.m, color = simset))+
    scale_color_manual(values=cbbPalette[c(6)])+
    theme_bw()+theme(text = element_text(size = 20))
    # facet_grid(gcm~.,scales = "free_x")
  plotName<-"Annualized.mx.s.box.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 3.5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
    # kgC / m2 / d
  #       *10000 m2 / ha
  #       *365 days/yr
  # kgC / ha / yr
  pA <- ggplot()+geom_boxplot(data=Annualized, aes(x=species2, y=nmin.m*3.65*10^6, color = simset))+
    scale_color_manual(values=cbbPalette[c(5,6)])+
    facet_grid(gcm~.,scales = "free_x")+
    theme_bw()+theme(text= element_text(size=25), strip.background = element_blank())   
  plotName<-"Annualized.nmin.box.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 9, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  # kgC / m2 / d
  #       *10000 m2 / ha
  #       *365 days/yr
  #       *1000 g/kg
  # gC / ha / yr
  # summary(Annualized$nlea.m*3.65*10^9)
  pA <- ggplot()+geom_boxplot(data=Annualized, aes(x=species2, y=nlea.m*3.65*10^6, color = simset))+
    scale_color_manual(values=cbbPalette[c(5,6)])+
    facet_grid(gcm~.,scales = "free_x")+
    theme_bw()+theme(text= element_text(size=25), strip.background = element_blank())
  plotName<-"Annualized.nlea.box.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 9, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  
}

#daily psi histogram
if(TRUE){
  
  max(Daily.stem.c$epv.psi)
  min(Daily.stem.c$epv.psi)
  max(Daily.stem.c$epv.psi_leaf)
  # print(Daily.stem.c[which(Daily.stem.c$epv.psi>-0.008),], max = 2000)
  summary(Daily.stem.c$m_psi_x)
  median(Daily.stem.c$m_psi_x)
  
  max(Daily.stem.l$epv.psi)
  min(Daily.stem.l$epv.psi)
  max(Daily.stem.l$epv.psi_leaf)
  # print(Daily.stem.l[which(Daily.stem.l$epv.psi>-0.008),], max = 2000)
  summary(Daily.stem.l$m_psi_x)
  median(Daily.stem.l$m_psi_x)
  
  max(Daily.stem.m$epv.psi)
  min(Daily.stem.m$epv.psi)
  max(Daily.stem.m$epv.psi_leaf)
  # print(Daily.stem.l[which(Daily.stem.l$epv.psi>-0.008),], max = 2000)
  summary(Daily.stem.m$m_psi_x)
  median(Daily.stem.m$m_psi_x)
  
  
  max(Daily.stem.p$epv.psi)
  min(Daily.stem.p$epv.psi)
  max(Daily.stem.p$epv.psi_leaf)
  # print(Daily.stem.p[which(Daily.stem.m$epv.psi>-0.002),], max = 2000)
  summary(Daily.stem.p$m_psi_x)
  median(Daily.stem.p$m_psi_x)
  
  max(Daily.stem.Arp$epv.psi)
  min(Daily.stem.Arp$epv.psi)
  max(Daily.stem.Arp$epv.psi_leaf)
  # print(Daily.stem.p[which(Daily.stem.m$epv.psi>-0.002),], max = 2000)
  summary(Daily.stem.Arp$m_psi_x)
  median(Daily.stem.Arp$m_psi_x)
  
  cbbPalette <- c("#E69F00", "#009E73", "#0072B2", "#56B4E9","#000000", "#D55E00", "#CC79A7")
  
  print(Daily.100 %>% group_by(simset, species) %>% 
    summarize(psi.m = mean(epv.psi),psi.sd = sd(epv.psi),
              psi_leaf.m = mean(epv.psi_leaf), psi_leaf.sd = sd(epv.psi_leaf),
              m_x.m = mean(m_psi_x), m_x.sd = sd(m_psi_x),
              count = n()), n= 25) 
  
  scale_x_continuous(limits = c(-17, 0))+
  pA <- ggplot(data=Daily.100)+scale_x_continuous(limits = c(-17, -0))+
    geom_histogram(aes(x=epv.psi*1000, fill=as.factor(simset)),bins = 100, alpha=0.3, position="identity" )+
    geom_vline(xintercept=-6, linetype=1, color="black", size=1)+
    geom_vline(xintercept=-17, linetype=2, color="black", size=1)+
    scale_fill_manual(values = cbbPalette[c(5,6)])+theme_bw()+
    facet_grid(.~gcm)+theme(text= element_text(size=25), strip.background = element_blank())
  plotName<-"daily.psi.close.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 12, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA <- ggplot(data=Daily.100)+
          geom_histogram(aes(x=epv.psi_leaf, color=as.factor(simset)),bins = 100)+
          theme_bw()
  plotName<-"daily.psi_leaf.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
                      scale = 1, height = 5, width = 7, units = c("in"),
                      dpi = 300, limitsize = TRUE)
               
  
}

#psi leaf by m_psi_x
if(TRUE){
  # erfc<-rlang::as_function(~1-erf(.x))
  # 
  # v<-rlang::as_function(~0.5 * erfc((.mu - .psi)/(sqrt(2.0)*.s)))
  # v(-3,-1.5,0.5)
  # erf()
  # f<-rlang::as_function(~-1.538*(.x)^(-0.402))  
  Annualized.stem.mc <- Annualized %>% filter(species %in% c("cordaites","macroneuropteris"), simset %in% c("stem_leafsoil") )
  
  pA<- ggplot(data=Annualized.stem.mc)+
    geom_point(aes(x=psi_leaf.m,y=1-m_psi_x.m, color = species))+
    theme_bw()+theme(text= element_text(size=25), strip.background = element_blank())
  plotName<-"Annualized.stem.mc.mpsix.by.psisoil.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  Daily.stem.m <- Daily %>% filter(species %in% c("macroneuropteris"), simset %in% c("stem_leafsoil") )
  
  pA<- ggplot(data=Daily.stem.m)+
    geom_point(aes(x=epv.psi,y=1-m_psi_x))+
    theme_bw()
  plotName<-"daily.mpsix.by.psisoil.cord.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<- ggplot(data=Daily.stem.c)+
    geom_point(aes(x=epv.psi,y=1-m_psi_x))+
    theme_bw()
  plotName<-"daily.mpsix.by.psisoil.cord.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)


  
  # scale_x_continuous(limits = c(-0.5, 0))+
  #   scale_y_continuous(limits = c(-0.5, 0))+

}

#psi leaf by psi soil
if(TRUE){
  max(Daily.stem.l$epv.psi)
  pA<- ggplot(data=Daily.stem.l)+
    geom_point(aes(x=epv.psi_leaf,y=epv.psi))+
    facet_grid(num.mprcp~gcm)+theme_bw()
  plotName<-"daily.psi.by.psisoil.lep.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
 
  # scale_x_continuous(limits = c(-0.5, 0))+
  #   scale_y_continuous(limits = c(-0.5, 0))+

  pA<- ggplot(data=Daily.stem.c)+
    geom_point(aes(x=epv.psi_leaf,y=epv.psi))+
    facet_grid(num.mprcp~gcm)+theme_bw()
  plotName<-"daily.psi.by.psisoil.cord.zoom.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<- ggplot(data=Daily.rev.l)+
    geom_point(aes(x=epv.psi_leaf,y=epv.psi))+
    facet_grid(num.mprcp~gcm)+theme_bw()
  plotName<-"daily.psi.by.psisoil.lepRev.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<- ggplot(data=Daily.rev.c)+
    geom_point(aes(x=epv.psi_leaf,y=epv.psi))+
    facet_grid(num.mprcp~gcm)+theme_bw()
  plotName<-"daily.psi.by.psisoil.cordRev.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
}

#mpsix by yearday
if(TRUE){
  cbbPalette <- c("#E69F00", "#009E73", "#0072B2", "#56B4E9","#000000", "#D55E00", "#CC79A7")
  
  pA <- ggplot(data=Annualized[which(Annualized$simset=="base"),],aes(x=simyearday, y=m_psi_x.m, color=gcm))+
    geom_point()+
    geom_errorbar(alpha=0.5, aes(x=simyearday, ymin=m_psi_x.m-m_psi_x.sd, ymax=m_psi_x.m+m_psi_x.sd ))+
    scale_color_manual(values = cbbPalette)+facet_grid(species~., scales = 'free_y')+theme_minimal()
  plotName<-"Annualized.base.m_psi_x_by_yearday.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)

  pA <- ggplot(data=Annualized[which(Annualized$simset=="stem_leafsoil" & Annualized$species %in% c("macroneuropteris","cordaites") & Annualized$gcm == "546_28"),])+
    geom_point(aes(x=simyearday, y=m_psi_x.m, color=species))+
    geom_line(data=Annualized[which(Annualized$simset=="stem_leafsoil" & Annualized$species =="macroneuropteris" & Annualized$gcm == "546_28"),],
                               size = 1.0, linetype=2,aes(x=simyearday, y=lswr.m))+
    geom_errorbar(alpha=0.5,aes(x=simyearday, ymin=m_psi_x.m-m_psi_x.sd, ymax=m_psi_x.m+m_psi_x.sd, color=species))+
    scale_color_manual(values = cbbPalette[c(6,2)])+
    theme_bw()+theme(text= element_text(size=25), strip.background = element_blank())
  plotName<-"Annualized.stem.m_psi_x_by_yearday.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 9, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  # 
  # Daily.stem.glac.l<-Daily.stem.l %>% filter(gcm %in% c("182_28"))
  # pA<- ggplot(data=Daily.stem.glac.l)+
  #   geom_point(aes(x=simyearday,y=m_psi_x))+
  #   facet_grid(num.mprcp~.,scales = "free_y")+theme_bw()
  # plotName<-"daily.mpsix.by.yearday.stem.glac.lep.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 7, width = 7, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  # 
  # Daily.stem.iglac.l<-Daily.stem.l %>% filter(gcm %in% c("546_28"))
  # pA<- ggplot(data=Daily.stem.iglac.l)+
  #   geom_point(aes(x=simyearday,y=m_psi_x))+
  #   facet_grid(num.mprcp~.,scales = "free_y")+theme_bw()
  # plotName<-"daily.mpsix.by.yearday.stem.iglac.lep.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 7, width = 7, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  # 
  # Daily.stem.glac.c<-Daily.stem.c %>% filter(gcm %in% c("182_28"))
  # pA<- ggplot(data=Daily.stem.glac.c)+
  #   geom_point(aes(x=simyearday,y=m_psi_x))+
  #   facet_grid(num.mprcp~.,scales = "free_y")+theme_bw()
  # plotName<-"daily.mpsix.by.yearday.stem.glac.cord.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 7, width = 7, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  # 
  # Daily.stem.iglac.c<-Daily.stem.c %>% filter(gcm %in% c("546_28"))
  # pA<- ggplot(data=Daily.stem.iglac.c)+
  #   geom_point(aes(x=simyearday,y=m_psi_x))+
  #   facet_grid(num.mprcp~.,scales = "free_y")+theme_bw()
  # plotName<-"daily.mpsix.by.yearday.stem.iglac.cord.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 7, width = 7, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  # 
  # 
  
}

#10 to 50 A etc by yearday
if(TRUE){
  ##################################################################  
  # [which(Daily.10to50.noRev.l.IG$simset=="stem"),]
  pA<-ggplot(data=Daily.10to50.noRev.l.IG[which(Daily.10to50.noRev.l.IG$simset=="stem"),])+
    geom_point(
      size=0.5 , alpha = 0.5, aes(x=simyearday, y=epv.proj_lai))+
    facet_grid(mprcp~simset, scales = "free_y")+theme_bw()+theme(text= element_text(size=15))
  plotName<-"Stem.10to50_Lep_IG_lai_free.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=Daily.10to50.noRev.l.IG[which(Daily.10to50.noRev.l.IG$simset=="stem"),])+
    geom_point(
      size=0.5 , alpha = 0.5, aes(x=simyearday, y=A.t))+
    facet_grid(mprcp~simset, scales = "free_y")+theme_bw()+theme(text= element_text(size=15))
  plotName<-"Stem.10to50_Lep_IG_A.t.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=Daily.10to50.noRev.l.IG[which(Daily.10to50.noRev.l.IG$simset=="stem"),])+
    geom_point(
      size=0.5 , alpha = 0.5, aes(x=simyearday, y=psn_sun.g))+
    facet_grid(mprcp~simset, scales = "free_y")+theme_bw()+theme(text= element_text(size=15))
  plotName<-"Stem.10to50_Lep_IG_psn_sun.g.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=Daily.10to50.noRev.l.IG[which(Daily.10to50.noRev.l.IG$simset=="stem"),])+
    geom_point(
      size=0.5 , alpha = 0.5, aes(x=simyearday, y=gls.t))+
    facet_grid(mprcp~simset, scales = "free_y")+theme_bw()+theme(text= element_text(size=15))
  plotName<-"Stem.10to50_Lep_IG_psn_sun.glst.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  ##################
  ##################
  pA<-ggplot(data=Daily.10to50.noRev.l.IG[which(Daily.10to50.noRev.l.IG$simset=="base"),])+
    geom_point(
      size=0.5 , alpha = 0.5, aes(x=simyearday, y=epv.proj_lai))+
    facet_grid(mprcp~simset, scales = "free_y")+theme_bw()+theme(text= element_text(size=15))
  plotName<-"Base.10to50_Lep_IG_lai_free.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=Daily.10to50.noRev.l.IG[which(Daily.10to50.noRev.l.IG$simset=="base"),])+
    geom_point(
      size=0.5 , alpha = 0.5, aes(x=simyearday, y=A.t))+
    facet_grid(mprcp~simset, scales = "free_y")+theme_bw()+theme(text= element_text(size=15))
  plotName<-"Base.10to50_Lep_IG_A.t.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=Daily.10to50.noRev.l.IG[which(Daily.10to50.noRev.l.IG$simset=="base"),])+
    geom_point(
      size=0.5 , alpha = 0.5, aes(x=simyearday, y=psn_sun.g))+
    facet_grid(mprcp~simset, scales = "free_y")+theme_bw()+theme(text= element_text(size=15))
  plotName<-"Base.10to50_Lep_IG_psn_sun.g.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=Daily.10to50.noRev.l.IG[which(Daily.10to50.noRev.l.IG$simset=="base"),])+
    geom_point(
      size=0.5 , alpha = 0.5, aes(x=simyearday, y=gls.t))+
    facet_grid(mprcp~simset, scales = "free_y")+theme_bw()+theme(text= element_text(size=15))
  plotName<-"Base.10to50_Lep_IG_psn_sun.glst.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  ############################################################
}

#A.t by gls.t
if(TRUE){
 
   # scale_color_gradient(low = 'blue', high = 'red')+
  assim_set <- c(1,3)
  col_set <- c("black","red")
  unique(Daily.100$species)
  Daily.100 <- Daily %>% filter(num.mprcp == 100)
  Daily.100.l <- Daily.100 %>% filter(species %in% c('lepidodendron'))
  Daily.100.c <- Daily.100 %>% filter(species %in% c('cordaites'))
  Daily.100.m <- Daily.100 %>% filter(species %in% c('macroneuropteris'))
  Daily.100.MA <- Daily.100 %>% filter(species %in% c('ebf', 'ebfMA'))
  Daily.100.Arp <- Daily.100 %>% filter(species %in% c('ebf', 'ebfArp'))
  
  # pA<-ggplot(data=Daily.100.l)+
  #   geom_point(alpha = 0.5, size = 1.5, aes(gls.t, gl.t, color = as.factor(A_type),  shape =as.factor(A_type)))+
  #   theme_bw()+scale_shape_manual(values=assim_set)+scale_color_manual(values=col_set)+
  #   facet_grid(gcm~simset)+theme(text= element_text(size=15))+
  #   labs(shape = "Assimilation Type", color = "Assimilation Type", x = "Total Stomatal Conductance", y = "Total Assimilation")
  # plotName<-"Daily.100.l.gl.t_by_gls.t_Atype_s.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 7, width = 8.5, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  # 
  # pA<-ggplot(data=Daily.100.l)+
  #   geom_point(alpha = 0.5, size = 1.5, aes(gl.t, gco2, color = as.factor(A_type),  shape =as.factor(A_type)))+
  #   theme_bw()+scale_shape_manual(values=assim_set)+scale_color_manual(values=col_set)+
  #   facet_grid(gcm~simset)+theme(text= element_text(size=15))+
  #   labs(shape = "Assimilation Type", color = "Assimilation Type", x = "Total Stomatal Conductance", y = "Total Assimilation")
  # plotName<-"Daily.100.l.gco2_by_gl.t_Atype_s.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 7, width = 8.5, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  # 
  # pA<-ggplot(data=Daily.100.l)+
  #   geom_point(alpha = 0.5, size = 1.5, aes(gco2, A.t, color = as.factor(A_type),  shape =as.factor(A_type)))+
  #   theme_bw()+scale_shape_manual(values=assim_set)+scale_color_manual(values=col_set)+
  #   facet_grid(gcm~simset)+theme(text= element_text(size=15))+
  #   labs(shape = "Assimilation Type", color = "Assimilation Type", x = "Total Stomatal Conductance", y = "Total Assimilation")
  # plotName<-"Daily.100.l.A.t_by_gco2_Atype_s.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 7, width = 8.5, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  # 
  # pA<-ggplot(data=Daily.100.l)+
  #   geom_point(alpha = 0.5, size = 1.5, aes(gl.t, A.t, color = as.factor(A_type),  shape =as.factor(A_type)))+
  #   theme_bw()+scale_shape_manual(values=assim_set)+scale_color_manual(values=col_set)+
  #   facet_grid(gcm~simset)+theme(text= element_text(size=15))+
  #   labs(shape = "Assimilation Type", color = "Assimilation Type", x = "Total Stomatal Conductance", y = "Total Assimilation")
  # plotName<-"Daily.100.l.A.tby_gl.t_Atype_s.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 7, width = 8.5, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=Daily.100.l)+
    geom_point(alpha = 0.5, size = 1.5, aes(gls.t,A.t, color = as.factor(A_type),  shape =as.factor(A_type)))+
    theme_bw()+scale_shape_manual(values=assim_set)+scale_color_manual(values=col_set)+
    facet_grid(gcm~simset)+theme(text= element_text(size=15))+
    labs(shape = "Assimilation Type", color = "Assimilation Type", x = "Total Stomatal Conductance", y = "Total Assimilation")
  plotName<-"Daily.100.l.A.t_by_gls.t_Atype_s.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 7, width = 8.5, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=Daily.100.c)+
    geom_point(alpha = 0.5, size = 1.5, aes(gls.t,A.t, color = as.factor(A_type),  shape =as.factor(A_type)))+
    theme_bw()+scale_shape_manual(values=assim_set)+scale_color_manual(values=col_set)+
    facet_grid(gcm~simset)+theme(text= element_text(size=15))+
    labs(shape = "Assimilation Type", color = "Assimilation Type", x = "Total Stomatal Conductance", y = "Total Assimilation")
  plotName<-"Daily.100.c.A.t_by_gls.t_Atype_s.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 7, width = 8.5, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=Daily.100.m)+
    geom_point(alpha = 0.5, size = 1.5, aes(gls.t,A.t, color = as.factor(A_type),  shape =as.factor(A_type)))+
    theme_bw()+scale_shape_manual(values=assim_set)+scale_color_manual(values=col_set)+
    facet_grid(gcm~simset)+theme(text= element_text(size=15))+
    labs(shape = "Assimilation Type", color = "Assimilation Type", x = "Total Stomatal Conductance", y = "Total Assimilation")
  plotName<-"Daily.100.m.A.t_by_gls.t_Atype_s.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 7, width = 8.5, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=Daily.100.Arp)+
    geom_point(alpha = 0.5, size = 1.5, aes(gls.t,A.t, color = as.factor(A_type),  shape =as.factor(A_type)))+
    theme_bw()+scale_shape_manual(values=assim_set)+scale_color_manual(values=col_set)+
    facet_grid(gcm~simset)+theme(text= element_text(size=15))+
    labs(shape = "Assimilation Type", color = "Assimilation Type", x = "Total Stomatal Conductance", y = "Total Assimilation")
  plotName<-"Daily.100.Arp.A.t_by_gls.t_Atype_s.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 7, width = 8.5, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=Daily.100.MA)+
    geom_point(alpha = 0.5, size = 1.5, aes(gls.t,A.t, color = as.factor(A_type),  shape =as.factor(A_type)))+
    theme_bw()+scale_shape_manual(values=assim_set)+scale_color_manual(values=col_set)+
    facet_grid(gcm~simset)+theme(text= element_text(size=15))+
    labs(shape = "Assimilation Type", color = "Assimilation Type", x = "Total Stomatal Conductance", y = "Total Assimilation")
  plotName<-"Daily.100.MA.A.t_by_gls.t_Atype_s.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 7, width = 8.5, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
}

#drydown assim vs glst, vegc vs glst
if(TRUE){
  Daily.l <- Daily %>% filter(species %in% c('lepidodendron'))
  Daily.l.I <- Daily.l %>% filter(gcm %in% c('546_28'))
  Daily.l.G <- Daily.l %>% filter(gcm %in% c('182_28'))
  # Daily.90 <- Daily %>% filter(num.mprcp == 90)
  # Daily.80 <- Daily %>% filter(num.mprcp == 80)
  # Daily.70 <- Daily %>% filter(num.mprcp == 70)
  # Daily.60 <- Daily %>% filter(num.mprcp == 60)
  # Daily.50 <- Daily %>% filter(num.mprcp == 50)
  # Daily.40 <- Daily %>% filter(num.mprcp == 40)
  # Daily.30 <- Daily %>% filter(num.mprcp == 30)
  # Daily.20 <- Daily %>% filter(num.mprcp == 20)
  # Daily.10 <- Daily %>% filter(num.mprcp == 10)
  
  
  pA<-ggplot(data=Daily.l.I)+
    geom_point(alpha = 0.5, size = 1.5, aes(gls.t,A.t, color = as.factor(A_type),  shape =as.factor(A_type)))+
    theme_bw()+scale_shape_manual(values=assim_set)+scale_color_manual(values=col_set)+
    facet_grid(num.mprcp~simset)+theme(text= element_text(size=15))+
    labs(shape = "Assimilation Type", color = "Assimilation Type", x = "Daily Stomatal Conductance", y = "Daily Assimilation")
  plotName<-"Daily.drydown.l.I.A.t_by_gls.t_Atype_s.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 7, width = 8.5, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=Daily.l.G)+
    geom_point(alpha = 0.5, size = 1.5, aes(gls.t,A.t, color = as.factor(A_type),  shape =as.factor(A_type)))+
    theme_bw()+scale_shape_manual(values=assim_set)+scale_color_manual(values=col_set)+
    facet_grid(num.mprcp~simset)+theme(text= element_text(size=15))+
    labs(shape = "Assimilation Type", color = "Assimilation Type", x = "Daily Stomatal Conductance", y = "Daily Assimilation")
  plotName<-"Daily.drydown.l.G.A.t_by_gls.t_Atype_s.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 7, width = 8.5, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  
  # pA<-ggplot(data=Daily.l.I)+
  #   geom_point(alpha = 0.5, size = 1.5, aes(gls.t,summary.vegc))+
  #   theme_bw()+
  #   facet_grid(num.mprcp~simset)+theme(text= element_text(size=15))+
  #   labs(shape = "Assimilation Type", color = "Assimilation Type", x = "Daily Stomatal Conductance", y = "Daily Assimilation")
  # plotName<-"Daily.drydown.l.I.vegc_by_gls.t.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 7, width = 8.5, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  # 
  # pA<-ggplot(data=Daily.l.G)+
  #   geom_point(alpha = 0.5, size = 1.5, aes(gls.t,summary.vegc))+
  #   theme_bw()+
  #   facet_grid(num.mprcp~simset)+theme(text= element_text(size=15))+
  #   labs(shape = "Assimilation Type", color = "Assimilation Type", x = "Daily Stomatal Conductance", y = "Daily Assimilation")
  # plotName<-"Daily.drydown.l.G.vegc_by_gls.t.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 7, width = 8.5, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  
  

}

#psi_ leaf by yearday 
if(TRUE){

  pA<- ggplot(data=Annualized)+
    geom_point(size = 2, aes(x=simyearday,y=psi_leaf.m, color = as.factor(simset)))+
    geom_errorbar(alpha = 0.33, aes(x=simyearday,ymin=psi_leaf.m-psi_leaf.sd, ymax=psi_leaf.m+psi_leaf.sd, color = as.factor(simset)))+
    geom_hline(yintercept=-2.5, linetype=2, color="black", size=1)+
    scale_color_manual(values = cbbPalette[c(5,6)])+facet_grid(species~gcm)+
    theme_bw()+theme(text= element_text(size=25), strip.background = element_blank())
  plotName<-"psi_leaf_by_yearday1.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 8.5, width = 11, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  groups(Annualized)
  Annualized.wiltingpoint.count<-Annualized %>% filter(psi_leaf.m>-2.5) %>% summarize(count=n()) %>% pivot_wider(id_cols = c("gcm","species"), names_from = simset, values_from = count)
  Annualized.wiltingpoint.count[which(Annualized.wiltingpoint.count$species %in% c("ebfArp","ebfMA") & Annualized.wiltingpoint.count$gcm =="546_28"),]$base<-284
  Annualized.wiltingpoint.count[which(Annualized.wiltingpoint.count$species %in% c("ebfArp","ebfMA") & Annualized.wiltingpoint.count$gcm =="182_28"),]$base<-247
  Annualized.wiltingpoint.count <- Annualized.wiltingpoint.count %>% filter(species != "ebf")
  Annualized.wiltingpoint.count$gday_diff <- Annualized.wiltingpoint.count$stem_leafsoil - Annualized.wiltingpoint.count$base
  write_csv(Annualized.wiltingpoint.count, "nonwilting_days.csv")
  #for drawing lines on the above plot
# air seeding
# Group     psi_p_50  sd_p pit_area
# * <chr>        <dbl> <dbl>    <dbl>
# 1 Cordaites    -3.43 0.869   0.0117
# 2 Lycopsid     -1.28 1.20    0.0266
# 3 Psaronius    -2.07 1.45    0.0211
# 
# implosion
# Group.x      psi_p50 psi_refill
# * <chr>          <dbl>      <dbl>
#   1 Cordaitalean   25.2    -0.00888
# 2 Lycopsid        6.30   -0.00457
# 3 Medullosan     20.8    -0.00387
# 4 Psaronius       2.57   -0.0044 *manually updated for b=75 um
# 5 Sphenophyte    12.5    -0.00545

# plotName<-"Daily.100.G.soilw"
# ggsave(plotName, plot = pA, device = "png", path = "~/Dropbox/Brent_Spiner/results_01262018/",
#        scale = 1, height = 1.4493, width = 2.7096, units = c("in"),
#        dpi = 300, limitsize = TRUE)


}

#psi soil by yearday
if(TRUE){
  cbbPalette <- c("#E69F00", "#009E73", "#0072B2", "#56B4E9","#000000", "#D55E00", "#CC79A7")
  pA <- ggplot(data=Annualized[which(Annualized$simset=="base"),],aes(x=simyearday, y=psi.m, color=gcm))+
  geom_point()+
  geom_errorbar(alpha=0.5, aes(x=simyearday, ymin=psi.m-psi.sd, ymax=psi.m+psi.sd ))+
  scale_color_manual(values = cbbPalette)+facet_grid(species~., scales = 'free_y')+theme_minimal()
plotName<-"Annualized.base.psi_by_yearday.png"
ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
       scale = 1, height = 5, width = 7, units = c("in"),
       dpi = 300, limitsize = TRUE)

pA <- ggplot(data=Annualized[which(Annualized$simset=="stem_leafsoil"),],aes(x=simyearday, y=psi.m, color=gcm))+
  geom_point()+
  geom_errorbar(alpha=0.5,aes(x=simyearday, ymin=psi.m-psi.sd, ymax=psi.m+psi.sd ))+
  scale_color_manual(values = cbbPalette)+facet_grid(species~., scales = 'free_y')+theme_minimal()
plotName<-"Annualized.stem.psi_by_yearday.png"
ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
       scale = 1, height = 5, width = 7, units = c("in"),
       dpi = 300, limitsize = TRUE)
}

#gls.t by yearday 
if(TRUE){
  # facet_grid(.~simset)+
  # 
  # pA <- ggplot(data=Daily.stem.c,aes(x=simyearday, y=gls.t))+
  #   geom_smooth(color = "black",se = TRUE, level = 0.99)+
  #   theme_minimal()+theme(text= element_text(size=15))
  # plotName<-"Daily.100.G.c.gls_by_yearday.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 5, width = 7, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  # 
  # pA <- ggplot(data=Daily.G.l,aes(x=simyearday, y=gls.t))+
  #   geom_smooth(color = "black",aes(x=simyearday, y=gls.t),se = TRUE)+
  #   facet_grid(.~simset)+theme_minimal()+theme(text= element_text(size=15))
  # plotName<-"Daily.100.G.l.gls_by_yearday.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 5, width = 7, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  
  pA <- ggplot(data=Annualized[which(Annualized$simset=="base"),],aes(x=simyearday, y=glst.m, color=gcm))+
    geom_point()+
    geom_errorbar(alpha=0.5, aes(x=simyearday, ymin=glst.m-glst.sd, ymax=glst.m+glst.sd ))+
    scale_color_manual(values = cbbPalette)+facet_grid(species~., scales = 'free_y')+theme_minimal()
  plotName<-"Annualized.base.gls_by_yearday.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA <- ggplot(data=Annualized[which(Annualized$simset=="stem_leafsoil"),],aes(x=simyearday, y=glst.m, color=gcm))+
    geom_point()+
    geom_errorbar(alpha=0.5,aes(x=simyearday, ymin=glst.m-glst.sd, ymax=glst.m+glst.sd ))+
    scale_color_manual(values = cbbPalette)+facet_grid(species~., scales = 'free_y')+theme_minimal()
  plotName<-"Annualized.stem.gls_by_yearday.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  cbbPalette <- c("#E69F00", "#009E73", "#0072B2", "#56B4E9","#000000", "#D55E00", "#CC79A7")
  
  
  pA <- ggplot(data=Daily[which(Daily$fossil=="fossil"),],aes(x=as.factor(simyearmonth), y=gls.t, color=simtype))+
    geom_boxplot()+facet_grid(co2~.)+scale_color_manual(values = c("#000000", "#D55E00"))+
    theme_minimal()
  plotName<-"gls_by_month.png"
  ggsave(plotName, plot = pA, device = "png", path = "/Users/willmatthaeus/Dropbox/Brent_Spiner/CESM2BGC/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA <- ggplot(data=Daily[which(Daily$fossil=="fossil"),],aes(x=as.factor(lon), y=gls.t, color=simtype))+
    geom_boxplot()+facet_grid(co2~species)+scale_color_manual(values = c("#000000", "#D55E00"))+
    theme_minimal()
  plotName<-"gls_by_lon_sp.png"
  ggsave(plotName, plot = pA, device = "png", path = "/Users/willmatthaeus/Dropbox/Brent_Spiner/CESM2BGC/",
         scale = 1, height = 6, width = 8, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  # unique(Daily.100.noRev$simset)
 # pA<-ggplot(data=Daily.100.noRev, aes(x=simyearday, y=gls.t))+
 #    geom_point(
 #      size=0.5 , alpha = 0.5, aes(color = species2, shape = species2))+
 #    geom_smooth(color = 'black')+
 #    facet_grid(climate~simset)+theme_minimal()+theme(text= element_text(size=15))
 #  plotName<-"Daily.100.glsT_by_yearday.png"
 #  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
 #         scale = 1, height = 5, width = 9, units = c("in"),
 #         dpi = 300, limitsize = TRUE)

  
}

#A.t by yearday 
if(TRUE){
  
  pA <- ggplot(data=Annualized[which(Annualized$simset=="base"),],aes(x=simyearday, y=At.m, color=gcm))+
    geom_point()+
    geom_errorbar(alpha=0.5, aes(x=simyearday, ymin=At.m-At.sd, ymax=At.m+At.sd ))+
    scale_color_manual(values = cbbPalette)+facet_grid(species~., scales = 'free_y')+theme_minimal()
  plotName<-"Annualized.base.A_by_yearday.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA <- ggplot(data=Annualized[which(Annualized$simset=="stem_leafsoil"),],aes(x=simyearday, y=At.m, color=gcm))+
    geom_point()+
    geom_errorbar(alpha=0.5,aes(x=simyearday, ymin=At.m-At.sd, ymax=At.m+At.sd ))+
    scale_color_manual(values = cbbPalette)+facet_grid(species~., scales = 'free_y')+theme_minimal()
  plotName<-"Annualized.stem.A_by_yearday.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  # pA<-ggplot(data=Daily.G.c,aes(x=simyearday, y=A.t))+
  #   geom_smooth()+
  #   facet_grid(.~simset)+theme_minimal()+theme(text= element_text(size=15),se = TRUE)
  # plotName<-"Daily.100.G.c.A.t_by_yearday.smooths.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 5, width = 7, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  # 
  # pA<-ggplot(data=Daily.G.l,aes(x=simyearday, y=A.t))+
  #    geom_smooth()+
  #    facet_grid(.~simset)+theme_minimal()+theme(text= element_text(size=15),se = TRUE)
  #  plotName<-"Daily.100.G.l.A.t_by_yearday.smooths.png"
  #  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #         scale = 1, height = 5, width = 7, units = c("in"),
  #         dpi = 300, limitsize = TRUE)
  # 
  # 
  #  pA<-ggplot(data=Daily.100.noRev, aes(x=simyearday, y=A.t))+
  #    geom_point(
  #      size=0.5 , alpha = 0.5, aes(color = species2, shape = species2))+
  #    geom_smooth(color = 'black')+
  #    facet_grid(climate~simset)+theme_minimal()+theme(text= element_text(size=15))
  #  plotName<-"Daily.100.At_by_yearday.png"
  #  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #         scale = 1, height = 5, width = 9, units = c("in"),
  #         dpi = 300, limitsize = TRUE)
   
  
}

#lai by yearday  or simday
if(TRUE){
  cbbPalette <- c("#E69F00", "#009E73", "#0072B2", "#56B4E9","#000000", "#D55E00", "#CC79A7")
  
  pA <- ggplot(data=Annualized[which(Annualized$simset=="base"),],aes(x=simyearday, y=lai.m, color=gcm))+
    geom_point()+
    geom_errorbar(alpha=0.5, aes(x=simyearday, ymin=lai.m-lai.sd, ymax=lai.m+lai.sd ))+
    scale_color_manual(values = cbbPalette)+facet_grid(species~., scales = 'free_y')+theme_minimal()
  plotName<-"Annualized.base.lai_by_yearday.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA <- ggplot(data=Annualized[which(Annualized$simset=="stem_leafsoil"),],aes(x=simyearday, y=lai.m, color=gcm))+
    geom_point()+
    geom_errorbar(alpha=0.5,aes(x=simyearday, ymin=lai.m-lai.sd, ymax=lai.m+lai.sd ))+
    scale_color_manual(values = cbbPalette)+facet_grid(species~., scales = 'free_y')+theme_minimal()
  plotName<-"Annualized.stem.lai_by_yearday.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  # pA<-ggplot(data=Daily.base.c)+
  #   geom_point(
  #     size=0.5 , alpha = 0.5, aes(x=simday, y=epv.proj_lai, color=as.factor(species)))+
  #   geom_smooth(color = "black",aes(x=simday, y=epv.proj_lai),method = "lm")+
  #   facet_grid(climate~.)+theme_minimal()+theme(text= element_text(size=15))
  # plotName<-"Daily.100.c.base.lai_simday.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 5, width = 10, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  # 
  # pA<-ggplot(data=Daily.100.noRev.cltf)+
  #   geom_point(
  #     size=0.5 , alpha = 0.5, aes(x=simyearday, y=epv.proj_lai, color=as.factor(species)))+
  #   geom_smooth(color = "black",aes(x=simyearday, y=epv.proj_lai))+
  #   facet_grid(climate~simset)+theme_minimal()+theme(text= element_text(size=15))
  # plotName<-"Daily.100.noRev.cltf.epv.proj_lai_by_yearday.summary.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 5, width = 7, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  # 
  # pA<-  ggplot(data=Daily.100.noRev.cltf[which(Daily.100.noRev.cltf$species!="lepidodendron"),])+
  #   geom_point(
  #     size=0.5 , alpha = 0.5, aes(x=simyearday, y=epv.proj_lai, color=as.factor(species)))+
  #   geom_smooth(color = "black",aes(x=simyearday, y=epv.proj_lai))+
  #   facet_grid(climate~simset)+theme_minimal()+theme(text= element_text(size=15))
  # plotName<-"Daily.100.noRev.ctf.epv.proj_lai_by_yearday.summary.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 5, width = 7, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  # 
  # Daily.100.noRev.c <- Daily.100.noRev.cltf %>% filter(species %in% "cordaites")
  # 
  # pA<-ggplot(data=Daily.100.noRev.c)+
  #   geom_point(
  #     size=0.1 , color = 'grey', alpha = 1, aes(x=simyearday, y=epv.proj_lai))+
  #   geom_smooth(color = 'black', aes(x=simyearday, y=epv.proj_lai))+
  #   geom_hline(yintercept=-3.43, linetype=2, color="red", size=0.5)+
  #   coord_cartesian(xlim = c(1, 365), ylim=c(-4.5,0))+
  #   facet_grid(climate~simset)+theme_minimal()+theme(text= element_text(size=15))
  # plotName<-"Daily.100.noRev.c.epv.proj_laiT_by_yearday.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 5, width = 6, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  # 
  # Daily.100.noRev.l <- Daily.100.noRev.cltf %>% filter(species %in% "lepidodendron")
  # pA<-ggplot(data=Daily.100.noRev.l)+
  #   geom_point(
  #     size=0.1 , color = 'grey', alpha = 1, aes(x=simyearday, y=epv.proj_lai))+
  #   geom_smooth(color = 'black', aes(x=simyearday, y=epv.proj_lai))+
  #   geom_hline(yintercept=-1.28, linetype=2, color="red", size=0.5)+
  #   coord_cartesian(xlim = c(1, 365), ylim=c(-4.5,0))+
  #   facet_grid(climate~simset)+theme_minimal()+theme(text= element_text(size=15))
  # plotName<-"Daily.100.noRev.l.epv.proj_lai_by_yearday.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 5, width = 6, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  # 
  # Daily.100.noRev.tf <- Daily.100.noRev.cltf %>% filter(species %in% "treefern")
  # pA<-ggplot(data=Daily.100.noRev.tf)+
  #   geom_point(
  #     size=0.1 , color = 'grey', alpha = 1, aes(x=simyearday, y=epv.proj_lai))+
  #   geom_smooth(color = 'black', aes(x=simyearday, y=epv.proj_lai))+
  #   geom_hline(yintercept=-2.07, linetype=2, color="red", size=0.5)+
  #   coord_cartesian(xlim = c(1, 365), ylim=c(-4.5,0))+
  #   facet_grid(climate~simset)+theme_minimal()+theme(text= element_text(size=15))
  # plotName<-"Daily.100.noRev.tf.epv.proj_lai_by_yearday.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 5, width = 6, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  
  
}

#vegc by yearday
if(TRUE){
  cbbPalette <- c("#E69F00", "#009E73", "#0072B2", "#56B4E9","#000000", "#D55E00", "#CC79A7")
  
  pA <- ggplot(data=Annualized[which(Annualized$simset=="base"),],aes(x=simyearday, y=vegc.m, color=gcm))+
    geom_point()+
    geom_errorbar(alpha=0.5, aes(x=simyearday, ymin=vegc.m-vegc.sd, ymax=vegc.m+vegc.sd ))+
    scale_color_manual(values = cbbPalette)+facet_grid(species~., scales = 'free_y')+theme_minimal()
  plotName<-"Annualized.base.vegc_by_yearday.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA <- ggplot(data=Annualized[which(Annualized$simset=="stem_leafsoil"),]$vegc.m,aes(x=simyearday, y=vegc.m, color=gcm))+
    geom_point()+
    # geom_errorbar(alpha=0.5,aes(x=simyearday, ymin=vegc.m-vegc.sd, ymax=vegc.m+vegc.sd ))+
    scale_color_manual(values = cbbPalette)+facet_grid(species~gcm, scales = 'free_y')+theme_minimal()
  plotName<-"Annualized.stem.vegc_by_yearday.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 9, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA <- ggplot()+
    geom_point(data=Daily[which(Daily$simset=="stem_leafsoil" & Daily$num.mprcp==100),],aes(x=simyearday, y=cs.deadcrootc, color=gcm))+
    scale_color_manual(values = cbbPalette)+facet_grid(species~gcm, scales = 'free_y')+theme_minimal()
  plotName<-"Daily.stem.vegc_by_yearday.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 9, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
}



#swo by yearday
if(TRUE){
  pA<-ggplot(data=Daily.100.noRev.cltf)+
    geom_point(
      size=0.5 , alpha = 0.5, aes(x=simyearday, y=wf.soilw_outflow, color=as.factor(species)))+
    geom_smooth(color = "black",aes(x=simyearday, y=wf.soilw_outflow))+
    facet_grid(climate~simset)+theme_minimal()+theme(text= element_text(size=15))
  plotName<-"Daily.100.noRev.cltf.SWO_by_yearday.summary.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  Daily.100.noRev.c <- Daily.100.noRev.cltf %>% filter(species %in% "cordaites")
  pA<-ggplot(data=Daily.100.noRev.c)+
    geom_point(
      size=0.1 , color = 'grey', alpha = 1, aes(x=simyearday, y=wf.soilw_outlow))+
    geom_smooth(color = 'black', aes(x=simyearday, y=wf.soilw_outlow))+
    geom_hline(yintercept=-3.43, linetype=2, color="red", size=0.5)+
    coord_cartesian(xlim = c(1, 365), ylim=c(-4.5,0))+
    facet_grid(climate~simset)+theme_minimal()+theme(text= element_text(size=15))
  plotName<-"Daily.100.noRev.c.SWO_by_yearday.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 6, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  Daily.100.noRev.l <- Daily.100.noRev.cltf %>% filter(species %in% "lepidodendron")
  pA<-ggplot(data=Daily.100.noRev.l)+
    geom_point(
      size=0.1 , color = 'grey', alpha = 1, aes(x=simyearday, y=wf.soilw_outlow))+
    geom_smooth(color = 'black', aes(x=simyearday, y=wf.soilw_outlow))+
    geom_hline(yintercept=-1.28, linetype=2, color="red", size=0.5)+
    coord_cartesian(xlim = c(1, 365), ylim=c(-4.5,0))+
    facet_grid(climate~simset)+theme_minimal()+theme(text= element_text(size=15))
  plotName<-"Daily.100.noRev.l.SWO_by_yearday.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 6, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  Daily.100.noRev.tf <- Daily.100.noRev.cltf %>% filter(species %in% "treefern")
  pA<-ggplot(data=Daily.100.noRev.tf)+
    geom_point(
      size=0.1 , color = 'grey', alpha = 1, aes(x=simyearday, y=wf.soilw_outlow))+
    geom_smooth(color = 'black', aes(x=simyearday, y=wf.soilw_outlow))+
    geom_hline(yintercept=-2.07, linetype=2, color="red", size=0.5)+
    coord_cartesian(xlim = c(1, 365), ylim=c(-4.5,0))+
    facet_grid(climate~simset)+theme_minimal()+theme(text= element_text(size=15))
  plotName<-"Daily.100.noRev.tf.SWO_by_yearday.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 6, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=Daily.100.noRev, aes(x=simyearday, y=wf.soilw_outflow))+
    geom_point(
      size=0.5 , alpha = 0.5, aes(color = species2, shape = species2))+
    geom_smooth(color = 'black')+
    facet_grid(climate~simset)+theme_minimal()+theme(text= element_text(size=15))
  plotName<-"Daily.100.swo_by_yearday.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 9, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=Daily.100.noRev[which(Daily.100.noRev$gcm == "546_28"),], aes(x=simyearday, y=wf.soilw_outflow/metv.prcp))+
    geom_point(
      size=0.5 , alpha = 0.5, aes(color = species2, shape = species2))+
    geom_smooth(color = 'black')+
    facet_grid(species2~simset)+theme_minimal()+theme(text= element_text(size=15))
  plotName<-"Daily.100.I.ROR_by_yearday.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 9, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=Daily.100.noRev[which(Daily.100.noRev$gcm == "182_28"),], aes(x=simyearday, y=wf.soilw_outflow/metv.prcp))+
    geom_point(
      size=0.5 , alpha = 0.5, aes(color = species2, shape = species2))+
    geom_smooth(color = 'black')+
    facet_grid(species2~simset)+theme_minimal()+theme(text= element_text(size=15))
  plotName<-"Daily.100.G.ROR_by_yearday.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 9, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
}



#difference plots
if(TRUE){
##########################################################################
#difference plots: SWO by yearday, NetA ..., WUEi ...
##!! spread data around simset as a key?
diff.Daily<- Daily.100 %>% dplyr::select(simset,climate,num.mprcp,species,simday,gls.t,wf.soilw_outflow,A.t,metv.prcp) %>%
  pivot_wider(id_cols = c(simset,climate,num.mprcp,species,simday),names_from=simset,
              values_from = c(gls.t,wf.soilw_outflow,A.t,metv.prcp))  
diff.Daily <- diff.Daily %>% filter(species %in% c("cordaites", "macroneuropteris")) 

diff.Daily <- diff.Daily%>%group_by(species,climate)
diff.Daily$simyearday <- (diff.Daily$simday%%365)+1
# unique(diff.Daily$prcp_base==diff.Daily$prcp_stem)
# diff.Daily.100<-diff.Daily %>% filter(num.mprcp %in% c("100"))

diff.Daily$diff.sw.out <- diff.Daily$wf.soilw_outflow_stem_leafsoil - diff.Daily$wf.soilw_outflow_base
diff.Daily$diff.A.t <- diff.Daily$A.t_stem_leafsoil - diff.Daily$A.t_base
diff.Daily$diff.gls.t <- diff.Daily$gls.t_stem_leafsoil - diff.Daily$gls.t_base

#######swo, precip kg/m2
# hist(diff.Daily.100$prcp_stem[which(diff.Daily.100$climate=="182_28_100%")])
# hist(diff.Daily.100$prcp_stem[which(diff.Daily.100$climate=="546_28_100%")])
diff.Daily <- diff.Daily %>% mutate(base.mean.glst = mean(gls.t_base), base.mean.swo = mean(wf.soilw_outflow_base),
                                    base.mean.At = mean(A.t_base))
diff.Daily <- diff.Daily %>% mutate(p.diff.glst = diff.gls.t/base.mean.glst, p.diff.swo = diff.sw.out/base.mean.swo, 
                                    p.diff.At = diff.A.t/base.mean.At)

# m.p.182.c<-mean(diff.Daily$metv.prcp_stem_leafsoil[which(diff.Daily$climate=="182_28_100%" & diff.Daily$species == "cordaites")])
# m.p.546.c<-mean(diff.Daily$metv.prcp_stem_leafsoil[which(diff.Daily$climate=="546_28_100%" & diff.Daily$species == "cordaites")])
# m.p.182.m<-mean(diff.Daily$metv.prcp_stem_leafsoil[which(diff.Daily$climate=="182_28_100%" & diff.Daily$species == "macroneuropteris")])
# m.p.546.m<-mean(diff.Daily$metv.prcp_stem_leafsoil[which(diff.Daily$climate=="546_28_100%" & diff.Daily$species == "macroneuropteris")])
# 
# diff.Daily$p.diff.sw.out <--9999999999999
# diff.Daily$p.diff.sw.out[which(diff.Daily$climate=="182_28_100%")] <-
#   diff.Daily$diff.sw.out[which(diff.Daily$climate=="182_28_100%")]/m.p.182
# diff.Daily$p.diff.sw.out[which(diff.Daily$climate=="546_28_100%")] <-
#   diff.Daily$diff.sw.out[which(diff.Daily$climate=="546_28_100%")]/m.p.546
# which(diff.Daily$p.diff.sw.out == -9999999999999)
cbbPalette <- c("#E69F00", "#009E73", "#0072B2", "#56B4E9","#000000","#000000", "#D55E00", "#CC79A7")

pA<-ggplot(data=diff.Daily)+
  geom_point(size=1.0 , color = 'grey', alpha = 1.0, 
             aes(x=simyearday, y=diff.sw.out,shape = as.factor(species)))+
  geom_smooth(aes(x=simyearday, y=diff.sw.out, color = as.factor(species)))+
  scale_color_manual(values=cbbPalette[c(2,4)])+facet_grid(climate~.)+theme_minimal()
plotName<-"Daily.100.diff.SWO_by_yearday.png"
ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
       scale = 1, height = 5, width = 7, units = c("in"),
       dpi = 300, limitsize = TRUE)

pA<-ggplot(data=diff.Daily)+
  geom_point(size=1.5 , fill = 'grey', alpha = 0.5, color = "black", stroke = .5,
             aes(x=simyearday,y=p.diff.swo,shape = as.factor(species)))+
  geom_smooth(aes(x=simyearday, y=p.diff.swo, color = as.factor(species)))+
  scale_shape_manual(values=c(21,24))+scale_color_manual(values=cbbPalette[c(2,4)])+
  facet_grid(climate~.)+theme_minimal()
plotName<-"Daily.100.p.diff.swo_by_yearday.png"
ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
       scale = 1, height = 5, width = 7, units = c("in"),
       dpi = 300, limitsize = TRUE)

####### Assim
# hist(diff.Daily.100$A.t_base)
# hist(diff.Daily.100$A.t_stem)
# hist(diff.Daily.100$simyearday[which(abs(diff.Daily.100$A.t_base) <= 0.01)])
# hist(diff.Daily.100$simyearday[which(abs(diff.Daily.100$A.t_stem) <= 0.01)])
pA<-ggplot(data=diff.Daily)+
  geom_point(size=1.50 , color = 'grey', alpha = 1, 
             aes(x=simyearday,y=diff.gls.t,shape = as.factor(species)))+
  geom_smooth(aes(x=simyearday, y=diff.gls.t, color = as.factor(species)))+
  scale_shape_manual(values=c(21,24))+scale_color_manual(values=cbbPalette[c(2,4)])+
  facet_grid(climate~.)+theme_minimal()
plotName<-"Daily.100.diff.glst_by_yearday.png"
ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
       scale = 1, height = 5, width = 7, units = c("in"),
       dpi = 300, limitsize = TRUE)

pA<-ggplot(data=diff.Daily)+
  geom_point(size=1.5 , fill = 'grey', alpha = 0.5, color = "black", stroke = .5,
             aes(x=simyearday,y=p.diff.glst,shape = as.factor(species)))+
  geom_smooth(aes(x=simyearday, y=p.diff.glst, color = as.factor(species)))+
  scale_shape_manual(values=c(21,24))+scale_color_manual(values=cbbPalette[c(2,4)])+
  facet_grid(climate~.)+theme_minimal()
plotName<-"Daily.100.p.diff.glst_by_yearday.png"
ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
       scale = 1, height = 5, width = 7, units = c("in"),
       dpi = 300, limitsize = TRUE)

pA<-ggplot(data=diff.Daily)+
  geom_point(size=1.50 , color = 'grey', alpha = 1, 
             aes(x=simyearday,y=p.diff.swo,shape = as.factor(species)))+
  geom_smooth(aes(x=simyearday, y=p.diff.swo, color = as.factor(species)))+
  facet_grid(climate~.)+theme_minimal()
plotName<-"Daily.100.p.diff.swo_by_yearday.png"
ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
       scale = 1, height = 5, width = 7, units = c("in"),
       dpi = 300, limitsize = TRUE)

Daily$epv.gl_s_sun
#monthly
diff.Monthly<- Daily %>% dplyr::select(simtype,co2,lat,lon,species,simmonth,epv.gl_s_sun,epv.gl_s_sun,epv.proj_lai) %>%
  pivot_wider(id_cols = c(simtype,co2,lon,species,simmonth),names_from=simtype,
              values_from = c(epv.gl_s_sun,epv.gl_s_sun,epv.proj_lai))  

write_csv(diff.Monthly,file = "CESM_to_BGC_gls_lai_070522.csv")

unique(diff.Monthly$species)
diff.Monthly$fossil<-"fossil"
diff.Monthly[which(diff.Monthly$species %in% c("ebf","enf","ebfArp","ebfMA")),]$fossil<-"modern"
diff.Monthly<-diff.Monthly[which(diff.Monthly$fossil=="fossil"),]

diff.Monthly$simyearmonth <- (diff.Monthly$simmonth%%12)+1
diff.Monthly$gls.t.diff <- diff.Monthly$gls.t_paleo-diff.Monthly$gls.t_stem
which(is.na(diff.Monthly$gls.t.diff))
which(is.na(diff.Monthly$simyearmonth))
# scale_color_manual(values = c("#000000", "#D55E00"))+
pA <- ggplot(data=diff.Monthly,aes(x=as.factor(simyearmonth), y=gls.t.diff), color = 'blue')+
  geom_boxplot()+facet_grid(co2~.)+
  theme_minimal()
plotName<-"glsDiff_by_month.png"
ggsave(plotName, plot = pA, device = "png", path = "/Users/willmatthaeus/Dropbox/Brent_Spiner/CESM2BGC/",
       scale = 1, height = 5, width = 7, units = c("in"),
       dpi = 300, limitsize = TRUE)

CAT <- rgb(red = 0,green=0,blue=0,maxColorValue = 255)
lowEP <- rgb(red = 104,green=249,blue=41,maxColorValue = 255)
lowCP <- rgb(red = 141,green=32,blue=230,maxColorValue = 255)
CPM <-   rgb(red = 242,green=164,blue=30,maxColorValue = 255)
lowWP <- rgb(red = 0,green=7,blue=244,maxColorValue = 255)
ARM <- rgb(red = 223,green=26,blue=138,maxColorValue = 255)

soph_palette <- c(ARM,lowWP, CPM, lowCP, lowEP, CAT)

pA <- ggplot(data=diff.Monthly,aes(x=as.factor(lon), y=gls.t.diff, color=as.factor(lon)))+
  geom_boxplot()+facet_grid(co2~species)+scale_color_manual(values = soph_palette)+
  theme_minimal()
plotName<-"glsDiff_by_lon_species.png"
ggsave(plotName, plot = pA, device = "png", path = "/Users/willmatthaeus/Dropbox/Brent_Spiner/CESM2BGC/",
       scale = 1, height = 6, width = 8, units = c("in"),
       dpi = 300, limitsize = TRUE)


}

#drydown, ENDSIM / mean.sim
if(TRUE){
  ###dump everything but last simday, vegc, lai
  #spread to individual sims
  # Daily$sim_simday<-paste(sep="_",Daily$sim,Daily$simday)
  
  simMean.Daily <- Daily%>%filter(species %in% c('cordaites', 'lepidodendron','macroneuropteris','treefern')) %>%
    group_by(simset,gcm,num.mprcp,species)%>%
    summarize(mean.lai = mean(epv.proj_lai), mean.vegc = mean(summary.vegc), mean.At = mean(A.t),
              mean.gls.t = mean(gls.t), mean.mpsix = mean(m_psi_x), mean.psi_leaf = mean(epv.psi_leaf),
              mean.vwc = mean(vwc), mean.prcp = mean(metv.prcp))
  
  simMedian.Daily <- Daily%>%filter(species %in% c('cordaites', 'lepidodendron','macroneuropteris','treefern'))%>%group_by(simset,gcm,num.mprcp,species)%>%
    summarize(median.lai = median(epv.proj_lai), median.vegc = median(summary.vegc), median.At = median(A.t),
              median.gls.t = median(gls.t), median.mpsix = median(m_psi_x), median.psi_leaf = median(epv.psi_leaf),
              median.vwc = median(vwc), median.prcp = median(metv.prcp))
  
  
  
  #%>%filter(median.lai!=0)
  Dead<-FALSE
  for(g in unique(simMean.Daily$gcm)){
  for(simtype in unique(simMean.Daily$simset)){
  for(s in unique(simMean.Daily$species)){
    Dead<-FALSE
    for(n in seq(100, 10, -10))
    {
 
      if(Dead){
        simMean.Daily <-simMean.Daily[-c(which(simMean.Daily$species==s & 
                             simMean.Daily$num.mprcp==n &
                             simMean.Daily$gcm==g &
                             simMean.Daily$simset==simtype)),]
          # simMean.Daily[which(simMean.Daily$species==s & 
          #                            simMean.Daily$num.mprcp==n &
          #                            simMean.Daily$gcm==g &
          #                            simMean.Daily$simset==simtype),]$mean.lai<-0
          #     
          #     simMean.Daily[which(simMean.Daily$species==s & 
          #                           simMean.Daily$num.mprcp==n &
          #                           simMean.Daily$gcm==g &
          #                           simMean.Daily$simset==simtype),]$mean.gls.t<-0
          #     
          #     simMean.Daily[which(simMean.Daily$species==s & 
          #                           simMean.Daily$num.mprcp==n &
          #                           simMean.Daily$gcm==g &
          #                           simMean.Daily$simset==simtype),]$mean.At<-0
          #     
          #     simMean.Daily[which(simMean.Daily$species==s & 
          #                           simMean.Daily$num.mprcp==n &
          #                           simMean.Daily$gcm==g &
          #                           simMean.Daily$simset==simtype),]$mean.vegc<-0
      }else{if(simMean.Daily[which(simMean.Daily$species==s & 
                             simMean.Daily$num.mprcp==n &
                             simMean.Daily$gcm==g &
                             simMean.Daily$simset==simtype),]$mean.lai==0){Dead<-TRUE}
      }    
    }
  }
  }
  }
  #simMean
  if(TRUE){
    cbbPalette <- c("#E69F00", "#009E73", "#0072B2", "#56B4E9","#000000", "#D55E00", "#CC79A7")
    # unique(simMean.Daily$species)
  pA<-ggplot(data=simMean.Daily)+
    geom_point(size=1.0 , 
               aes(x=as.numeric(num.mprcp), y=mean.prcp,
                   color=as.factor(species)))+
    geom_line(size=1.0 , 
              aes(x=as.numeric(num.mprcp), y=mean.prcp,
                  color=as.factor(species)))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm)+theme_minimal()
  plotName<-"simave.drydown.prcp.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  #
  pA<-ggplot(data=simMean.Daily)+
    geom_point(size=1.0 , 
               aes(x=as.numeric(num.mprcp), y=mean.lai,
                   color=as.factor(species)))+
    geom_line(size=1.0 , 
               aes(x=as.numeric(num.mprcp), y=mean.lai,
                   color=as.factor(species)))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm)+theme_minimal()
  plotName<-"simave.drydown.lai.nofilter.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  #
  pA<-ggplot(data=simMean.Daily)+
    geom_point(size=1.0 , 
               aes(x=as.numeric(num.mprcp), y=mean.vegc,
                   color=as.factor(species)))+
    geom_line(size=1.0 , 
              aes(x=as.numeric(num.mprcp), y=mean.vegc,
                  color=as.factor(species)))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm)+theme_minimal()
  plotName<-"simave.drydown.vegc.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  #
  pA<-ggplot(data=simMean.Daily)+
    geom_point(size=1.0 , 
               aes(x=as.numeric(num.mprcp), y=mean.At,
                   color=as.factor(species)))+
    geom_line(size=1.0 , 
              aes(x=as.numeric(num.mprcp), y=mean.At,
                  color=as.factor(species)))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm)+theme_minimal()
  plotName<-"simave.drydown.At.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  #
  pA<-ggplot(data=simMean.Daily)+
    geom_point(size=1.0 , 
               aes(x=as.numeric(num.mprcp), y=mean.gls.t,
                   color=as.factor(species)))+
    geom_line(size=1.0 , 
              aes(x=as.numeric(num.mprcp), y=mean.gls.t,
                  color=as.factor(species)))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm, scales = "free_y")+theme_minimal()
  plotName<-"simave.drydown.glst.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  #
  
  pA<-ggplot(data=simMean.Daily)+
    geom_point(size=1.0 , 
               aes(x=as.numeric(num.mprcp), y=mean.mpsix,
                   color=as.factor(species)))+
    geom_line(size=1.0 , 
              aes(x=as.numeric(num.mprcp), y=mean.mpsix,
                  color=as.factor(species)))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm, scales = "free_y")+theme_minimal()
  plotName<-"simave.drydown.mpsix.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  #
  pA<-ggplot(data=simMean.Daily)+
    geom_point(size=1.0 , 
               aes(x=as.numeric(num.mprcp), y=mean.psi_leaf,
                   color=as.factor(species)))+
    geom_line(size=1.0 , 
              aes(x=as.numeric(num.mprcp), y=mean.psi_leaf,
                  color=as.factor(species)))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm, scales = "free_y")+theme_minimal()
  plotName<-"simave.drydown.psi_leaf.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  #
  pA<-ggplot(data=simMean.Daily)+
    geom_point(size=1.0 , 
               aes(x=as.numeric(num.mprcp), y=mean.vwc,
                   color=as.factor(species)))+
    geom_line(size=1.0 , 
              aes(x=as.numeric(num.mprcp), y=mean.vwc,
                  color=as.factor(species)))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm, scales = "free_y")+theme_minimal()
  plotName<-"simave.drydown.vwc.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  #
  }
  
  #simMedian
  if(TRUE){
    pA<-ggplot(data=simMedian.Daily)+
      geom_point(size=1.0 , 
                 aes(x=as.numeric(num.mprcp), y=median.prcp,
                     shape = as.factor(species),color=as.factor(species)))+
      geom_line(size=1.0 , 
                aes(x=as.numeric(num.mprcp), y=median.prcp,
                    color=as.factor(species)))+
      scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
      facet_grid(simset~gcm)+theme_minimal()
    plotName<-"simmed.drydown.prcp.png"
    ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
           scale = 1, height = 5, width = 7, units = c("in"),
           dpi = 300, limitsize = TRUE)
    #
    pA<-ggplot(data=simMedian.Daily)+
      geom_point(size=1.0 , 
                 aes(x=as.numeric(num.mprcp), y=median.lai,
                     shape = as.factor(species),color=as.factor(species)))+
      geom_line(size=1.0 , 
                aes(x=as.numeric(num.mprcp), y=median.lai,
                    color=as.factor(species)))+
      scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
      facet_grid(simset~gcm)+theme_minimal()
    plotName<-"simmed.drydown.lai.png"
    ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
           scale = 1, height = 5, width = 7, units = c("in"),
           dpi = 300, limitsize = TRUE)
    #
    pA<-ggplot(data=simMedian.Daily)+
      geom_point(size=1.0 , 
                 aes(x=as.numeric(num.mprcp), y=median.vegc,
                     shape = as.factor(species),color=as.factor(species)))+
      geom_line(size=1.0 , 
                aes(x=as.numeric(num.mprcp), y=median.vegc,
                    color=as.factor(species)))+
      scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
      facet_grid(simset~gcm)+theme_minimal()
    plotName<-"simmed.drydown.vegc.png"
    ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
           scale = 1, height = 5, width = 7, units = c("in"),
           dpi = 300, limitsize = TRUE)
    #
    pA<-ggplot(data=simMedian.Daily.noRev.cltf)+
      geom_point(size=1.0 , 
                 aes(x=as.numeric(num.mprcp), y=median.At,
                     shape = as.factor(species),color=as.factor(species)))+
      geom_line(size=1.0 , 
                aes(x=as.numeric(num.mprcp), y=median.At,
                    color=as.factor(species)))+
      scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
      facet_grid(simset~gcm)+theme_minimal()
    plotName<-"simave.drydown.noRev.cltf.At.med.png"
    ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
           scale = 1, height = 5, width = 7, units = c("in"),
           dpi = 300, limitsize = TRUE)
    #
    pA<-ggplot(data=simMedian.Daily.noRev.cltf)+
      geom_point(size=1.0 , 
                 aes(x=as.numeric(num.mprcp), y=median.gls.t,
                     shape = as.factor(species),color=as.factor(species)))+
      geom_line(size=1.0 , 
                aes(x=as.numeric(num.mprcp), y=median.gls.t,
                    color=as.factor(species)))+
      scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
      facet_grid(simset~gcm, scales = "free_y")+theme_minimal()
    plotName<-"simave.drydown.noRev.cltf.glst.med.png"
    ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
           scale = 1, height = 5, width = 7, units = c("in"),
           dpi = 300, limitsize = TRUE)
    #
    
    pA<-ggplot(data=simMedian.Daily.noRev.cltf)+
      geom_point(size=1.0 , 
                 aes(x=as.numeric(num.mprcp), y=median.mpsix,
                     shape = as.factor(species),color=as.factor(species)))+
      geom_line(size=1.0 , 
                aes(x=as.numeric(num.mprcp), y=median.mpsix,
                    color=as.factor(species)))+
      scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
      facet_grid(simset~gcm, scales = "free_y")+theme_minimal()
    plotName<-"simave.drydown.noRev.cltf.mpsix.med.png"
    ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
           scale = 1, height = 5, width = 7, units = c("in"),
           dpi = 300, limitsize = TRUE)
    #
    pA<-ggplot(data=simMedian.Daily.noRev.cltf)+
      geom_point(size=1.0 , 
                 aes(x=as.numeric(num.mprcp), y=median.psi_leaf,
                     shape = as.factor(species),color=as.factor(species)))+
      geom_line(size=1.0 , 
                aes(x=as.numeric(num.mprcp), y=median.psi_leaf,
                    color=as.factor(species)))+
      scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
      facet_grid(simset~gcm, scales = "free_y")+theme_minimal()
    plotName<-"simave.drydown.noRev.cltf.psi_leaf.med.png"
    ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
           scale = 1, height = 5, width = 7, units = c("in"),
           dpi = 300, limitsize = TRUE)
    #
    pA<-ggplot(data=simMedian.Daily.noRev.cltf)+
      geom_point(size=1.0 , 
                 aes(x=as.numeric(num.mprcp), y=median.vwc,
                     shape = as.factor(species),color=as.factor(species)))+
      geom_line(size=1.0 , 
                aes(x=as.numeric(num.mprcp), y=median.vwc,
                    color=as.factor(species)))+
      scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
      facet_grid(simset~gcm, scales = "free_y")+theme_minimal()
    plotName<-"simave.drydown.noRev.cltf.vwc.med.png"
    ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
           scale = 1, height = 5, width = 7, units = c("in"),
           dpi = 300, limitsize = TRUE)
    #
  }
  
  #sim End
  if(TRUE){
  # unique(Daily.noRev.cltf$sim)
  maxDays <-Daily %>% group_by(name) %>% summarize(maxDay = max(simday))#18250
  unique(maxDays$maxDay)
  lastyear.Daily <- Daily %>% filter(simday<=18250 & simday>(18250-365))
  lastyear.Daily %>% filter(num.mprcp==100, simset=="stem_leafsoil") %>% group_by(species) %>%  summarise(max.lsw_ratio = max(lsw_ratio))
  
  lastyear.Daily.clmtf <- lastyear.Daily %>% 
    filter(species %in% c("cordaites","lepidodendron","macroneuropteris","treefern"), simset != "stemRev")
  
  lastyear.Daily.clmtf$lai_ls_ratio <- lastyear.Daily.clmtf$cs.leafc/lastyear.Daily.clmtf$cs.livestemc
    
  
  #congo
  #check which precip sim matches
  #extrapolate runoff 4x10^6 km2 basin
  
  #appalachian basin
  #ask which precip sim and which assemblage
  #extrapolate runoff to 2.0957x10^6 km2
  
  hist(lastyear.Daily.clmtf$A.t, breaks = 100)
  
  endSim <- lastyear.Daily.clmtf %>% group_by(gcm, simset, species, num.mprcp) %>%
        summarize(lai.m = mean(epv.proj_lai), lai.sd = sd(epv.proj_lai),
                  vc.m = mean(vegc), vc.sd = sd(vegc),
                  lvc.m = mean(live.vegc), lvc.sd = sd(live.vegc),
                  lvstemc.m = mean(cs.livestemc), lvstemc.sd = sd(cs.livestemc),
                  lsr.m = mean(lai_ls_ratio), lsr.sd = sd(lai_ls_ratio),
                  wuei.m = mean(WUEi.t), wuei.sd = sd(WUEi.t),
                  wuee.m = mean(WUEe.t), wuee.sd = sd(WUEe.t),
                  nmin.m = mean(epv.daily_net_nmin), nmin.sd = sd(epv.daily_net_nmin),
                  nlea.m = mean(nf.sminn_leached), nlea.sd = sd(nf.sminn_leached),
                  at.m = mean(A.t), at.sd = sd(A.t), gs.m = mean(gls.t), gs.sd = sd(gls.t),
                  gs_days = sum(gls.t>0), at_days = sum(A.t>3))
  
  # groups(lastyear.Daily)
  lastyear.C.summary <- lastyear.Daily %>% filter(num.mprcp==100) %>% group_by(simset, gcm, species) %>%
    summarize(rootc.m.m = mean(rootc), stemc.m.m = mean(stemc), 
                                         leafc.m.m = mean(cs.leafc), lai.m.m = mean(epv.proj_lai))
  
  lastyear.C.summary$leafc.m.m <- lastyear.C.summary$leafc.m.m*1000
  # print(lastyear.C.summary,n = 22)
  
  lastyear.C.summary.w <- lastyear.C.summary%>%pivot_wider(id_cols = c(gcm, simset, species),
                                                               names_from = simset, values_from = c(rootc.m.m, stemc.m.m,leafc.m.m, lai.m.m),)
  
  ebf_546<-lastyear.C.summary.w[which(lastyear.C.summary.w$species=="ebf" & lastyear.C.summary.w$gcm =="546_28"),]
  ebf_182<-lastyear.C.summary.w[which(lastyear.C.summary.w$species=="ebf" & lastyear.C.summary.w$gcm =="182_28"),]
  
  lastyear.C.summary.w[which(lastyear.C.summary.w$species %in% c("ebfArp","ebfMA") & lastyear.C.summary.w$gcm =="546_28"),]$rootc.m.m_base<-43.0
  lastyear.C.summary.w[which(lastyear.C.summary.w$species %in% c("ebfArp","ebfMA") & lastyear.C.summary.w$gcm =="182_28"),]$rootc.m.m_base<-24.8

  lastyear.C.summary.w[which(lastyear.C.summary.w$species %in% c("ebfArp","ebfMA") & lastyear.C.summary.w$gcm =="546_28"),]$stemc.m.m_base<-143.0
  lastyear.C.summary.w[which(lastyear.C.summary.w$species %in% c("ebfArp","ebfMA") & lastyear.C.summary.w$gcm =="182_28"),]$stemc.m.m_base<-82.1
  
  lastyear.C.summary.w[which(lastyear.C.summary.w$species %in% c("ebfArp","ebfMA") & lastyear.C.summary.w$gcm =="546_28"),]$leafc.m.m_base<-173.0
  lastyear.C.summary.w[which(lastyear.C.summary.w$species %in% c("ebfArp","ebfMA") & lastyear.C.summary.w$gcm =="182_28"),]$leafc.m.m_base<-122.0
  
  lastyear.C.summary.w[which(lastyear.C.summary.w$species %in% c("ebfArp","ebfMA") & lastyear.C.summary.w$gcm =="546_28"),]$lai.m.m_base<-2.08
  lastyear.C.summary.w[which(lastyear.C.summary.w$species %in% c("ebfArp","ebfMA") & lastyear.C.summary.w$gcm =="182_28"),]$lai.m.m_base<-1.46
  
  lastyear.C.summary.w<-lastyear.C.summary.w %>% filter(species  %in% c("ebfArp","ebfMA","cordaites","treefern","macroneuropteris","lepidodendron"))
  
  lastyear.C.summary.w$d.root <- lastyear.C.summary.w$rootc.m.m_stem_leafsoil-lastyear.C.summary.w$rootc.m.m_base
  lastyear.C.summary.w$d.stem <- lastyear.C.summary.w$stemc.m.m_stem_leafsoil-lastyear.C.summary.w$stemc.m.m_base
  lastyear.C.summary.w$d.leaf <- lastyear.C.summary.w$leafc.m.m_stem_leafsoil-lastyear.C.summary.w$leafc.m.m_base
  lastyear.C.summary.w$d.lai <- lastyear.C.summary.w$lai.m.m_stem_leafsoil-lastyear.C.summary.w$lai.m.m_base
  
  C.table <-lastyear.C.summary.w %>% select(species, lai.m.m_base, d.lai, leafc.m.m_base, d.leaf, stemc.m.m_base, d.stem, rootc.m.m_base, d.root)
  write_csv(C.table,"C_table.csv")
  # pA<-ggplot(data=endSim.Daily)+
  #   geom_point(size=1.0 , 
  #              aes(x=as.numeric(num.mprcp), y=metv.prcp,
  #                  shape = as.factor(species),color=as.factor(species)))+
  #   geom_line(size=1.0 , 
  #             aes(x=as.numeric(num.mprcp), y=metv.prcp,
  #                 color=as.factor(species)))+
  #   scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
  #   facet_grid(simset~gcm)+theme_minimal()
  # plotName<-"simave.drydown.noRev.cltf.prcp.end.png"
  # ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
  #        scale = 1, height = 5, width = 7, units = c("in"),
  #        dpi = 300, limitsize = TRUE)
  #
  cbbPalette <- c("#E69F00", "#009E73", "#0072B2", "#56B4E9","#000000", "#D55E00", "#CC79A7")
  # unique(simMean.Daily$species)

  pA<-ggplot(data=endSim)+
    geom_point( 
       aes(x=num.mprcp, y=lai.m, color = species))+
          geom_errorbar(width = 0.5, aes(x=num.mprcp, 
                ymin=lai.m-lai.sd, ymax=lai.m+lai.sd, color = species))+
    geom_line(size=1.0 , 
              aes(x=num.mprcp, y=lai.m,
                  color=species))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm)+
    theme_bw()+theme(text= element_text(size=15), strip.background = element_blank(), legend.text = element_blank())
  plotName<-"lastyear.drydown.lai.clmtf.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width =7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  pA<-ggplot(data=endSim)+
    geom_point( 
      aes(x=num.mprcp, y=vc.m, color = species))+
    geom_errorbar(width = 0.5, aes(x=num.mprcp, 
                                   ymin=vc.m-vc.sd, ymax=vc.m+vc.sd, color = species))+
    geom_line(size=1.0 , 
              aes(x=num.mprcp, y=vc.m,
                  color=species))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm)+
    theme_minimal()
  plotName<-"lastyear.drydown.vc.clmtf.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=endSim)+
    geom_point( 
      aes(x=num.mprcp, y=lvc.m, color = species))+
    geom_errorbar(width = 0.5, aes(x=num.mprcp, 
                                   ymin=lvc.m-lvc.sd, ymax=lvc.m+lvc.sd, color = species))+
    geom_line(size=1.0 , 
              aes(x=num.mprcp, y=lvc.m,
                  color=species))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm)+
    theme_minimal()
  plotName<-"lastyear.drydown.lvc.clmtf.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=endSim)+
    geom_point( 
      aes(x=num.mprcp, y=lvstemc.m, color = species))+
    geom_errorbar(width = 0.5, aes(x=num.mprcp, 
                                   ymin=lvstemc.m-lvstemc.sd, ymax=lvstemc.m+lvstemc.sd, color = species))+
    geom_line(size=1.0 , 
              aes(x=num.mprcp, y=lvstemc.m,
                  color=species))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm)+
    theme_minimal()
  plotName<-"lastyear.drydown.lvstemc.clmtf.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=endSim)+
    geom_point( 
      aes(x=num.mprcp, y=lsr.m, color = species))+
    geom_errorbar(width = 0.5, aes(x=num.mprcp, 
                                   ymin=lsr.m-lsr.sd, ymax=lsr.m+lsr.sd, color = species))+
    geom_line(size=1.0 , 
              aes(x=num.mprcp, y=lsr.m,
                  color=species))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm)+
    theme_minimal()
  plotName<-"lastyear.drydown.lsr.clmtf.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=endSim)+
    geom_point( 
      aes(x=num.mprcp, y=wuei.m, color = species))+
    geom_errorbar(width = 0.5, aes(x=num.mprcp, 
                                   ymin=wuei.m-wuei.sd, ymax=wuei.m+wuei.sd, color = species))+
    geom_line(size=1.0 , 
              aes(x=num.mprcp, y=wuei.m,
                  color=species))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm)+
    theme_minimal()
  plotName<-"lastyear.drydown.wuei.clmtf.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=endSim)+
    geom_point( 
      aes(x=num.mprcp, y=wuei.m, color = species))+
    geom_errorbar(width = 0.5, aes(x=num.mprcp, 
                                   ymin=wuei.m-wuei.sd, ymax=wuei.m+wuei.sd, color = species))+
    geom_line(size=1.0 , 
              aes(x=num.mprcp, y=wuei.m,
                  color=species))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm)+
    theme_minimal()
  plotName<-"lastyear.drydown.wuei.clmtf.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=endSim)+
    geom_point( 
      aes(x=num.mprcp, y=wuee.m, color = species))+
    geom_errorbar(width = 0.5, aes(x=num.mprcp, 
                                   ymin=wuee.m-wuee.sd, ymax=wuee.m+wuee.sd, color = species))+
    geom_line(size=1.0 , 
              aes(x=num.mprcp, y=wuee.m,
                  color=species))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm)+
    theme_minimal()
  plotName<-"lastyear.drydown.wuee.clmtf.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=endSim)+
    geom_point( 
      aes(x=num.mprcp, y=nmin.m, color = species))+
    geom_errorbar(width = 0.5, aes(x=num.mprcp, 
                                   ymin=nmin.m-nmin.sd, ymax=nmin.m+nmin.sd, color = species))+
    geom_line(size=1.0 , 
              aes(x=num.mprcp, y=nmin.m,
                  color=species))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm)+
    theme_minimal()
  plotName<-"lastyear.drydown.nmin.clmtf.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=endSim)+
    geom_point( 
      aes(x=num.mprcp, y=nlea.m, color = species))+
    geom_errorbar(width = 0.5, aes(x=num.mprcp, 
                                   ymin=nlea.m-nlea.sd, ymax=nlea.m+nlea.sd, color = species))+
    geom_line(size=1.0 , 
              aes(x=num.mprcp, y=nlea.m,
                  color=species))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm)+
    theme_minimal()
  plotName<-"lastyear.drydown.nlea.clmtf.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=endSim)+
    geom_point( 
      aes(x=num.mprcp, y=at.m, color = species))+
    geom_errorbar(width = 0.5, aes(x=num.mprcp, 
                                   ymin=at.m-at.sd, ymax=at.m+at.sd, color = species))+
    geom_line(size=1.0 , 
              aes(x=num.mprcp, y=at.m,
                  color=species))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm)+
    theme_minimal()
  plotName<-"lastyear.drydown.at.clmtf.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=endSim)+
    geom_point( 
      aes(x=num.mprcp, y=gs.m, color = species))+
    geom_errorbar(width = 0.5, aes(x=num.mprcp, 
                                   ymin=gs.m-gs.sd, ymax=gs.m+gs.sd, color = species))+
    geom_line(size=1.0 , 
              aes(x=num.mprcp, y=gs.m,
                  color=species))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm)+
    theme_minimal()
  plotName<-"lastyear.drydown.gs.clmtf.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  pA<-ggplot(data=endSim)+
    geom_point( 
      aes(x=num.mprcp, y=at_days, color = species))+
      geom_line(size=1.0 , 
              aes(x=num.mprcp, y=at_days,
                  color=species))+
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_color_manual(values=cbbPalette)+facet_grid(simset~gcm)+
    theme_minimal()
  plotName<-"lastyear.drydown.atdays.clmtf.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  }

}


########met file summaries
################################################################
################################################################
if(TRUE){
  mtc43Head<- c("Yearday",
                "Tmax",
                "Tmin",
                "Tday",
                "Prcp",
                "VPD",
                "Srad",
                "Daylen")
  GLAC<-read_tsv("/Users/willmatthaeus/Dropbox/Brent_Spiner/stem_bgc/stem_model/metdata/GLAC_182_28_A.mtc43",skip = 1,col_names = FALSE)
  colnames(GLAC)<-mtc43Head
  summary(GLAC)

  Anz.GLAC  <- GLAC %>% group_by(Yearday) %>% summarize(tm.m = mean(Tday), tm.sd = sd(Tday),
                                  pr.m = mean(Prcp), pr.sd = sd(Prcp))
  Anz.GLAC$gcm <- "182_28"
                          
  IGLAC5<-read_tsv("/Users/willmatthaeus/Dropbox/Brent_Spiner/stem_bgc/stem_model/metdata/IGLAC_546_28_A.mtc43",skip = 1,col_names = FALSE)
  colnames(IGLAC5)<-mtc43Head
  summary(IGLAC5)
  
  Anz.IGLAC5  <- IGLAC5 %>% group_by(Yearday) %>% summarize(tm.m = mean(Tday), tm.sd = sd(Tday),
                                                        pr.m = mean(Prcp), pr.sd = sd(Prcp))
  Anz.IGLAC5$gcm <- "546_28"
  
  IGLAC6<-read_tsv("/Users/willmatthaeus/Dropbox/Brent_Spiner/stem_bgc/stem_model/metdata/IGLAC_600_21_A.mtc43",skip = 1,col_names = FALSE)
  colnames(IGLAC6)<-mtc43Head
  summary(IGLAC6)
  
  Anz.IGLAC6  <- IGLAC6 %>% group_by(Yearday) %>% summarize(tm.m = mean(Tday), tm.sd = sd(Tday),
                                                            pr.m = mean(Prcp), pr.sd = sd(Prcp))
  
  Anz.IGLAC6$gcm <- "600_21"
  
  mets <- rbind(Anz.GLAC,Anz.IGLAC5, Anz.IGLAC6)
  
  mets0 <- mets
  mets0$Yearday <- mets0$Yearday-365
  
  mets2 <- mets
  mets2$Yearday <- mets2$Yearday+365
  
  metsC <- rbind(mets0, mets, mets2)
  # scale_color_manual(values=cbbPalette)+
  
  pA<-ggplot(data=metsC)+coord_cartesian(xlim=c(0,365))+
    geom_point(aes(x=Yearday, y=tm.m, shape = gcm))+
    geom_smooth(se = FALSE, size=1, method = "gam", formula = y ~ s(x, k = 75), aes(x=Yearday, y=tm.m, color = gcm))+
    theme_minimal()+scale_shape_manual(values=c(16,17,0))+
    theme(text= element_text(size=20))
  plotName<-"TdayM.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  pA<-ggplot(data=metsC)+coord_cartesian(xlim=c(0,365))+
    geom_point(aes(x=Yearday, y=pr.m, shape = gcm))+
    geom_smooth(size=1 , method = "gam", formula = y ~ s(x, k = 75), aes(x=Yearday, y=pr.m, color = gcm))+
    theme_minimal()+scale_shape_manual(values=c(16,17,0))+
    theme(text= element_text(size=20))
  plotName<-"PrcpM.png"
  ggsave(plotName, plot = pA, device = "png", path = "~/Desktop/prelim stem figures/",
         scale = 1, height = 5, width = 7, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
}













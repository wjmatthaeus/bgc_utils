##data inputs generated using binary_grid_data_processing.R

pkgs <- c("readr", "ggplot2", "dplyr", "readxl",
          "tibble", "tidyr", "stringr", "raster",
          "tidyterra","scales","paletteer","plotbiomes",
          "viridis")
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE)

# install.packages("patchwork", "rio", "vegan", "gratia", "mgcViz", "DHARMa")
# devtools::install_github("valentinitnelav/plotbiomes")
# library(plotbiomes)
# whittaker_base_plot()

##workd for Ecosystem Engineers in Deep Time special issue paper
location<-"/path/to/EcosystemEngineers2024/"
setwd(location)

######Model update (stem2) run on SIM data from CESM. Statistics calculated per simulation on daily data using binary_grid_data_processing_kodiak.R

#data import and basic variables from .RDS created by "binary_grid_data..R" can be skipped by reading in EEDT24_data.csv
if(TRUE){

annual_means <- readRDS(file="path/to/binary_grid_data_processing/my.RDS")
#
 unique(annual_means$Genus)

annual_means <- annual_means %>%
  mutate(Tracheid_type = case_when(Genus == "arthropitys" ~ "cbp",
                                   Genus == "cordaites" ~ "cbp",
                                   Genus ==  "macroneuropteris" ~ "cbp",
                                   Genus == "arthropitysS" ~ "scalariform",
                                   Genus ==  "lepidodendron" ~ "scalariform",
                                   Genus ==  "spenophyllum" ~ "scalariform",
                                   Genus ==  "treefern" ~ "scalariform"))%>%
  rename(Model = `Empirical Model`)

annual_means <- annual_means %>%
  mutate(Group = case_when(Genus == "arthropitys" ~ "Sphenophyte (R)",
                                   Genus == "cordaites" ~ "Early Diverging Conifer",
                                   Genus ==  "macroneuropteris" ~ "Medullosan",
                                   Genus == "arthropitysS" ~ "Sphenophyte (M)",
                                   Genus ==  "lepidodendron" ~ "Arborescent Lycopsid",
                                   Genus ==  "spenophyllum" ~ "Sphenophyte (V)",
                                   Genus ==  "treefern" ~ "Tree Fern"))
 

annual_means$co2_f <- as.factor(annual_means$co2_mean)
#area calculations
    if(TRUE){
      unique(paste0(annual_means$Lat, annual_means$Lon,annual_means$Genus,annual_means$co2_f)%>%duplicated())
      testLoc<-unique(as.data.frame(cbind(annual_means$Lon,annual_means$Lat)))
      colnames(testLoc)<-c("Lon","Lat")
      testLoc$Area<-2*2*111*111.321*(cos((testLoc$Lat+1)*pi/180)+cos((testLoc$Lat-1)*pi/180))/2
      hist(testLoc$Area)
      total_area<-sum(testLoc$Area)
      formatC(total_area, format = "e", digits = 2)

      annual_means$Area<-2*2*111*111.321*(cos((annual_means$Lat+1)*pi/180)+cos((annual_means$Lat-1)*pi/180))/2
      hist(annual_means$Area)
    }

thisDF <- annual_means %>%
  dplyr::select(Lon, Lat, Model, co2_f, Genus) %>%
  filter(Model == "Eud", co2_f=="560",Genus=="cordaites") %>%
  dplyr::select(Lon, Lat)
thisDF$mask<-1
m<-thisDF%>%select(c("Lon","Lat","mask")) %>% 
  as_spatraster(xycols=1:2, crs = "" ,digits=1)

m_fort<-m%>%raster()%>%rasterToPolygons(dissolve=TRUE,n=16)%>%fortify()

#

all_genera<- unique(annual_means$Genus)
core_genera <- c("cordaites" ,    
                 "lepidodendron" ,   "macroneuropteris", 
                 "treefern" )
sphenophyte_genera <- c("arthropitys" ,     "arthropitysS"  , "spenophyllum")

core_means <- annual_means %>% filter(Genus %in% core_genera, tmin_min>-40)
sphen_means <- annual_means %>% filter(Genus %in% sphenophyte_genera)

}


#gs unit conversion
if(TRUE){
  df<-read_excel("EE24_WJM_Table1.xlsx")
  #(m/s)->(mmol/s) implicit per meter2
  #1mm/m2 = 1 L
  #1m/m2 = 1000L = 5.535×10^7 mmol (millimoles)
  #so... 5.535×10^7 mmol Water per 1 m
  tempC <- 20
  conv_factor_mmol <-  101325/(8.314510*(273.16+tempC)) *1000
  df$Gmax_mmol <-df$Gmax * conv_factor 
  
  df %>% select(Group, Gmax_mmol)
  }
  
#pbdb data for core groups
if(TRUE){

lep<-read_csv("pbdb_scotese_penn/pbdb_data_lep_scotese.csv",skip = 19,col_names = T)
lep_mega<- lep %>%  filter(!grepl("spora|spori",identified_name))
lep_mega$Group <-"Arborescent Lycopsid"
# unique(lycopsid_mega$identified_name) 
conifers<-read_csv("pbdb_scotese_penn/pbdb_data_pinophyta_scotese.csv",skip = 19,col_names = T)
# unique(conifers$identified_name) 
conifers_mega<- conifers %>%  filter(!grepl("spora|spori",identified_name))
conifers_mega$Group<-"Early Diverging Conifer"
#
medullosans <- read_csv("pbdb_scotese_penn/pbdb_data_medullosaceae_scotese.csv",skip = 19,col_names = T)
medullosans_mega<- medullosans %>%  filter(!grepl("spora|spori",identified_name))
# unique(medullosans_mega$identified_name) 
medullosans_mega$Group<-"Medullosan"
# medullosans_mega[which(medullosans_mega$paleolat)]
#
treeferns <- read_csv("pbdb_scotese_penn/pbdb_data_marattiales_scotese.csv",skip = 19,col_names = T)
treeferns_mega<- treeferns %>%  filter(!grepl("spora|spori",identified_name))
# unique(treeferns_mega$identified_name) 
treeferns_mega$Group<-"Tree Fern"

pbdb_data<-rbind( conifers_mega, medullosans_mega, treeferns_mega, lep_mega)
pbdb_data <- pbdb_data %>% unique()%>% select(Group,paleolng2, paleolat2)%>%rename(Lon=paleolng2,Lat=paleolat2)

#all data
pbdb <- read_csv("pbdb_scotese_penn/pbdb_data_all_scotese.csv",skip = 19,col_names = T)
pbdb_mega<- pbdb %>%  filter(!grepl("spora|spori",identified_name))%>%
  unique()%>%select(paleolng2,paleolat2)%>%rename(Lon=paleolng2,Lat=paleolat2)
}

#impact  of biosphere metrics (land area and soil water)
if(TRUE){
  # calculate global and grid-wise area, check distributions are the same
 
core_means %>% group_by(co2_f, Genus) %>% summarize(count = n())

modern_land_area<-144680000#from Ding et al. 2020: 10.1029/2020EF001618
modern_veg_land_area<-103900000
formatC(modern_veg_land_area, format = 'e', digits = 2)
paleo_per_mod_total_area <- total_area/modern_land_area

core_means%>%filter(epv.proj_lai_median>0.1)%>%
    group_by(co2_f,Genus)%>%
    summarize(
              vegArea=sum(Area)/1e6,
              vegArea_P_paleoLandArea=100*sum(Area)/total_area,
              vegArea_P_modVegArea=100*sum(Area)/modern_veg_land_area,
              # count=n()
              )


lep_sw <- core_means %>% filter(Genus == "lepidodendron") %>% 
  group_by(co2_f)%>%
  summarize(
    lepSwMean=mean(soilw_mean),
    lepSwoMean=mean(soilw_outflow_mean)
  )

sw_280 <- core_means %>% filter(co2_f == "280") %>% 
  group_by(Group)%>%
  summarize(
    lowSwMean=mean(soilw_mean),
    lowSwoMean=mean(soilw_outflow_mean)
  )

 sw_levels <- c("Arborescent Lycopsid" , "Tree Fern" ,"Early Diverging Conifer")
 sw_levels <- c("lepidodendron","treefern","cordaites")

         

core_means%>%filter(Genus %in% sw_levels)%>%
  pivot_wider(id_cols = c("co2_f","Lat","Lon"),names_from  = c("Genus"), values_from = c("soilw_outflow_mean","soilw_mean","Area"))%>%
  group_by(co2_f)%>%
  filter(soilw_outflow_mean_lepidodendron>0)%>%
  mutate(swo_tf_p_lep_local = 100*((soilw_outflow_mean_treefern-soilw_outflow_mean_lepidodendron)/soilw_outflow_mean_lepidodendron)*Area_lepidodendron,
         sw_tf_p_lep_local = 100*((soilw_mean_treefern-soilw_mean_lepidodendron)/soilw_mean_lepidodendron)*Area_lepidodendron,
         swo_cord_p_lep_local = 100*((soilw_outflow_mean_cordaites-soilw_outflow_mean_lepidodendron)/soilw_outflow_mean_lepidodendron)*Area_lepidodendron,
         sw_cord_p_lep_local = 100*((soilw_mean_cordaites-soilw_mean_lepidodendron)/soilw_mean_lepidodendron)*Area_lepidodendron,
         )%>%
  summarize(to_tf_swo = mean(swo_tf_p_lep_local,na.rm=T)/mean(Area_lepidodendron),
            to_tf_sw = mean(sw_tf_p_lep_local,na.rm=T)/mean(Area_lepidodendron),
            to_cord_swo = mean(swo_cord_p_lep_local,na.rm=T)/mean(Area_lepidodendron),
            to_cord_sw = mean(sw_cord_p_lep_local,na.rm=T)/mean(Area_lepidodendron))

  
core_means%>%filter(Genus %in% sw_levels)%>%
  pivot_wider(id_cols = c("Genus","Lat","Lon"),names_from  = c("co2_f"), values_from = c("soilw_outflow_mean","soilw_mean","Area"))%>%
  group_by(Genus)%>%
  mutate(Genus= factor(Genus, levels=sw_levels))%>%filter(soilw_outflow_mean_280>0)%>%
  mutate(swo_p_280_local = 100*((soilw_outflow_mean_560-soilw_outflow_mean_280)/soilw_outflow_mean_280)*Area_280,
         sw_p_280_local = 100*((soilw_mean_560-soilw_mean_280)/soilw_mean_280)*Area_280)%>%
  summarize(  sw_co2_2x = mean(sw_p_280_local,na.rm=T)/mean(Area_280),
              swo_co2_2x = mean(swo_p_280_local,na.rm=T)/mean(Area_280))



}

#vegetation cover by Genus
if(TRUE){

  hist(annual_means[which(annual_means$tmin_min< -40),]$tmin_sd, breaks=100)
 
  
  pA <- core_means%>%ggplot()+
    geom_tile(aes(x=Lon, y=Lat, fill=epv.proj_lai_median))+ 
    geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
    geom_point(data=pbdb_data, aes(x=Lon+160, y=Lat), color="red", size = 0.75)+
    scale_fill_gradientn(colors=viridis(6),limits=c(0.1,max(core_means$epv.proj_lai_median)))+theme_minimal()+ ##options A-H #values = v, breaks=l,
    theme(text = element_text(size=30), 
          legend.position = 'bottom', legend.key.width = unit(2,'in'),
          axis.text.x=element_blank())+#, legend.key.width = unit(1,'in')
    guides(fill=guide_colorbar(title=""))+
    labs(x="")+
    coord_map(projection="mollweide",ylim = c(-90,75), xlim = c(60,345))+
    scale_y_continuous(breaks = (-2:2) * 30)+
    facet_grid(Group~co2_f)
  
  plotName<-paste("core_CESM_St2_lai.png")
  ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
         scale = 1, height = 22, width = 17, units = c("in"),
         dpi = 300, limitsize = FALSE)

  
  ### sphen zero lai
  pA <- sphen_means %>% filter(epv.proj_lai_median==0)%>%
            ggplot()+
            geom_tile(aes(x=Lon, y=Lat, fill=epv.proj_lai_median), color = 'chocolate')+ #size=.25*((75-abs(Lat))/75)),,shape=15
            geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
            scale_fill_gradientn(colors=viridis(6),limits=c(0.1,max(core_means$epv.proj_lai_median)))+theme_minimal()+ ##options A-H #values = v, breaks=l,
            theme(text = element_text(size=30), 
                  legend.position = 'bottom', legend.key.width = unit(2,'in'),
                  axis.text.x=element_blank())+#, legend.key.width = unit(1,'in')
            guides(fill=guide_colorbar(title=""))+
            labs(x="")+
            coord_map(projection="mollweide",ylim = c(-90,75), xlim = c(60,345))+
            scale_y_continuous(breaks = (-2:2) * 30)+
            facet_grid(Genus~co2_f)
 
  plotName<-paste("sphen_CESM_St2_zeroLai.png")
  ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
         scale = 1, height = 22, width = 17, units = c("in"),
         dpi = 300, limitsize = FALSE)
  

  ## all fossils from pdbd
  pA <- annual_means%>%filter(Group=="Early Diverging Conifer")%>%ggplot()+
    geom_tile(aes(x=Lon, y=Lat, fill=prcp_mean*365/10))+ 
    geom_point(data=pbdb_mega, aes(x=Lon+160, y=Lat), color="green",size=0.75)+#size=.25*((75-abs(Lat))/75)),,shape=15
    geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
    scale_fill_gradientn(colors=plasma(6))+theme_minimal()+ ##options A-H #values = v, breaks=l,
    theme(text = element_text(size=30), 
          legend.position = 'bottom', legend.key.width = unit(1,'in'),
          axis.text.x=element_blank())+#, legend.key.width = unit(1,'in')
    coord_map(projection="mollweide",ylim = c(-90,75), xlim = c(60,345))+
    guides(fill=guide_colorbar(title="Simulated Mean\nAnnual Precipitation (cm)"))+
    labs(x="")+
    facet_grid(.~co2_f)
  
  plotName<-paste("pbdb_all_CESM_MAP.png")
  ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
         scale = 1, height = 7, width = 14, units = c("in"),
         dpi = 300, limitsize = FALSE)
  
  ####vegc
  v <- seq(0.0,1.0,0.1)
  q <- floor(unname(quantile(annual_means$vegc, probs = v)))
  # l <- seq(q[1],q[length(q)],1)
  o <- c(0.1,q[length(q)-1])
  
  clr<-c(viridis(length(l)))
  
  pA <- ggplot(data=annual_means[which(annual_means$vegc>0),])+
    geom_tile(aes(x=Lon, y=Lat, fill=vegc))+ #size=.25*((75-abs(Lat))/75)),,shape=15
    geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
    scale_fill_gradientn(colors=viridis(10), limits=o)+theme_minimal()+ #options A-H #values = v, breaks=l,
    theme(text = element_text(size=30), legend.position = 'bottom', legend.key.width = unit(2,'in'))+#, 
    coord_map(projection="mollweide",ylim = c(-90,75), xlim = c(60,345))+
    facet_grid(Genus~co2xmodern)
  
  plotName<-paste("CESM_St2_vegc.png")
  ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024//",
         scale = 1, height = 160, width = 15, units = c("in"),
         dpi = 300, limitsize = FALSE)
  
 
  for(g in genera){
   #unify scales
   # g <- "treefern"
   temp <- annual_means%>%filter(Genus==g)
   temp_cdf<-ecdf(temp$epv.proj_lai)
   temp_v<-unlist(lapply(l, temp_cdf))
   
   pA <- ggplot(data=temp)+
     geom_tile(aes(x=Lon, y=Lat, fill=epv.proj_lai))+ #size=.25*((75-abs(Lat))/75)),,shape=15
     geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
     scale_fill_gradientn(colors=clr,values = temp_v, breaks=l)+theme_minimal()+ #options A-H #, limits=o
     theme(text = element_text(size=30), legend.position = 'bottom', legend.key.width = unit(5,'cm'))+#
     coord_map(projection="mollweide",ylim = c(-90,75), xlim = c(60,345))+
     facet_wrap(.~co2xmodern, ncol=1)
   
   plotName<-paste(g,"_CESM_St2_lai.png")
   ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024//",
          scale = 1, height = 25, width = 15, units = c("in"),
          dpi = 300, limitsize = TRUE)
 }
  
}    

# soil water, runoff by community change
if(TRUE){
 
  #separate by taxa and co2
  #absolute water
  if(TRUE){
    runoff <- core_means %>% filter(Genus%in%c("lepidodendron","treefern"))%>%
      dplyr::mutate(
        m3_total_runoff = soilw_outflow_mean*365/1000*Area,
        m3_max_runoff = soilw_outflow_max/1000*Area)
    
   
    pA <- runoff%>%ggplot()+
      geom_tile(aes(x=Lon, y=Lat, fill=m3_max_runoff))+ 
      scale_fill_gradientn(colors=rev(plasma(6)))+
      geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
      theme_minimal()+ ##options A-H #values = v, breaks=l,
      theme(text = element_text(size=30), 
            legend.position = 'bottom', legend.key.width = unit(2,'in'),
            axis.text.x=element_blank())+#, legend.key.width = unit(1,'in')
      guides(fill=guide_colorbar(title=""))+
      labs(x="")+
      coord_map(projection="mollweide",ylim = c(-90,75), xlim = c(100,295))+
      scale_y_continuous(breaks = (-2:2) * 30)+
      facet_grid(Genus~co2_f)
    
    plotName<-paste("core_CESM_St2_m3_max_runoff.png")
    ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
           scale = 1, height = 22, width = 17, units = c("in"),
           dpi = 300, limitsize = FALSE)
    
  }
  
  
  #percent change sw and runoff
  if(TRUE){
  water_pc <- core_means %>% 
          tidyr::pivot_wider(id_cols = c("Lon","Lat","co2_f"),
                             names_from = c("Genus"), 
                             values_from = c(soilw_mean, soilw_outflow_mean
                                             ))%>%
          dplyr::mutate(PC_soilw_lepToTF = 100*(soilw_mean_treefern-soilw_mean_lepidodendron)/soilw_mean_lepidodendron,
                        PC_soilw_lepToCord = 100*(soilw_mean_cordaites-soilw_mean_lepidodendron)/soilw_mean_lepidodendron,
                        PC_runoff_lepToTF = 100*(soilw_outflow_mean_treefern-soilw_outflow_mean_lepidodendron)/soilw_outflow_mean_lepidodendron,
                        PC_runoff_lepToCord = 100*(soilw_outflow_mean_cordaites-soilw_outflow_mean_lepidodendron)/soilw_outflow_mean_lepidodendron
                        )
  
  pA <- cc%>%ggplot()+
    geom_tile(aes(x=Lon, y=Lat, fill=PC_runoff_lepToTF))+ 
    geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
    scale_fill_gradientn(colors=magma(6))+theme_minimal()+ ##options A-H #values = v, breaks=l,
    theme(text = element_text(size=30), 
          legend.position = 'bottom', legend.key.width = unit(2,'in'),
          axis.text.x=element_blank())+#, legend.key.width = unit(1,'in')
    guides(fill=guide_colorbar(title=""))+
    labs(x="")+
    coord_map(projection="mollweide",ylim = c(-90,75), xlim = c(60,345))+
    scale_y_continuous(breaks = (-2:2) * 30)+
    facet_grid(.~co2_f)
  
  
  plotName<-paste("core_CESM_St2_PC_runoff_lepToTF.png")
  ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
         scale = 1, height = 22, width = 17, units = c("in"),
         dpi = 300, limitsize = FALSE)
  }
  
  #absolute water change due to plants
  if(TRUE){
  water <- core_means %>% 
    tidyr::pivot_wider(id_cols = c("Lon","Lat","co2_f"),
                       names_from = c("Genus"), 
                       values_from = c(soilw_mean, soilw_outflow_max, soilw_outflow_mean,Area
                       ))%>%
    dplyr::mutate(
                  m3_total_runoff_lepToTF = (soilw_outflow_mean_treefern-soilw_outflow_mean_lepidodendron)*365/1000*Area_treefern,
                  m3_max_runoff_lepToTF = (soilw_outflow_max_treefern-soilw_outflow_max_lepidodendron)/1000*Area_treefern
    )
  hist(water$Area_cordaites)
  #water histograms
  if(TRUE){
  pA<-ggplot()+geom_histogram(data=water, aes(x=m3_total_runoff_lepToTF),fill="black",bins = 55)+
    geom_vline(xintercept =0,color="white")+theme_minimal()+
    theme(plot.background = element_blank(), panel.background = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), panel.grid = element_blank(),
          axis.line = element_blank(), panel.grid.major = element_blank())
  # c(min(water$m3_total_runoff_lepToTF, na.rm = F), 0, max(water$m3_total_runoff_lepToTF, na.rm = F)
  plotName<-paste("core_CESM_St2_m3_total_runoff_lepToTF_hist.png")
  ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
         scale = 1, height = 3, width = 5, units = c("in"),
         dpi = 300, limitsize = FALSE)
  
  pA<-ggplot()+geom_histogram(data=water, aes(x=m3_max_runoff_lepToTF),fill="black", bins = 55)+
    geom_vline(xintercept =0,color="white")+theme_minimal()+
    theme(plot.background = element_blank(), panel.background = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), panel.grid = element_blank(),
          axis.line = element_blank(), panel.grid.major = element_blank())
  plotName<-paste("core_CESM_St2_m3_max_runoff_lepToTF_hist.png")
  ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
         scale = 1, height = 3, width = 5, units = c("in"),
         dpi = 300, limitsize = FALSE)
  }

  #  runoff colors
  if(TRUE){
  colors<-c("#DE4968FF","#DE4968FF", "#FE9F6DFF", "#FCFDBFFF","white","lightsteelblue1","lightsteelblue3", "#0D0887FF")
  colors_ext<-c(colors[1],colors,colors[6])
  
  this_var <-water$m3_total_runoff_lepToTF
  this_sd<-sd(this_var, na.rm=TRUE)
  this_mean <- mean(this_var, na.rm=TRUE)
  this_max <- max(this_var, na.rm=T)
  this_min <- min(this_var, na.rm = T)
  val <- scales::rescale(x = c(this_min,this_mean-2*this_sd, this_mean-this_sd, 
                               this_mean,0, 
                               this_mean+this_sd,this_mean+2*this_sd,this_max),
                 to = c(0,1), 
                 from = c(this_min,this_max))
  }
 
  pA <- water%>%ggplot()+
    geom_tile(aes(x=Lon, y=Lat, fill=m3_total_runoff_lepToTF))+ 
    geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
    scale_fill_gradientn(colors = colors, values = val)+theme_minimal()+ ##options A-H #values = v, breaks=l,
    theme(text = element_text(size=40), 
          legend.position = 'bottom', legend.key.width = unit(2,'in'),
          axis.text.x=element_blank())+#, legend.key.width = unit(1,'in')
    guides(fill=guide_colorbar(title=""))+
    labs(x="")+
    coord_map(projection="mollweide",ylim = c(-90,75), xlim = c(100,295))+
    scale_y_continuous(breaks = (-2:2) * 30)+
    facet_grid(.~co2_f)
  
  plotName<-paste("core_CESM_St2_m3_total_runoff_lepToTF.png")
  ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
         scale = 1, height = 22, width = 17, units = c("in"),
         dpi = 300, limitsize = FALSE)
  
  }
  
  #absolute water change due to co2 Jay Lumen -- System; Dark Light -- A*S*Y*S and Dominik Schwarz
  if(TRUE){
    
    water_co2 <- core_means %>%filter(Genus %in% c("lepidodendron","treefern")) %>% 
      tidyr::pivot_wider(id_cols = c("Lon","Lat","Genus"),
                         names_from = c("co2_f"), 
                         values_from = c(soilw_mean, soilw_outflow_max, soilw_outflow_mean,Area
                         ))%>%
      dplyr::mutate(
        m3_total_runoff_560to280 = (soilw_outflow_mean_280-soilw_outflow_mean_560)*365/1000*Area_560,
        m3_max_runoff_560to280 = (soilw_outflow_max_280-soilw_outflow_max_560)/1000*Area_560
      )
   water_co2 <- water_co2 %>% filter(m3_total_runoff_560to280 > -100000)
   
   summary(water_co2)
  
    if(TRUE){
      pA<-ggplot()+geom_histogram(data=water_co2, aes(x=m3_total_runoff_560to280),fill="black",bins = 55)+
        geom_vline(xintercept =0,color="white")+theme_minimal()+
        theme(plot.background = element_blank(), panel.background = element_blank(),
              axis.title.y = element_blank(), axis.text.y = element_blank(),
              axis.ticks.y = element_blank(), panel.grid = element_blank(),
              axis.line = element_blank(), panel.grid.major = element_blank())
  
      plotName<-paste("core_CESM_St2_m3_total_runoff_560to280_hist.png")
      ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
             scale = 1, height = 3, width = 5, units = c("in"),
             dpi = 300, limitsize = FALSE)
      
      pA<-ggplot()+geom_histogram(data=water_co2, aes(x=m3_max_runoff_560to280),fill="black", bins = 55)+
        geom_vline(xintercept =0,color="white")+theme_minimal()+
        theme(plot.background = element_blank(), panel.background = element_blank(),
              axis.title.y = element_blank(), axis.text.y = element_blank(),
              axis.ticks.y = element_blank(), panel.grid = element_blank(),
              axis.line = element_blank(), panel.grid.major = element_blank())
      plotName<-paste("core_CESM_St2_m3_max_runoff_560to280_hist.png")
      ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
             scale = 1, height = 3, width = 5, units = c("in"),
             dpi = 300, limitsize = FALSE)
    }
    # runoff colors
    if(TRUE){
      colors<-c("#DE4968FF","#DE4968FF", "#FE9F6DFF", "#FCFDBFFF","white","lightsteelblue1","lightsteelblue3", "#0D0887FF")
      
      this_var <- water_co2$m3_max_runoff_560to280
      this_sd<-sd(this_var, na.rm=TRUE)
      this_mean <- mean(this_var, na.rm=TRUE)
      this_max <- max(this_var, na.rm=T)
      this_min <- min(this_var, na.rm = T)
      val <- scales::rescale(x = c(this_min,this_mean-2*this_sd, this_mean-this_sd, 
                                   this_mean,
                                   0, 
                                   this_mean+this_sd,this_mean+2*this_sd,this_max),
                             to = c(0,1), 
                             from = c(this_min,this_max))
    }
    
    pA <- water_co2%>%ggplot()+
      geom_tile(aes(x=Lon, y=Lat, fill=m3_max_runoff_560to280))+ 
      geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
      scale_fill_gradientn(colors=colors, values = val)+theme_minimal()+ ##options A-H #values = v, breaks=l,
      theme(text = element_text(size=40), 
            legend.position = 'bottom', legend.key.width = unit(2,'in'),
            axis.text.x=element_blank())+#, legend.key.width = unit(1,'in')
      guides(fill=guide_colorbar(title=""))+
      labs(x="")+
      coord_map(projection="mollweide",ylim = c(-90,75), xlim = c(100,295))+
      scale_y_continuous(breaks = (-2:2) * 30)+
      facet_grid(.~Genus)
    
    plotName<-paste("core_CESM_St2_m3_max_runoff_560to280.png")
    ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
           scale = 1, height = 22, width = 17, units = c("in"),
           dpi = 300, limitsize = FALSE)
    
  }

  #interaction between genus and co2
  if(TRUE){
  water_intGenusCo2 <- core_means %>% filter(Genus %in% c("lepidodendron","treefern")) %>%
    tidyr::pivot_wider(id_cols = c("Lon","Lat","Genus"),
                       names_from = c("co2_f"), 
                       values_from = c(soilw_mean, soilw_outflow_max, soilw_outflow_mean,Area
                       ))%>%
    tidyr::pivot_wider(id_cols = c("Lon","Lat"),
                       names_from = c("Genus"), 
                       values_from = c(soilw_outflow_max_280, soilw_outflow_mean_280,Area_560,
                                       soilw_outflow_max_560, soilw_outflow_mean_560
                       ))%>%
    dplyr::mutate(
      m3_total_runoff_Lep560toTf280 = 
        (soilw_outflow_mean_280_treefern-soilw_outflow_mean_560_lepidodendron)*365/1000*Area_560_treefern,
      m3_max_runoff_Lep560toTf280 = 
        (soilw_outflow_max_280_treefern-soilw_outflow_max_560_lepidodendron)/1000*Area_560_treefern,
      m3_total_runoff_Tf560toLep280 = 
        (soilw_outflow_mean_280_lepidodendron-soilw_outflow_mean_560_treefern)*365/1000*Area_560_treefern,
      m3_max_runoff_Tf560toLep280 = 
        (soilw_outflow_max_280_lepidodendron-soilw_outflow_max_560_treefern)/1000*Area_560_treefern)%>%
    pivot_longer(cols = c("m3_total_runoff_Lep560toTf280","m3_max_runoff_Lep560toTf280",
                         "m3_total_runoff_Tf560toLep280","m3_max_runoff_Tf560toLep280"),
                 names_to = "contrast"
                )
  
  
  
  for(i in unique(water_intGenusCo2$contrast)){
    temp <- water_intGenusCo2%>%filter(contrast == i)
     
    if(TRUE){
      colors<-c("#DE4968FF","#DE4968FF", "#FE9F6DFF", "#FCFDBFFF","white","lightsteelblue1","lightsteelblue3", "#0D0887FF")
      
      this_var <- temp$value
      this_sd<-sd(this_var, na.rm=TRUE)
      this_mean <- mean(this_var, na.rm=TRUE)
      this_max <- max(this_var, na.rm=T)
      this_min <- min(this_var, na.rm = T)
      val <- scales::rescale(x = c(this_min,this_mean-2*this_sd, this_mean-this_sd, 
                                   this_mean,
                                   0, 
                                   this_mean+this_sd,this_mean+2*this_sd,this_max),
                             to = c(0,1), 
                             from = c(this_min,this_max))
    }
    
  pA <- temp%>%ggplot()+
    geom_tile(aes(x=Lon, y=Lat, fill=value))+ 
    geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
    scale_fill_gradientn(colors=colors , values = val) +theme_minimal()+ ##options A-H #values = v, breaks=l,
    theme(text = element_text(size=60), 
          legend.position = 'bottom', legend.key.width = unit(2,'in'),
          axis.text.x=element_blank())+#, legend.key.width = unit(1,'in')
    guides(fill=guide_colorbar(title=""))+
    labs(x="")+
    coord_map(projection="mollweide",ylim = c(-90,75), xlim = c(100,295))+
    scale_y_continuous(breaks = (-2:2) * 30)+
    facet_grid(.~contrast)
  
  plotName<-paste("core_CESM_St2_",i,".png")
  ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
         scale = 1, height = 22, width = 17, units = c("in"),
         dpi = 300, limitsize = FALSE)
  }
  
  
}

} 

#change in LAI due to CO2
if(TRUE){
  lai_change_Co2 <- core_means %>%
    tidyr::pivot_wider(id_cols = c("Lon","Lat","Group"),
                       names_from = c("co2_f"), 
                       values_from = c(epv.proj_lai_mean))%>%
    dplyr::mutate(
      lai_560to280 = `280`-`560`      )
    
  
 
  if(TRUE){
    delLai_colors<-c("brown","#FE9F6DFF","goldenrod2","white","white","lightgreen","green", "darkgreen")

    delLai_colors<-c("brown","brown","white","darkgreen", "darkgreen")
    this_var <- lai_change_Co2$lai_560to280
    this_sd<-sd(this_var, na.rm=TRUE)
    this_mean <- mean(this_var, na.rm=TRUE)
    this_max <- max(this_var, na.rm=T)
    this_min <- min(this_var, na.rm = T)
    val <- scales::rescale(x = c(this_min,this_mean-2*this_sd, this_mean-this_sd, 
                                 0, 
                                 this_mean,
                                 this_mean+this_sd,this_mean+2*this_sd,this_max),
                           to = c(0,1), 
                           from = c(this_min,this_max))
    val <- scales::rescale(x = c(this_min,this_mean-4*this_sd,  
                                 0, 
                                 this_mean+4*this_sd,this_max),
                           to = c(0,1), 
                           from = c(this_min,this_max))
    }
  
  pA <- lai_change_Co2%>%ggplot()+
    geom_tile(aes(x=Lon, y=Lat, fill=lai_560to280))+ 
    geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
    scale_fill_gradientn(colors=delLai_colors , values = val) +theme_minimal()+ ##options A-H #values = v, breaks=l,
    theme(text = element_text(size=60), 
          legend.position = 'bottom', legend.key.width = unit(2,'in'),
          axis.text.x=element_blank(),
          strip.text.y =element_text(size=30))+#, legend.key.width = unit(1,'in')
    guides(fill=guide_colorbar(title=""))+
    labs(x="")+
    coord_map(projection="mollweide",ylim = c(-90,75), xlim = c(100,295))+
    scale_y_continuous(breaks = (-2:2) * 30)+
    facet_grid(Group~.)
  
  plotName<-paste("delLai_byGroup_560to280.png")
  ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
         scale = 1, height = 22, width = 17, units = c("in"),
         dpi = 300, limitsize = FALSE)
  
  
  
  for(i in unique(lai_change_Co2$Genus)){
    temp <- lai_change_Co2%>%filter(Genus == i)
    
  
    
    pA <- temp%>%ggplot()+
      geom_tile(aes(x=Lon, y=Lat, fill=lai_560to280))+ 
      geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
      scale_fill_gradientn(colors=colors , values = val) +theme_minimal()+ ##options A-H #values = v, breaks=l,
      theme(text = element_text(size=60), 
            legend.position = 'bottom', legend.key.width = unit(2,'in'),
            axis.text.x=element_blank())+#, legend.key.width = unit(1,'in')
      guides(fill=guide_colorbar(title=""))+
      labs(x="")+
      coord_map(projection="mollweide",ylim = c(-90,75), xlim = c(100,295))+
      scale_y_continuous(breaks = (-2:2) * 30)+
      facet_grid(.~contrast)
    
    plotName<-paste("core_CESM_St2_",i,".png")
    ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
           scale = 1, height = 22, width = 17, units = c("in"),
           dpi = 300, limitsize = FALSE)
  }
  
  
}

#precip stats
if(TRUE){
   precipstats <- core_means %>% filter(Genus == "treefern") %>%
     group_by(co2_f,Lat,Lon)%>%
     mutate(prcp_cv = prcp_sd/prcp_mean)%>%
    select(prcp_mean, prcp_max, prcp_cv, Area)%>%
     # mutate(prcp_mean=prcp_mean*Area, 
     #        prcp_max=prcp_max*Area)%>%
     mutate(prcp_mean=prcp_mean*0.1*365,
            prcp_max=prcp_max*0.1)%>%
   tidyr::pivot_longer(names_to = "stat", names_prefix = "prcp_",
                       cols=c("prcp_mean",
                              "prcp_max",
                              "prcp_cv"))%>%
     mutate( 
            stat = factor(stat, levels=c("mean", "cv","max")))%>%
   tidyr::pivot_wider(id_cols = c("Lon","Lat","stat"),
                      names_from = c("co2_f"), 
                      values_from = value) %>%
     dplyr::mutate(
       delCO2560to280 = `280` - `560` ,
       delCO2280to560 = `560` - `280` )
     

   
   prcp_mean <- precipstats %>% filter(stat=="mean")
   prcp_max <- precipstats %>% filter(stat=="max")
   prcp_cv <- precipstats %>% filter(stat=="cv")
   
  
   # colors_topext<-rev(plasma(6))
   # colors_topext<-c(colors_topext, 
   #                  colors_topext[6],colors_topext[6],colors_topext[6],colors_topext[6],
   #                  colors_topext[6],colors_topext[6],colors_topext[6],colors_topext[6])
 ##############################################differences in precip plots  
   if(TRUE){
     colors<-c("red","red","white","blue","blue")
     # colors_ext<-c(colors[1],colors,colors[6])
     
     this_var <-prcp_mean$delCO2280to560
     this_sd<-sd(this_var, na.rm=TRUE)
     this_mean <- mean(this_var, na.rm=TRUE)
     this_max <- max(this_var, na.rm=T)
     this_min <- min(this_var, na.rm = T)
     val <- scales::rescale(x = c(this_min,this_mean-3*this_sd, 
                                  0, 
                                  this_mean+3*this_sd,this_max),
                            to = c(0,1), 
                            from = c(this_min,this_max))
   }
   
   pA <- ggplot()+
     geom_tile(data = prcp_mean, aes(x=Lon, y=Lat, fill=delCO2280to560))+ 
     scale_fill_gradientn(colors=colors,values = val)+#,breaks=seq(from=50,to=600, by=50))+
     geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
     theme_minimal()+ ##options A-H #values = v, breaks=l,
     theme(text = element_text(size=60), 
           # legend.position = 'bottom', legend.key.width = unit(2,'in'),
           legend.key.height = unit(0.75,'in'),
           axis.text.x=element_blank())+#, legend.key.width = unit(1,'in')
     guides(fill=guide_colorbar(title=""))+
     labs(x="")+
     coord_map(projection="mollweide",ylim = c(-90,75), xlim = c(100,295))+
     scale_y_continuous(breaks = (-2:2) * 30)
   
   
   plotName<-paste("dprecip_mean_cm.png")
   ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
          scale = 1, height = 10, width = 17, units = c("in"),
          dpi = 300, limitsize = FALSE)
  
   if(TRUE){
     colors<-c("red","red","white","blue","blue")
     # colors_ext<-c(colors[1],colors,colors[6])
     
     this_var <-prcp_cv$delCO2280to560
     this_sd<-sd(this_var, na.rm=TRUE)
     this_mean <- mean(this_var, na.rm=TRUE)
     this_max <- max(this_var, na.rm=T)
     this_min <- min(this_var, na.rm = T)
     val <- scales::rescale(x = c(this_min,this_mean-3*this_sd, 
                                  0, 
                                  this_mean+3*this_sd,this_max),
                            to = c(0,1), 
                            from = c(this_min,this_max))
   }
   
   pA <- ggplot()+  geom_tile(data= prcp_cv, aes(x=Lon, y=Lat, fill=delCO2280to560))+ 
     scale_fill_gradientn(colors=colors, values = val)+
     geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
     theme_minimal()+ ##options A-H #values = v, breaks=l,
     theme(text = element_text(size=60), 
           # legend.position = 'bottom', legend.key.width = unit(2,'in'),
           legend.key.height = unit(0.75,'in'),
           axis.text.x=element_blank())+#, legend.key.width = unit(1,'in')
     guides(fill=guide_colorbar(title=""))+
     labs(x="")+
     coord_map(projection="mollweide",ylim = c(-90,75), xlim = c(100,295))+
     scale_y_continuous(breaks = (-2:2) * 30)
   
   
   plotName<-paste("dprecip_cv_cm.png")
   ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
          scale = 1, height = 10, width = 17, units = c("in"),
          dpi = 300, limitsize = FALSE)
     
   if(TRUE){
     colors<-c("red","red","white","blue","blue")
     # colors_ext<-c(colors[1],colors,colors[6])
     
     this_var <-prcp_max$delCO2280to560
     this_sd<-sd(this_var, na.rm=TRUE)
     this_mean <- mean(this_var, na.rm=TRUE)
     this_max <- max(this_var, na.rm=T)
     this_min <- min(this_var, na.rm = T)
     val <- scales::rescale(x = c(this_min,this_mean-3*this_sd, 
                                  0, 
                                  this_mean+3*this_sd,this_max),
                            to = c(0,1), 
                            from = c(this_min,this_max))
   }
   
   pA <- ggplot()+ geom_tile(data= prcp_max, aes(x=Lon, y=Lat, fill=delCO2280to560))+ 
     scale_fill_gradientn(colors=colors, values = val)+
     geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
     theme_minimal()+ ##options A-H #values = v, breaks=l,
     theme(text = element_text(size=60), 
           # legend.position = 'bottom', legend.key.width = unit(2,'in'),
           legend.key.height = unit(0.75,'in'),
           axis.text.x=element_blank())+#, legend.key.width = unit(1,'in')
     guides(fill=guide_colorbar(title=""))+
     labs(x="")+
     coord_map(projection="mollweide",ylim = c(-90,75), xlim = c(100,295))+
     scale_y_continuous(breaks = (-2:2) * 30)
   
   
   plotName<-paste("dprecip_max_cm.png")
   ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
          scale = 1, height = 10, width = 17, units = c("in"),
          dpi = 300, limitsize = FALSE)
     
     
   
  
   # prcp_m3 = value*(1/100)*Area
   
   #one plot same scale, problem
   if(FALSE){ 
     pA <- precipstats%>%ggplot()+
     geom_tile(aes(x=Lon, y=Lat, fill=prcp_m3))+ 
     geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
     scale_fill_gradientn(colors=rev(plasma(6)))+theme_minimal()+ ##options A-H #values = v, breaks=l,
     theme(text = element_text(size=30), 
           legend.position = 'bottom', legend.key.width = unit(2,'in'),
           axis.text.x=element_blank())+#, legend.key.width = unit(1,'in')
     guides(fill=guide_colorbar(title=""))+
     labs(x="")+
     coord_map(projection="mollweide",ylim = c(-90,75), xlim = c(60,345))+
     scale_y_continuous(breaks = (-2:2) * 30)+
     facet_grid(stat~co2_f, scales = "free_y")
   
   plotName<-paste("precip_stats.png")
   ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
          scale = 1, height = 22, width = 17, units = c("in"),
          dpi = 300, limitsize = FALSE)}
   
   #fix using ggnewscale
   if(TRUE){
     prcp_mean <- precipstats %>% filter(stat=="mean")
     prcp_max <- precipstats %>% filter(stat=="max")
     prcp_cv <- precipstats %>% filter(stat=="cv")
     
     library(ggnewscale)
     colors_topext<-rev(plasma(6))
     colors_topext<-c(colors_topext, 
                      colors_topext[6],colors_topext[6],colors_topext[6],colors_topext[6],
                      colors_topext[6],colors_topext[6],colors_topext[6],colors_topext[6])
     
     pA <- precipstats%>%filter(stat!="min")%>%ggplot()+
       geom_tile(data= prcp_mean, aes(x=Lon, y=Lat, fill=value))+ 
       scale_fill_gradientn(colors=colors_topext,breaks=seq(from=50,to=600, by=50))+
       new_scale_fill()+
       geom_tile(data= prcp_cv, aes(x=Lon, y=Lat, fill=value))+ 
       scale_fill_gradientn(colors=rev(rocket(6)))+
       new_scale_fill()+
       geom_tile(data= prcp_max, aes(x=Lon, y=Lat, fill=value))+ 
       scale_fill_gradientn(colors=rev(plasma(6)))+
       geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
       theme_minimal()+ ##options A-H #values = v, breaks=l,
       theme(text = element_text(size=30), 
             # legend.position = 'bottom', legend.key.width = unit(2,'in'),
             legend.key.height = unit(0.75,'in'),
             axis.text.x=element_blank())+#, legend.key.width = unit(1,'in')
       guides(fill=guide_colorbar(title=""))+
       labs(x="")+
       coord_map(projection="mollweide",ylim = c(-90,75), xlim = c(60,345))+
       scale_y_continuous(breaks = (-2:2) * 30)+
       facet_grid(stat~co2_f)
     
     plotName<-paste("precip_stats_cm.png")
     ggsave(plotName, plot = pA, device = "png", path = "/EcosystemEngineers2024/",
            scale = 1, height = 22, width = 17, units = c("in"),
            dpi = 300, limitsize = FALSE)
  
   
   }
   
}    

#low vs met
if(TRUE){
  # annual_means[grep(pattern="_min",names(annual_means))]
  install.packages("ggthemes")
  library(ggthemes)
  lowLAI <- annual_means%>%filter(epv.proj_lai_median<1)
  epcs<-read_excel(path = "EE24_WJM_Table1.xlsx")
  pA <- suffDay_means%>%ggplot()+
    geom_boxplot(aes(y=soil_psi_min, x = Group))+
    geom_tile(data=epcs, aes(x=Group, y=ΨC), fill="red")+
    geom_hline(aes(yintercept =-2),color="green")+
    coord_cartesian(ylim = c(-9,0.1))+
    theme_classic()+
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 60, vjust = .9, hjust=.9))+
    labs(x="Plant Group",y="Minimum Soil Water Potential (Ψs ;MPa)")+
    facet_grid(.~co2_f)
   
  plotName<-paste("suffDayl_soilPsiMin.png")
  ggsave(plotName, plot = pA, device = "png", path = "/Users/willmatthaeus/Dropbox/bgc_work/stem_bgc/EcosystemEngineers2024/",
         scale = 1, height = 8.5, width = 11, units = c("in"),
         dpi = 300, limitsize = FALSE)
  
  
  pA <- suffDay_means%>%ggplot()+
    geom_point(aes(x=soil_psi_median, y = gl_s_sun_median+gl_s_shade_median))+
    # coord_cartesian(ylim = c(-9,0.1))+
    theme_classic()+
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 60, vjust = .9, hjust=.9),
          strip.text.x = element_text(size = 10))+
    labs(y="Stomatal Conductance",x="Median Soil Water Potential")+
    facet_grid(co2_f~Group,)
  
  plotName<-paste("suffDayl_soilpsimedVSglsmed.png")
  ggsave(plotName, plot = pA, device = "png", path = "/Users/willmatthaeus/Dropbox/bgc_work/stem_bgc/EcosystemEngineers2024/",
         scale = 1, height = 8.5, width = 11, units = c("in"),
         dpi = 300, limitsize = FALSE)
  
  enough_day <-annual_means%>%filter(epv.proj_lai_median>0)%>%group_by(co2_f)%>%summarize(minDayl = min(dayl_median))
  
  suffDay_means <- annual_means %>% group_by(co2_f) %>% 
            left_join(y = enough_day,by = "co2_f")%>%
            filter(dayl_median>minDayl)
  
  
}

##survivors
if(TRUE){
  '%!in%' <- function(x,y)!('%in%'(x,y))

  
  survivor_summary<-annual_means %>% group_by(Tracheid_type, co2xmodern) %>% filter(epv.proj_lai>0.1) %>%
    dplyr::summarise(survivors = list(unique(Genus)), n_survivors = length(unlist(survivors)))
  
  survivor_spatial_summary <- annual_means %>% group_by(Lon, Lat, co2xmodern) %>% filter(epv.proj_lai>0.1) %>%
    dplyr::summarise(survivors = list(unique(Genus)), n_survivors = length(unlist(survivors))) 
  
  survivor_spatial_summary_tt <- annual_means %>% group_by(Lon, Lat, Tracheid_type, co2xmodern) %>% filter(epv.proj_lai>0.1) %>%
    dplyr::summarise(survivors = list(unique(Genus)), n_survivors = length(unlist(survivors))) 
  
   survivor_spatial_summary <- survivor_spatial_summary %>% mutate(survivor_f = as.factor(paste(unlist(survivors), collapse = ", ")))
  survivor_spatial_summary_tt <- survivor_spatial_summary_tt %>% mutate(survivor_f = as.factor(paste(unlist(survivors), collapse = ", ")))
  
  
  survivor_spatial_summary %>% group_by(survivor_f, co2xmodern) %>% 
    summarize(freq = n()) %>% write_tsv(file = "CESM_St2_surivor_assemblasges.tsv")
  
  survivor_spatial_summary_tt %>% group_by(survivor_f, co2xmodern) %>% 
    summarize(freq = n()) %>% write_tsv(file = "CESM_St2_surivor_tt_assemblasges.tsv")
  
  
  mask <- annual_means %>% 
    dplyr::select(Lon, Lat, Model, Fc, Genus) %>% 
    filter(Model == "Eud", Fc=="05",Genus=="arthropitys") %>% 
    dplyr::select(Lon, Lat)
  
  mask$vegland <- 1

  m<-mask%>%as_spatraster(xycols=1:2, crs = "" ,digits=1)
  m_fort<-m%>%raster()%>%rasterToPolygons(dissolve=TRUE,n=16)%>%fortify()

  
  
  pA <- ggplot(data=survivor_spatial_summary)+
    geom_point(aes(x=Lon, y=Lat, color=n_survivors),shape=15, size=2.5)+
    geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
    viridis::scale_color_viridis()+theme_minimal()+
    theme(text = element_text(size=30))+
    facet_grid(co2xmodern~.)
  
  plotName<-"survivors.png"
  # plotName<-"Stem2_paleo_nppSubDom_Cover_m20.png"
  ggsave(plotName, plot = pA, device = "png", path = location,
         scale = 1, height = 25, width = 15, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  ####
  pA <- ggplot(data=survivor_spatial_summary)+
    geom_point(aes(x=Lon, y=Lat, color=n_survivors),shape=15, size=2.5)+
    geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
    viridis::scale_color_viridis()+theme_minimal()+
    theme(text = element_text(size=30))+
    facet_wrap(.~co2xmodern,ncol = 1)
    
  plotName<-"survivors.png"
  # plotName<-"Stem2_paleo_nppSubDom_Cover_m20.png"
  ggsave(plotName, plot = pA, device = "png", path = location,
         scale = 1, height = 25, width = 15, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  ##try same idea but with discrete community values
  # unique(survivor_spatial_summary$survivor_f)
  thisCoverCbPalette <- c("darkgreen","#E69F00", "#009E73","#0072B2","black","grey", "#56B4E9")
  survivor_spatial_summary$survivor_f <- factor(survivor_spatial_summary$survivor_f, 
            levels = c("arthropitys, lepidodendron, spenophyllum, treefern",
                       "arthropitys, treefern",
                       "arthropitys, spenophyllum, treefern",
                       "arthropitys, spenophyllum",
                       "arthropitys",
                       "treefern",
                       "arthropitys, lepidodendron, treefern" ))
  names(thisCoverCbPalette)<-levels(survivor_spatial_summary$survivor_f)
  summary(survivor_spatial_summary$Lon)
  pA <- ggplot(data=survivor_spatial_summary)+
    geom_tile(aes(x=Lon, y=Lat, fill=survivor_f))+ #size=.25*((75-abs(Lat))/75)),,shape=15
    geom_path(data=m_fort, aes(x = long, y = lat, group = group),linewidth =1, col="black")+
    scale_fill_manual(values = thisCoverCbPalette)+theme_minimal()+
    theme(text = element_text(size=30), legend.position = 'none')+
    coord_map(projection="mollweide",ylim = c(-90,75), xlim = c(60,345))+
    facet_wrap(Fc~Model,ncol = 1)
  
  plotName<-"FC_Survivors_Moll.png"
  ggsave(plotName, plot = pA, device = "png", path = "/Users/willmatthaeus/Dropbox/bgc_work/stem_bgc/FC_Displacement/",
         scale = 1, height = 25, width = 15, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  
  prcpCbPalette <- c("#F0E442","#56B4E9")
  pA <- ggplot()+
    geom_point(data=thisCover, aes(x=LON, y=LAT, fill=Mean_prcp*365/10),shape=22, size=4, stroke = 0)+
    scale_fill_gradient(low = "#F0E442", high = "#000099",breaks = seq(50,350,50))+theme_bw()+
    theme(text = element_text(size=30))+
    geom_text(data=example_locations,aes(x=LON, y=LAT,label=Name),hjust=1.2, vjust=-1., size=8)+
    geom_point(data=example_locations, aes(x=LON, y=LAT), color="black", size=11, shape=0, stroke=2)+
    guides(fill=guide_legend(title="Mean Ann. Prcp"))
  
  
  plotName<-"Stem2_precip_ex.png"
  ggsave(plotName, plot = pA, device = "png", path = "/Users/willmatthaeus/Dropbox/bgc_work/stem2/",
         scale = 1, height = 10, width = 15, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  
  pA <- ggplot()+
    geom_point(data=mean_mets_G, aes(x=LON, y=LAT, fill=survivor_number),shape=22, size=4, stroke = 0)+
    scale_fill_viridis()+theme_bw()+
    theme(text = element_text(size=30))+
    geom_text(data=example_locations,aes(x=LON, y=LAT,label=Name),hjust=1.2, vjust=-1., size=8)+
    geom_point(data=example_locations, aes(x=LON, y=LAT), color="black", size=11, shape=0, stroke=2)+
    guides(fill=guide_legend(title="Surviving PFTs\n(LAI >0.1)"))
  
  
  plotName<-"survivors_Glacial.png"
  ggsave(plotName, plot = pA, device = "png", path = "/Users/willmatthaeus/Dropbox/bgc_work/stem2/",
         scale = 1, height = 10, width = 15, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  pA <- ggplot()+
    geom_point(data=mean_mets_IG, aes(x=LON, y=LAT, fill=survivor_number),shape=22, size=4, stroke = 0)+
    scale_fill_viridis()+theme_bw()+
    theme(text = element_text(size=30))+
    geom_text(data=example_locations,aes(x=LON, y=LAT,label=Name),hjust=1.2, vjust=-1., size=8)+
    geom_point(data=example_locations, aes(x=LON, y=LAT), color="black", size=11, shape=0, stroke=2)+
    guides(fill=guide_legend(title="Surviving PFTs\n(LAI >0.1)"))
  
  
  plotName<-"survivors_InterGlacial.png"
  ggsave(plotName, plot = pA, device = "png", path = "/Users/willmatthaeus/Dropbox/bgc_work/stem2/",
         scale = 1, height = 10, width = 15, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  range(diffs$diff_LAI, na.rm = TRUE)
  pA <- ggplot()+
    geom_point(data=diffs, aes(x=LON, y=LAT, fill=diff_LAI),shape=22, size=4, stroke = 0)+
    scale_fill_gradient()+theme_bw()+
    theme(text = element_text(size=30))+
    geom_text(data=example_locations,aes(x=LON, y=LAT,label=Name),hjust=1.2, vjust=-1., size=8)+
    geom_point(data=example_locations, aes(x=LON, y=LAT), color="black", size=11, shape=0, stroke=2)+
    guides(fill=guide_legend(title="LAI differences"))+facet_grid(Species~CO2)
  
  
  plotName<-"model_diff_LAI_rect.png"
  ggsave(plotName, plot = pA, device = "png", path = "/Users/willmatthaeus/Dropbox/bgc_work/stem2/",
         scale = 1, height = 20, width = 30, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
}

#community
if(TRUE){
  library(data.table)
  annual_mean_dt <- as.data.table(annual_means)
  annual_mean_dt <- annual_mean_dt %>% mutate(loc = paste(Lat,Lon))
  dt_1x <- annual_mean_dt %>% filter(co2xmodern=="1x")
  dt_1x[dt_1x[, .I[epv.proj_lai == max(epv.proj_lai)], by=loc]$V1]
}

##Whitaker charts with CESM colored by survival
if(TRUE){
  
  clr <- viridis(6)
  clr <- c(clr, "#FDE725FF","#FDE725FF","#FDE725FF")
  mean_mets <- core_means %>% filter(epv.proj_lai_median>0.1) %>% 
    mutate(loc = paste(Lat,Lon))%>%
    group_by(loc,co2_f,Genus)%>%
    summarize(precp_cm=mean(prcp_mean)*365/10, temp_c=mean(tavg_mean), lai = mean(epv.proj_lai_median))
  
  for(g in unique(mean_mets$Genus)){
    temp<-mean_mets%>%filter(Genus==g)
    
    pA<-ggplot()+
      geom_polygon(data = Whittaker_biomes,
                   aes(x    = temp_c,
                       y    = precp_cm,
                       fill = biome),
                   # adjust polygon borders
                   colour = "gray98",
                   size   = 1) +
      geom_point(data=temp, aes(x=temp_c, y=precp_cm), size=0.5)+ 
      scale_color_gradientn(colors = clr)+facet_grid(co2_f~.)
    
    plotName<-paste0(g,"_Whitaker_InterGlacialCESM_survivors.png")
    # plotName<-"Stem2_paleo_nppSubDom_Cover_m20.png"
    ggsave(plotName, plot = pA, device = "png", path = "/Users/willmatthaeus/Dropbox/bgc_work/stem_bgc/EcosystemEngineers2024/",
           scale = 1, height = 7.5, width = 8, units = c("in"),
           dpi = 300, limitsize = TRUE)
  }
  
  
  #survivors
   mean_mets <- annual_means %>% group_by(Lon, Lat, co2xmodern) %>% filter(epv.proj_lai>0.1) %>%
     mutate(loc = paste(Lat,Lon))%>%
     dplyr::summarize(precp_cm=mean(prcp)*365/10, temp_c=mean(tavg), 
                      lai = mean(epv.proj_lai), vegc = mean(vegc),
                      survivors = list(unique(Genus)), n_survivors = length(unlist(survivors)))
    

   ###
   pA<-ggplot()+
     geom_polygon(data = Whittaker_biomes,
                  aes(x    = temp_c,
                      y    = precp_cm,
                      fill = biome),
                  # adjust polygon borders
                  colour = "gray98",
                  size   = 1) +
     geom_point(data=mean_mets, aes(x=temp_c, y=precp_cm, color = n_survivors), size=0.5)+ 
     scale_color_gradientn(colors = clr)+facet_grid(co2xmodern~.)
   
   plotName<-"Whitaker_InterGlacialCESM_survivors.png"
   # plotName<-"Stem2_paleo_nppSubDom_Cover_m20.png"
   ggsave(plotName, plot = pA, device = "png", path = "/Users/willmatthaeus/Dropbox/bgc_work/stem_bgc/EcosystemEngineers2024/",
          scale = 1, height = 7.5, width = 8, units = c("in"),
          dpi = 300, limitsize = TRUE)
   
}


;load "/Users/wjmatthaeus/miniconda3/pkgs/ncl-6.6.2-h7fcb0a1_21/lib/ncarg/nclex/gsun/gsn_code.ncl"
;load "/Users/wjmatthaeus/miniconda3/pkgs/ncl-6.6.2-h7fcb0a1_21/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "/Users/wjmatthaeus/miniconda3/pkgs/ncl-6.6.2-h7fcb0a1_21/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "/Users/wjmatthaeus/miniconda3/pkgs/ncl-6.6.2-h7fcb0a1_21/lib/ncarg/nclscripts/csm/shea_util.ncl"

;load "/home/matthaeusw/bin/ncl_utils/ncarg/nclex/gsun/gsn_code.ncl"
;load "/home/matthaeusw/bin/ncl_utils/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "/home/matthaeusw/bin/ncl_utils/ncarg/nclscripts/csm/contributed.ncl"
;load "/home/matthaeusw/bin/ncl_utils/ncarg/nclscripts/csm/shea_util.ncl"

;*******************************************************
;WJ Matthaeus 04-2022
;
;ncl code to take cesm climate data and a land surface mask
;and output a set of single-location mtc43 files for
;input to x-BGC 
;*******************************************************

;*******************************************************
;if you're using anaconda remember to activate environment
;
;conda activate ncl_stable
;*******************************************************
begin
;*******************************************************
; open two files and load in monthly climatologies
;*******************************************************

;IGLAC_546_28_99_-83_walchia_MEDCN_MEDg_Vr
;
; = addfile("/data/matthaeusw/CESM2BGC_270121_source/PV_LPIA5ig_B1850C5CN_f19_g16_1xC_WISO.006b.clm2.cam.h1.1502-1551.dailytseries.nc","r")
;cesm_in = addfile("/data/matthaeusw/CESM2BGC_270121_source/PV_LPIA5ig_B1850C5CN_f19_g16_2xC_WISO.006b.clm2.cam.h1.1102-1151.dailytseries.nc","r")
cesm_in = addfile("/data/matthaeusw/CESM2BGC_270121_source/PV_LPIA5ig_B1850C5CN_f19_g16_1xC_WISO.006b.clm2.cam.h1.1502-1551.dailytseries.nc","r")
lat = cesm_in->lat
lon = cesm_in->lon
time = cesm_in->time
nlat = dimsizes(lat)
nlon = dimsizes(lon)
ntim = dimsizes(time)

PRECT = cesm_in->PRECT
FSDS=cesm_in->FSDS

TSA=cesm_in->TSA
TREFMXAV=cesm_in->TREFMXAV
TREFMNAV=cesm_in->TREFMNAV
Q2M=cesm_in->Q2M
DAYLf=cesm_in->DAYL
DAYL = floattoint(DAYLf)

; 96 lats
; 144 lons
; 18250 times

;convert precip from m/s to cm/d
; m/s * 100 cm/m * 86400 s/d
PRECT_CM = PRECT*86400*100
;convert temps to C
TSA_C = TSA-273.5
TREFMXAV_C = TREFMXAV-273.5
TREFMNAV_C = TREFMNAV-273.5
;calcualte VPD (based on svp, Arden-Buck Eq for vapor over water)

svp = 611.21*exp(((TSA_C*(18.678 - TSA_C/234.5))/(257.14 + TSA_C)))
vp = 611.21*exp(((TREFMNAV_C*(18.678 - TREFMNAV_C/234.5))/(257.14 + TREFMNAV_C)))
VPD = svp - vp
DAYOFYEAR = mod(ispan(1, ntim, 1),366)


;missing daylen!!!!!!!!!!!!!!!!!!

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;tests and notes
;these look reasonable 
;res=True
;
;res@gsnHistogramBinIntervals = (/0,250,500,750,1000,1250,1500,1750,2000/)
;wks = gsn_open_wks("png","histo_vpd")
;plot = gsn_histogram(wks,ndtooned(VPD),res)
;
;res@gsnHistogramBinIntervals := (/0,2.5,5.0,7.5,10.0,12.5,15.0,17.5,20.0/)
;wks = gsn_open_wks("png","histo_prect")
;plot = gsn_histogram(wks,ndtooned(PRECT),res)
;
;res@gsnHistogramBinIntervals := (/0,2.5,5.0,7.5,10.0,12.5,15.0,17.5,20.0/)
;wks = gsn_open_wks("png","histo_prect")
;plot = gsn_histogram(wks,ndtooned(PRECT),res)
;
;res@gsnHistogramBinIntervals := (/0,2.5,5.0,7.5,10.0,12.5,15.0,17.5,20.0/)
;wks = gsn_open_wks("png","histo_sw_cm")
;plot = gsn_histogram(wks,ndtooned(sw),res)
;
;res@gsnHistogramBinIntervals := (/-50,-25,0,25,50/)
;wks = gsn_open_wks("png","histo_tsa_c")
;plot = gsn_histogram(wks,ndtooned(TSA_K),res)
;
;res@gsnHistogramBinIntervals := (/-50,-25,0,25,50/)
;wks = gsn_open_wks("png","histo_TREFMXAV_c")
;plot = gsn_histogram(wks,ndtooned(TREFMXAV),res)
;
;res@gsnHistogramBinIntervals := (/-50,-25,0,25,50/)
;wks = gsn_open_wks("png","histo_TREFMNAV_c")
;plot = gsn_histogram(wks,ndtooned(TREFMNAV),res)
;
;aveP := avg(PRECT)
;aveF := avg(FSDS)
;aveV := avg(VPD)
;aveP := dim_avg_n(PRECT,(/0,1/))
;printVarSummary(aveP)

;new((/ntim, nlon, nlat/), float)
; ws_soilw@name = "ws_soilw"
; ws_soilw!0 = "time"
; ws_soilw!1 = "lon"
; ws_soilw!2 = "lat"
; filevardef(f, "ws_soilw", typeof(ws_soilw), getvardims(ws_soilw))
; filevarattdef(f,"ws_soilw", ws_soilw)

;calculate vpd

;columns in .dayout.ascii
;ncol = 60
;
;;loop over all dayout.ascii that uses glac_co2 , lon, lat, species to assign subsets of variables
;;not all locations are represented
;;track number of present files (not ocean,ice in grid)
;not_missing = 0

;err = NhlGetErrorObjectId()
;setvalues err
;  "errPrint" : "False" ; turn off all error msg printing
;end setvalues

;setvalues err
;  "errPrint" : "True" ; turn off all error msg printing
;end setvalues


;printVarSummary(DAYOFYEAR)
;printVarSummary(TREFMXAV_C(:,i,j))
;printVarSummary(TREFMNAV_C(:,i,j))
;printVarSummary(TSA_C(:,i,j))
;printVarSummary(PRECT_CM(:,i,j))
;printVarSummary(VPD(:,i,j))
;printVarSummary(FSDS(:,i,j))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;tests
fil = addfile("/data/matthaeusw/ncl_work/CESM/surfdata_ocnlndice_LPIA5ig_200611.nc","r") 
;lmask = dim_avg_n_Wrap(fil->mask,0)
lmask = fil->mask
lat = fil->lat
lon = fil->lon

;i=15
;j=15
;out_list = [/TREFMXAV_C,TREFMNAV_C,TSA_C,PRECT_CM,VPD,FSDS/]
;system("rm /data/matthaeusw/CESM2BGC_270121/allMets_CESM2BGC_270121/*.mtc43")
;;system("mkdir /data/matthaeusw/CESM2BGC_270121/allMets_CESM2BGC_270121")
;system("cd /data/matthaeusw/CESM2BGC_270121/allMets_CESM2BGC_270121")

do i = 0, nlat-1
;   this_lon = floattoint(lon(i))
       this_lat = lat(i)
  do j = 0, nlon-1
        this_lon = lon(j)
        ;test for LMASK here
        if (lmask({this_lat},{this_lon}).eq.1) then;landmask.. check for ice-free land locations
;            print("here")
             out_list = [/DAYOFYEAR,TREFMXAV_C(:,i,j),TREFMNAV_C(:,i,j),TSA_C(:,i,j),PRECT_CM(:,i,j),VPD(:,i,j),FSDS(:,i,j),DAYL(:,i,j)/]
            fname = this_lon+"_"+this_lat+".mtc43NH"
;           write_table(fname, "w", out_list, "%d%0.2f%0.2f%0.2f%0.2f%0.2f%0.2f")
            write_table(fname, "w", out_list, "%d %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %d")    
        end if
;outDF2<-outDF2%>%mutate(Yearday=as.integer(Yearday),`Tmax (deg C)`=sprintf("%0.2f",`Tmax (deg C)`),`Tmin (deg C)`=sprintf("%0.2f",`Tmin (deg C)`),
;    `Tday (deg C)`=sprintf("%0.2f",`Tday (deg C)`),`Prcp (cm)`=sprintf("%0.2f",`Prcp (cm)`),VPD=sprintf("%0.2f",VPD),`Srad (W/m^2)`=sprintf("%0.2f",`Srad (W/m^2)`),
;    `Daylen (Sec)`=sprintf("%0.2f",as.numeric(`Daylen (Sec)`)))
 end do
end do

delete(fil)
delete(cesm_in)

;bash commands to add header cesm_to_mtc43_NOTES.sh

;===================================================================
;Generate time and space coordinate variable.
; time = ispan(1, ntim, 1)
; time@name = "time"
; time!0 = "time"
; filevardef(f, "time", typeof(time), getvardims(time))
; filevarattdef(f,"time", time)
; Write the data out.
; f->time  = (/time/)
; 
;
;
;take lon/lat from landmask readin above
;  lon@name = "lon"
; lon!0 = "lon"
; filevardef(f, "lon", typeof(lon), getvardims(lon))
; filevarattdef(f,"lon", lon)
;  f->lon  = (/lon/)
; 
;
;  lat@name = "lat"
; lat!0 = "lat"
; filevardef(f, "lat", typeof(lat), getvardims(lat))
; filevarattdef(f,"lat", lat)
;  f->lat  = (/lat/)
;
; 
;ws_soilw = new((/ntim, nlon, nlat/), float)
; ws_soilw@name = "ws_soilw"
; ws_soilw!0 = "time"
; ws_soilw!1 = "lon"
; ws_soilw!2 = "lat"
; filevardef(f, "ws_soilw", typeof(ws_soilw), getvardims(ws_soilw))
; filevarattdef(f,"ws_soilw", ws_soilw)
;
;epv_proj_lai = new((/ntim, nlon, nlat/), float)
; epv_proj_lai@name = "epv_proj_lai"
; epv_proj_lai!0 = "time"
; epv_proj_lai!1 = "lon"
; epv_proj_lai!2 = "lat"
; filevardef(f, "epv_proj_lai", typeof(epv_proj_lai), getvardims(epv_proj_lai))
; filevarattdef(f,"epv_proj_lai", epv_proj_lai)
;
;m_psi  = new((/ntim, nlon, nlat/), float)
; m_psi@name = "m_psi"
; m_psi!0 = "time"
; m_psi!1 = "lon"
; m_psi!2 = "lat"
; filevardef(f, "m_psi", typeof(m_psi), getvardims(m_psi))
; filevarattdef(f,"m_psi", m_psi)
;
;epv_psi = new((/ntim, nlon, nlat/), float)
; epv_psi@name = "epv_psi"
; epv_psi!0 = "time"
; epv_psi!1 = "lon"
; epv_psi!2 = "lat"
; filevardef(f, "epv_psi", typeof(epv_psi), getvardims(epv_psi))
; filevarattdef(f,"epv_psi", epv_psi)
;
;leaf_psi  = new((/ntim, nlon, nlat/), float)
; leaf_psi@name = "leaf_psi"
; leaf_psi!0 = "time"
; leaf_psi!1 = "lon"
; leaf_psi!2 = "lat"
; filevardef(f, "leaf_psi", typeof(leaf_psi), getvardims(leaf_psi))
; filevarattdef(f,"leaf_psi", leaf_psi)
;
;epv_gl_s_sun  = new((/ntim, nlon, nlat/), float)
; epv_gl_s_sun@name = "epv_gl_s_sun"
; epv_gl_s_sun!0 = "time"
; epv_gl_s_sun!1 = "lon"
; epv_gl_s_sun!2 = "lat"
; filevardef(f, "epv_gl_s_sun", typeof(epv_gl_s_sun), getvardims(epv_gl_s_sun))
; filevarattdef(f,"epv_gl_s_sun", epv_gl_s_sun)
;
;wf_soilw_outflow  = new((/ntim, nlon, nlat/), float)
; wf_soilw_outflow@name = "wf_soilw_outflow"
; wf_soilw_outflow!0 = "time"
; wf_soilw_outflow!1 = "lon"
; wf_soilw_outflow!2 = "lat"
; filevardef(f, "wf_soilw_outflow", typeof(wf_soilw_outflow), getvardims(wf_soilw_outflow))
; filevarattdef(f,"wf_soilw_outflow", wf_soilw_outflow)
;
;nleached_snk = new((/ntim, nlon, nlat/), float)
; nleached_snk@name = "nleached_snk"
; nleached_snk!0 = "time"
; nleached_snk!1 = "lon"
; nleached_snk!2 = "lat"
; filevardef(f, "nleached_snk", typeof(nleached_snk), getvardims(nleached_snk))
; filevarattdef(f,"nleached_snk", nleached_snk)
; 
;
;
;
;;test assignment.. seems good
;delete(f)
;

;NOTES
;fil2 = addfile("/data/matthaeusw/LMASK_182.28_41.nc","r") 
;lmask2 = dim_avg_n_Wrap(fil2->LMASK,0)
;lat = fil2->lat
;lon = fil2->lon
;
;fn = "/data/matthaeusw/test_nc4_multiudim.nc"
;f = addfile(fn,"r")
;vn = getfilevarnames(f)
;t = f->time
;sw = f->ws_soilw
;l = f->lon
;;
;;
;sw_195_27 = sw(:,{this_lon},{this_lat})
;;sw_195_27 = sw(:,:,this_lat)
;
;lmask_ext_time = conform_dims((/18250,90,180/))
;sw_vegland = where(lmask2.eq.1,sw,sw@_FillValue); ice free land only
;
;wks  = gsn_open_wks("png","sw_195_27")                 ; send graphics to PNG file
;res                 = True         
;res@tiMainString    = "27N 15E vs 23.72S 149W"  ; title
;res@tiXAxisString   = "SW"          
;x =  sw_195_27&time
;y = sw_195_27
;plot = gsn_xy(wks,sw_195_27&time,sw_195_27,res)

;f->ws_soilw(:,{this_lon},{this_lat})            = data(:,0);;;;;;;;;
;f->ws_snoww(:,{this_lon},{this_lat})           = data(:,1)
;f->wf_canopyw_evap(:,{this_lon},{this_lat})    = data(:,2)
;f->wf_snoww_subl (:,{this_lon},{this_lat})     = data(:,3)
;f->wf_soilw_evap(:,{this_lon},{this_lat})      = data(:,4)
;f->wf_soilw_trans(:,{this_lon},{this_lat})     = data(:,5)
;f->wf_soilw_outflow(:,{this_lon},{this_lat})   = data(:,6);;;;;;;;
;f->cs_leafc(:,{this_lon},{this_lat})           = data(:,7)
;f->cs_frootc(:,{this_lon},{this_lat})          = data(:,8)
;f->cs_livestemc(:,{this_lon},{this_lat})       = data(:,9)
;f->cs_deadstemc(:,{this_lon},{this_lat})       = data(:,10)
;f->cs_livecrootc(:,{this_lon},{this_lat})      = data(:,11)
;f->cs_deadcrootc(:,{this_lon},{this_lat})      = data(:,12)
;f->epv_proj_lai(:,{this_lon},{this_lat})       = data(:,13);;;;;;;;
;f->epv_psi(:,{this_lon},{this_lat})            = data(:,14);;;;;;;;
;f->epv_daily_net_nmin(:,{this_lon},{this_lat}) = data(:,15)
;f->epv_gl_s_sun(:,{this_lon},{this_lat})       = data(:,16);;;;;;;;
;f->epv_gl_s_shade(:,{this_lon},{this_lat})     = data(:,17)
;f->epv_m_Kl(:,{this_lon},{this_lat})           = data(:,18)
;f->psn_sun_lnc(:,{this_lon},{this_lat})        = data(:,19)
;f->psn_sun_g(:,{this_lon},{this_lat})          = data(:,20)
;f->psn_sun_Ci(:,{this_lon},{this_lat})         = data(:,21)
;f->psn_sun_Ca(:,{this_lon},{this_lat})         = data(:,22)
;f->psn_sun_A(:,{this_lon},{this_lat})          = data(:,23)
;f->psn_shade_lnc(:,{this_lon},{this_lat})      = data(:,24)
;f->psn_shade_g(:,{this_lon},{this_lat})        = data(:,25)
;f->psn_shade_Ci(:,{this_lon},{this_lat})       = data(:,26)
;f->psn_shade_Ca(:,{this_lon},{this_lat})       = data(:,27)
;f->psn_shade_A(:,{this_lon},{this_lat})        = data(:,28)
;f->summary_daily_npp(:,{this_lon},{this_lat})  = data(:,29)
;f->summary_daily_nep(:,{this_lon},{this_lat})  = data(:,30)
;f->summary_daily_nee(:,{this_lon},{this_lat})  = data(:,31)
;f->summary_daily_gpp(:,{this_lon},{this_lat})  = data(:,32)
;f->summary_daily_mr(:,{this_lon},{this_lat})   = data(:,33)
;f->summary_daily_gr(:,{this_lon},{this_lat})   = data(:,34)
;f->summary_daily_hr(:,{this_lon},{this_lat})   = data(:,35)
;f->summary_vegc(:,{this_lon},{this_lat})       = data(:,36)
;f->summary_litrc(:,{this_lon},{this_lat})      = data(:,37)
;f->summary_soilc(:,{this_lon},{this_lat})      = data(:,38)
;f->summary_totalc(:,{this_lon},{this_lat})     = data(:,39)
;f->metv_prcp(:,{this_lon},{this_lat})          = data(:,40)
;f->metv_tmax(:,{this_lon},{this_lat})          = data(:,41)
;f->metv_tmin(:,{this_lon},{this_lat})          = data(:,42)
;f->metv_tavg(:,{this_lon},{this_lat})          = data(:,43)
;f->metv_tday(:,{this_lon},{this_lat})          = data(:,44)
;f->metv_tnight(:,{this_lon},{this_lat})        = data(:,45)
;f->metv_tsoil(:,{this_lon},{this_lat})         = data(:,46)
;f->metv_vpd(:,{this_lon},{this_lat})           = data(:,47)
;f->cs_litr1c(:,{this_lon},{this_lat})          = data(:,48)
;f->cs_litr2c(:,{this_lon},{this_lat})          = data(:,49)
;f->cs_litr3c(:,{this_lon},{this_lat})          = data(:,50)
;f->cs_litr4c(:,{this_lon},{this_lat})          = data(:,51)
;f->cs_soil1c(:,{this_lon},{this_lat})          = data(:,52)
;f->cs_soil2c(:,{this_lon},{this_lat})          = data(:,53)
;f->cs_soil3c(:,{this_lon},{this_lat})          = data(:,54)
;f->assim_sun(:,{this_lon},{this_lat})          = data(:,55)
;f->assim_shade(:,{this_lon},{this_lat})        = data(:,56)
;f->nleached_snk(:,{this_lon},{this_lat})       = data(:,57);;;;;;;;;
;f->sminn_leached(:,{this_lon},{this_lat})      = data(:,58)
;f->m_tmin(:,{this_lon},{this_lat})             = data(:,59)
;f->m_psi(:,{this_lon},{this_lat})              = data(:,60);;;;;;;;
;f->m_co2(:,{this_lon},{this_lat})              = data(:,61)
;f->m_ppfd_sun(:,{this_lon},{this_lat})         = data(:,62)
;f->m_ppfd_shade(:,{this_lon},{this_lat})       = data(:,63)
;f->m_vpd(:,{this_lon},{this_lat})              = data(:,64)
;f->m_final_sun(:,{this_lon},{this_lat})        = data(:,65)
;f->m_final_shade(:,{this_lon},{this_lat})      = data(:,66)
;f->leaf_psi(:,{this_lon},{this_lat})           = data(:,67);;;;;;;;
;f->m_stem_embolism(:,{this_lon},{this_lat})    = data(:,68)






end

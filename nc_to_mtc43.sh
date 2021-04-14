########################################################################################################################
#shell script to port data
#from GCM .nc output
#to BGC .mtc43
#
#not configured to run as a whole
#several snippets of useful auxiliary code are found throughout
#
#this software comes with no guarantees
#built in OSX, may require troubleshooting on other *nix systems
#
#W J Matthaeus 4.13.21
########################################################################################################################

########set manually
GLACIATION=GLAC
CO2=182
O2=28

GLACIATION=GLAC
CO2=200
O2=21

GLACIATION=GLAC
CO2=600
O2=21

GLACIATION=GLAC
CO2=546
O2=28

#####

GLACIATION=IGLAC
CO2=546
O2=28

GLACIATION=IGLAC
CO2=600
O2=21

##GCM output file names
#LMASK_182.28_41.nc
#TSMIN_182.28_41.2x2.nc
#FRSA_Bilinear_toGrid.nc
#PRECIP_182.28_41.2x2.nc
#TS_182.28_41.2x2.nc
#TSMAX_182.28_41.2x2.nc

#manually make lat.txt and lon.txt 
#by copying from ncdump -h and rearranging into column
#e.g. with find replace "\s+" with "\n", and "," with nothing
touch ${GLACIATION}_${CO2}_${O2}_lat.txt
touch ${GLACIATION}_${CO2}_${O2}_lon.txt
open -a BBEdit ${GLACIATION}_${CO2}_${O2}_lat.txt
open -a BBEdit ${GLACIATION}_${CO2}_${O2}_lon.txt
ncdump -v lat LMASK_${CO2}.${O2}_41.nc
ncdump -v lon LMASK_${CO2}.${O2}_41.nc


##### make landmask text file
#extract LAT,LON pairs that are code 1 (or code for land)
#write them into LAMASK_${CO2}_${O2}_LON_LAT.txt
# CO2=546
# O2=28
mkdir landmask
cd landmask
cp ../${GLACIATION}_${CO2}_${O2}_lon.txt ../${GLACIATION}_${CO2}_${O2}_lat.txt .
#examples for testing
# LATS=11
# LONS=277
touch LAMASK_${GLACIATION}_${CO2}_${O2}_LON_LAT.txt
while read LATS; do
    while read LONS; do
       cdo -s --no_warnings sellonlatbox,${LONS},${LONS},${LATS},${LATS} ../LMASK_${CO2}.${O2}_41.nc LMASK_${GLACIATION}_${CO2}_${O2}_${LONS}_${LATS}.nc
       #ncdump LMASK_${CO2}_${O2}_${LONS}_${LATS}.nc | sed -n -e '/data/,/}/p'| sed '1,6d' | grep ";" | tail -1 | cut -d ' ' -f3 | uniq >> LMASK_${CO2}_${O2}_${LONS}_${LATS}_uniqByLL.txt
       MASK="$(ncdump LMASK_${GLACIATION}_${CO2}_${O2}_${LONS}_${LATS}.nc | sed -n -e '/data/,/}/p'| sed '1,6d' | grep ";" | tail -1 | cut -d ' ' -f3 | uniq)"
       rm LMASK_${GLACIATION}_${CO2}_${O2}_${LONS}_${LATS}.nc
       #echo $MASK
       if [ "${MASK}" = "1" ] 
       then
            #could do the rest of the stuff here in stead of writing a file... or in addition to
            echo "${MASK} ${LONS} ${LATS}">>LAMASK_${GLACIATION}_${CO2}_${O2}_LON_LAT.txt
       fi
    done < ${GLACIATION}_${CO2}_${O2}_lon.txt
done < ${GLACIATION}_${CO2}_${O2}_lat.txt

cat LAMASK_${GLACIATION}_${CO2}_${O2}_LON_LAT.txt | cut -d' ' -f2,3 > ${GLACIATION}_${CO2}_${O2}_LON_LAT.txt
cd ..

#some error checking
#wc -l ${CO2}_${O2}_LON_LAT.txt
#4211 600_21_LON_LAT.txt

##check to make sure coordinate systems etc are the same and correct
##ncdump -h DAYLEN_${CO2}.${O2}_41.2x2.nc
# e.g. variable name change, and remap
# to text file "changeFRSAtime"
# "time@long_name="time";
# time@units="days";
# time@FORTRAN_format="f9.3";"
# for YEAR in $(seq 41 50); do
# 
#     ncap2 -S ../../../changeFRSAtime FRSA_${CO2}.${O2}_${YEAR}.T31.nc FRSA_${CO2}.${O2}_${YEAR}_renamed.nc
# 
#     #regrid FRSA to PRECIP grid
#     cdo griddes PRECIP_${CO2}.${O2}_${YEAR}.2x2.nc > myGridDef
#     cdo remapbil,myGriddef FRSA_${CO2}.${O2}_${YEAR}_renamed.nc FRSA_${CO2}.${O2}_${YEAR}_Bilinear_toGrid.nc
#     rm FRSA_${CO2}.${O2}_${YEAR}_renamed.nc
# 
# done
# # 

#day lengths were created separately, vary with latitude only
cd ..
mkdir daylengths
cd daylengths
cp ../../../../myDAYLEN_Bilinear_toGrid.nc .
while read LAT; do
 cdo sellonlatbox,171,171,${LAT},${LAT} myDAYLEN_Bilinear_toGrid.nc DAYLEN_${GLACIATION}_${CO2}_${O2}_${LAT}.nc
 ncdump DAYLEN_${GLACIATION}_${CO2}_${O2}_${LAT}.nc | sed -n -e '/data/,/}/p'| sed '1,6d' | grep ";" | cut -d ' ' -f3  > LAT_${GLACIATION}_${CO2}_${O2}_${LAT}_DAYLEN.txt
done < ../${GLACIATION}_${CO2}_${O2}_lat.txt
cd ..
# check newlines
# ncdump DAYLEN_${LAT}.nc | sed -n -e '/data/,/}/p'| sed '1,6d' | grep ";" | cut -d ' ' -f3 | perl -pe 'chomp if eof' |perl -pe 's/\n/\t/g' > LAT_${LAT}_DAYLEN.txt
# tab separated
# ncdump DAYLEN_${LAT}.nc | sed -n -e '/data/,/}/p'| sed '1,6d' | grep ";" | cut -d ' ' -f3 | perl -pe 'chomp if eof' |perl -pe 's/\n/\t/g' > LAT_${LAT}_DAYLEN.txt



############################################################
#loops below
#output text of the required variable to text
#extra info in this fomat...
#process using clt to just array
#just the data in one column
# for all the coordinates in  LON_LAT.txt
#separate steps can be broken out by uncommenting inner loops
############################################################

#run in main co2 o2 dir
#intermediate 'done' points left in for troubleshooting.
while read LON_LAT; do
 LON="$(echo ${LON_LAT} | cut -d' ' -f1)"
 LAT="$(echo ${LON_LAT} | cut -d' ' -f2)"
# echo ${LON}
# echo ${LAT}
#done < LON_LAT.txt
 mkdir ${LON}_${LAT}
 cd ${LON}_${LAT}

for YEAR in $(seq 41 50); do


 cdo -s --no_warnings sellonlatbox,${LON},${LON},${LAT},${LAT} ../PRECIP_${CO2}.${O2}_${YEAR}.2x2.nc PRECIP_${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR}
 cdo -s --no_warnings sellonlatbox,${LON},${LON},${LAT},${LAT} ../FRSA_${CO2}.${O2}_${YEAR}_Bilinear_toGrid.nc FRSA_${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR}
 cdo -s --no_warnings sellonlatbox,${LON},${LON},${LAT},${LAT} ../TS_${CO2}.${O2}_${YEAR}.2x2.nc TS_${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR}
 cdo -s --no_warnings sellonlatbox,${LON},${LON},${LAT},${LAT} ../TSMIN_${CO2}.${O2}_${YEAR}.2x2.nc TSMIN_${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR}
 cdo -s --no_warnings sellonlatbox,${LON},${LON},${LAT},${LAT} ../TSMAX_${CO2}.${O2}_${YEAR}.2x2.nc TSMAX_${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR}

# done < LON_LAT.txt
# while read LON_LAT; do
#  LON="$(echo ${LON_LAT} | cut -d' ' -f1)"
#  LAT="$(echo ${LON_LAT} | cut -d' ' -f2)"

 ncdump -v PRECIP PRECIP_${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR} | sed -n -e '/data/,/;/p' | sed '1,3d' | cut -d',' -f1 | sed 's/;//' | sed 's/^  //' > ${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR}__5_PRECIP
 ncdump -v FRSA FRSA_${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR} | sed -n -e '/data/,/;/p'| sed '1,3d' | cut -d',' -f1 | sed 's/;//' | sed 's/^  //'  > ${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR}__7_FRSA
 ncdump -v TS TS_${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR} | sed -n -e '/data/,/;/p'| sed '1,3d' | cut -d',' -f1 | sed 's/;//' | sed 's/^  //'  > ${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR}__4_TS
 ncdump -v TS_M TSMIN_${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR} | sed -n -e '/data/,/;/p'| sed '1,3d' | cut -d',' -f1 | sed 's/;//' | sed 's/^  //'  > ${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR}__3_TSMIN
 ncdump -v TS_M TSMAX_${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR} | sed -n -e '/data/,/;/p'| sed '1,3d' | cut -d',' -f1 | sed 's/;//' | sed 's/^  //'  > ${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR}__2_TSMAX

# done < LON_LAT.txt
# while read LON_LAT; do
#     LON="$(echo ${LON_LAT} | cut -d' ' -f1)"
#     LAT="$(echo ${LON_LAT} | cut -d' ' -f2)"

#assemble, uses weird number in name for order... lol
 paste ${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR}__* > ${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR}_no_YEARDAY_VPD.tsv
##add yearday
 awk '{ print NR"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6 }' ${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR}_no_YEARDAY_VPD.tsv > ${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR}_noVPD.tsv
# done < LON_LAT.txt
# while read LON_LAT; do
#     LON="$(echo ${LON_LAT} | cut -d' ' -f1)"
#     LAT="$(echo ${LON_LAT} | cut -d' ' -f2)"
##############Change path
 Rscript --vanilla /path/to/makingMtclim.R ${GLACIATION}_${CO2}_${O2}_${LON}_${LAT}_${YEAR}_noVPD.tsv

done
 cd ..
done < landmask/${GLACIATION}_${CO2}_${O2}_LON_LAT.txt


### remove everything you just wrote
# while read LON_LAT; do
#      LON="$(echo ${LON_LAT} | cut -d' ' -f1)"
#      LAT="$(echo ${LON_LAT} | cut -d' ' -f2)"
#      rm -r ${LON}_${LAT}
# done < ${GLACIATION}_${CO2}_${O2}_LON_LAT.txt

## plant sims
# these variables set superdirectory,target simulation set, column in daily output files,
#and output filenames
DATE=01142019
OUTVAR=soilw_outflow
OUTVARNUM=7
GLACIATION=GLAC
CO2=182
O2=28
cd $DATA


for i in grid_${GLACIATION}_${CO2}_${O2}_${SIMDATE}/*/outputs/*.dayout.ascii; do 
    awk -v var="$OUTVARNUM" 'FNR==NR {sum+=$var; i+=1; next}; {M2+=($var-(sum/i))^2; M3 +=($var-(sum/i))^3; M4 +=($var-(sum/i))^4;}
    END {if (M2!=0) {print FILENAME , NR, i , M2 , M3, M4, sum/i, M3/(i*sqrt(M2/i)^3) , M4/(i*(M2/i)^2)-3} else {print FILENAME , NR, i , M2 , M3, M4, sum/i}}' $i $i  >> Stats8_${OUTVAR}_${GLACIATION}_${CO2}_${O2}_${DATE}.txt
done



GLACIATION=IGLAC
CO2=546
O2=28
cd $DATA


for i in grid_${GLACIATION}_${CO2}_${O2}_${SIMDATE}/*/outputs/*.dayout.ascii; do 
    awk -v var="$OUTVARNUM" 'FNR==NR {sum+=$var; i+=1; next}; {M2+=($var-(sum/i))^2; M3 +=($var-(sum/i))^3; M4 +=($var-(sum/i))^4;}
    END {if (M2!=0) {print FILENAME , NR, i , M2 , M3, M4, sum/i, M3/(i*sqrt(M2/i)^3) , M4/(i*(M2/i)^2)-3} else {print FILENAME , NR, i , M2 , M3, M4, sum/i}}' $i $i  >> Stats8_${OUTVAR}_${GLACIATION}_${CO2}_${O2}_${DATE}.txt
done



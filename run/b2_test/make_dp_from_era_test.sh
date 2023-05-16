#!/bin/bash
#Script to dgrib era5 data - from Sean W. Freeman, augmented by Peter J. Marinescu

# Path to dgrib code
cd /home6/grleung1/basinlcc-rams/bin.dp.grib1/

# ERA5 Data input
IN_FILE="/nobackupp12/grleung1/basinlcc/b2_test/ERA5/ERA5_20190911.grib"

# DP_FILE Output directory needs trailing slash
OUT_DIR="/nobackupp12/grleung1/basinlcc/b2_test/dp/"

YEAR=2019 # YYYY
MONTH=09 # MM
START_DAY=11 # DD
END_DAY=13 # DD

#### UPDATE INFOMRATION ABOVE HERE ^^^^^

for day in $(seq -f "%02g" $START_DAY $END_DAY)
do
    for hour in $(seq -f "%02g" 0 23)
    do
        # -h 0 current analysis time (forecase hour) 
        # -d specifies the input date to pull from the file
        # -t 6 specifies to use ERA5 with specific humidity
        # -f specifies in the input grid file
        ./dgrib-6.3.02 -t 6 -d $YEAR$MONTH$day$hour -h 0 -f $IN_FILE
        mv dp-p$YEAR-$MONTH-$day-"$hour"00* $OUT_DIR
    done
done

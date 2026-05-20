#!/bin/bash

less P18.onlyPASS.maf | sed 's/tumor/tumor_18/g' | sed 's/normal/normal_18/g' > P18.onlyPASS_labled.maf 
less P19.onlyPASS.maf | sed 's/tumor/tumor_19/g' | sed 's/normal/normal_19/g' > P19.onlyPASS_labled.maf

sed -e "1,2d" P19.onlyPASS_labled.maf > P19.onlyPASS_labled-2.maf

cat P18.onlyPASS_labled.maf P19.onlyPASS_labled-2.maf > Both.onlyPASS_labled.maf

rm P19.onlyPASS_labled-2.maf
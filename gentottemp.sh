grep "TOTAL TEMP" $1 |perl -ane  '$i++; print "$i @F[-2]\n"' > tottemp

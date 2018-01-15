perl -ane '$i = 1 if /STATES \[/; $i = 0 if /DECOMPOSED_/; print "@F[4]\n" if $i' $1 > state
perl -ne 'if(/parent/) {print "\n";} else {chomp; print "$_ ";}' state > dooooooooooo
tail -n +2 dooooooooooo > state
pr -t -m -w 2000 timestep state > dooooooooooo
perl -ane 'for(@F) {print "$_ " ;} print "\n"' dooooooooooo > state.dat

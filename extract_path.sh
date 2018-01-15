perl -ane '$i = 1 if /STATES \[/; $i = 0 if /DECOMPOSED_/; print "@F[6]\n" if $i' $1 > path
perl -ne 'if(/shell/) {print "\n";} else {chomp; print "$_ ";}' path > dooooooooooo
tail -n +2 dooooooooooo > path
pr -t -m -w 2000 timestep path > dooooooooooo
perl -ane 'for(@F) {print "$_ " ;} print "\n"' dooooooooooo > path.dat

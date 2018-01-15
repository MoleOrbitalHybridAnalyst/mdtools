perl -ane '$i = 1 if /STATES \[/; $i = 0 if /DECOMPOSED_/; print "@F[2]\n" if $i' $1 > shell
perl -ne 'if(/id/) {print "\n";} else {chomp; print "$_ ";}' shell > dooooooooooo
tail -n +2 dooooooooooo > shell
pr -t -m -w 2000 timestep shell > dooooooooooo
perl -ane 'for(@F) {print "$_ " ;} print "\n"' dooooooooooo > shell.dat

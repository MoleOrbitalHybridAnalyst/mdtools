# read evb.out
# use ~/share/read_ci_fromevb.py 
# all the hard-coded things are in ~/share/read_ci_fromevb.py
evbout=$1
# do eigen
grep -A 1 EIGEN_VEC $evbout |grep -v EIGEN|grep -v "\-\-" > eigen_vector
grep  TIMESTEP $evbout |cut -d' ' -f2 > timestep
pr -t -m -w 2000 timestep eigen_vector > doo
perl -ane 'for(@F) {print "$_ " ;} print "\n"' doo > eigen_vector.dat
echo "done with eigen_vector.dat"
# do state
perl -ane '$i = 1 if /STATES \[/; $i = 0 if /DECOMPOSED_/; print "@F[4]\n" if $i' $evbout > state
perl -ne 'if(/parent/) {print "\n";} else {chomp; print "$_ ";}' state > doo
tail -n +2 doo > state
pr -t -m -w 500 timestep state > doo
perl -ane 'for(@F) {print "$_ " ;} print "\n"' doo > state.dat
echo "done with state.dat"
# do shell
perl -ane '$i = 1 if /STATES \[/; $i = 0 if /DECOMPOSED_/; print "@F[2]\n" if $i' $evbout > shell
perl -ne 'if(/id/) {print "\n";} else {chomp; print "$_ ";}' shell > doo
tail -n +2 doo > shell
pr -t -m -w 500 timestep shell > doo
perl -ane 'for(@F) {print "$_ " ;} print "\n"' doo > shell.dat
echo "done with shell.dat"
# do path
perl -ane '$i = 1 if /STATES \[/; $i = 0 if /DECOMPOSED_/; print "@F[6]\n" if $i' $evbout > path
perl -ne 'if(/shell/) {print "\n";} else {chomp; print "$_ ";}' path > doo
tail -n +2 doo > path
pr -t -m -w 500 timestep path > doo
perl -ane 'for(@F) {print "$_ " ;} print "\n"' doo > path.dat
echo "done with path.dat"
python3 ~/share/read_ci_fromevb.py 
echo "ci_evb.dat generated"

grep -A 1 EIGEN_VEC $1 |grep -v EIGEN|grep -v "\-\-" > eigen_vector
grep TIMESTEP $1 | cut -d' ' -f2 > timestep
pr -t -m -w 2000 timestep eigen_vector > dooooooooooo
perl -ane 'for(@F) {print "$_ ";} print "\n"' dooooooooooo > eigen_vector.dat

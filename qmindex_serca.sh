grep "MM_INDEX" qm_kind_cell_link |perl -ane '$i++; $_=$_-1 for @F[1..@F-1]; print join " ", @F[1..@F-1]; print "\n"; last if $i>4;'

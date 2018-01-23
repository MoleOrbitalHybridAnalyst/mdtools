#cat $1|perl -ane 'print " -block @F[0]/COLVAR -bxy 1:3 "'|xargs xmgrace
ls $1*/COLVAR |perl -ane 'print " -block $F[0] -bxy 1:3 "'|xargs xmgrace

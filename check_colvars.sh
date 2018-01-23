# sh this.sh dir_list column
if [ $2 -eq 2 ] 
then
   perl -ane 'print "$F[0]/COLVAR "' $1|xargs xmgrace -legend load
else
   perl -ane 'print " -block $F[0]/COLVAR -bxy 1:'$2'"' $1 |xargs xmgrace 
fi

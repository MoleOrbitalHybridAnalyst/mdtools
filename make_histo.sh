for i in `cat dir`
do
   perl ~/share/makehisto2.pl 20 1 $i/COLVAR > histos/$i.histo
done

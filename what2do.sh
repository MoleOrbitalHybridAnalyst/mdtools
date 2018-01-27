# sh this.sh dir.1d dir.2d
for i in `cat $1`
do
   # --- find the largest and smallest dir that has COLVAR
   max=125; min=125
   max=`grep --perl "^$i[-|$]" $2| while read j; do cd $j; if [ -e COLVAR ]; then echo $j; fi; cd ..; done|cut -d'-' -f2 | while read sol; do if [ $sol != $i ]; then if [ $sol -gt $max ]; then max=$sol; echo $max; fi; fi; done|tail -n1`
   min=`grep --perl "^$i[-|$]" $2| while read j; do cd $j; if [ -e COLVAR ]; then echo $j; fi; cd ..; done|cut -d'-' -f2 | while read sol; do if [ $sol != $i ]; then if [ $min -gt $sol ]; then min=$sol; echo $min; fi; fi; done|tail -n1`
   if [ ! $max ] 
   then
      max=125
   fi
   if [ ! $min ] 
   then 
      min=125
   fi
   max_dir=$i; min_dir=$i
   if [ $max -eq $min ]
   then
      if [ ! -e $i/COLVAR ]
      then
         echo "sh sdo.sh $i"
         continue
      fi
   else
      max_dir=$i-$max
      min_dir=$i-$min
   fi
   max_next=$(expr $max + 5)
   min_next=$(expr $min - 5)
   test_dir="$i-$max_next"
   # check whehter the next step dir exists
   if [ -e $test_dir ]
   then
      cd $test_dir
      # check whether inputs are prepared in the next step dir
      if [ -e restart ]
      then
         tmp=`squeue -u chhli -n $test_dir|wc -l`
         if [  $tmp -eq 1 ]
         then
            printf "sh sdo.sh $test_dir; "
         else
            if [ $tmp -gt 2 ]
            then
               echo "more than 1 $test_dir!!!!"
            fi
         fi
      else
         printf "sh copy_restart.sh $max_dir $test_dir; "
      fi
      cd ..
   else
      printf "mkdir $test_dir; "
   fi
   test_dir="$i-$min_next"
   if [ -e $test_dir ]
   then
      cd $test_dir
      if [ -e restart ]
      then
         tmp=`squeue -u chhli -n $test_dir|wc -l`
         if [ $tmp -eq 1 ]
         then
            printf "sh sdo.sh $test_dir; "
         else
            if [ $tmp -gt 2 ]
            then
               echo "more than 1 $test_dir!!!!"
            fi
         fi
      else
         printf "sh copy_restart.sh $min_dir $test_dir; "
      fi
      cd ..
   else
      printf "mkdir $test_dir"
   fi
   echo ""
done

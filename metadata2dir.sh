cat $1 | cut -d ' ' -f1  | perl  -F"\/" -ane 'print join "\/", @F[2..$#F-1]; print "\n"'

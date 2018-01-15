set nf [molinfo 0 get numframes]
for {set i 0} {$i < $nf} {incr i} {
   set a [atomselect 0 "all" frame $i]
   $a writepdb "[lindex $argv 0]/${i}.pdb"
}
exit

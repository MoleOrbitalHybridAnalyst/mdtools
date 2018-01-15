#calculate min tip distances between TM1,2 and TM7,8 ; TM4,5 and TM10,11
#read trj from command line arguments
set nf [molinfo 0 get numframes]
set fp [open [lindex $argv 16] w]
set tm12 [atomselect 0 "name CA and (resid >=[lindex $argv 0] and resid <=[lindex $argv 1] or resid >=[lindex $argv 2] and resid <=[lindex $argv 3])"]
set tm78 [atomselect 0 "name CA and (resid >=[lindex $argv 4] and resid <=[lindex $argv 5] or resid >=[lindex $argv 6] and resid <=[lindex $argv 7])"]
set tm45 [atomselect 0 "name CA and (resid >=[lindex $argv 8] and resid <=[lindex $argv 9] or resid >=[lindex $argv 10] and resid <=[lindex $argv 11])"]
set tm1011 [atomselect 0 "name CA and (resid >=[lindex $argv 12] and resid <=[lindex $argv 13] or resid >=[lindex $argv 14] and resid <=[lindex $argv 15])"]
for {set i 0} {$i < $nf} {incr i} {
    #update selection
    $tm12 frame $i
    $tm78 frame $i
    $tm45 frame $i
    $tm1011 frame $i
    #intial min
    set min1 1000
    set min2 1000
    set x12s [$tm12 get {x y z}] 
    set x78s [$tm78 get {x y z}] 
    set x45s [$tm45 get {x y z}] 
    set x1011s [$tm1011 get {x y z}] 
    foreach x12 $x12s {
        foreach x78 $x78s {
            set dist [vecdist $x12 $x78]
            if {$dist < $min1} {
                set min1 $dist
            }
        }
    }
    foreach x45 $x45s {
        foreach x1011 $x1011s {
            set dist [vecdist $x45 $x1011]
            if {$dist < $min2} {
                set min2 $dist
            }
        }
    }
    puts $fp "$i\t$min1\t$min2"
}
exit

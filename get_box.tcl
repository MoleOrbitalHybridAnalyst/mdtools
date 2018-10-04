set nf [molinfo 0 get numframes]
set fp [open [lindex $argv 0] w]
for {set i 0} {$i < $nf} {incr i} {
   animate goto $i
   set pbc_info [lindex [pbc get] 0]
   puts $fp "$i [lindex $pbc_info 0] [lindex $pbc_info 1] [lindex $pbc_info 2]"
}
close $fp
exit

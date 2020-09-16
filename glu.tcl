mol representation DynamicBonds 1.3 0.1 12.0
mol representation VDW 0.2 12.0
mol showrep 0 0 off
mol addrep 0
mol modselect 1 0 protein
mol modstyle 1 0 Licorice
mol addrep 0
mol modselect 2 0 not protein and within 8 of (protein and sidechain and oxygen)
mol modstyle 2 0 DynamicBonds
mol selupdate 2 0 1
mol addrep 0
mol modselect 3 0 not protein and within 8 of (protein and sidechain and oxygen)
mol modstyle 3 0 VDW
mol selupdate 3 0 1
mol addrep 0
mol modselect 4 0 serial [expr [molinfo 0 get numatoms] + 1]
display projection orthographic
#color Display Background white
#color Labels Bonds black
source ~/share/vmdtools.tcl

set nf [molinfo 0 get numframes]
for {set i 0} {$i < $nf} {incr i} {
   animate goto $i
   set mm [measure minmax [atomselect 0 all]]
   set xmean [expr ([lindex [lindex $mm 0] 0] + [lindex [lindex $mm 1] 0]) / 2]
   set ymean [expr ([lindex [lindex $mm 0] 1] + [lindex [lindex $mm 1] 1]) / 2]
   set zmean [expr ([lindex [lindex $mm 0] 2] + [lindex [lindex $mm 1] 2]) / 2]
   set cen [geo_center [atomselect 0 "name OE OC"]]
   set shift {}
   lappend shift [expr $xmean - [lindex $cen 0]]
   lappend shift [expr $ymean - [lindex $cen 1]]
   lappend shift [expr $zmean - [lindex $cen 2]]
   [atomselect 0 all] moveby $shift
   pbc wrap
}

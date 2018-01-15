mol new PepTXc.pdb
#0 trr ;1 xrystall
set fp [open "rmsd.dat" w]
set nf [molinfo 0 get numframes]
for {set i 0} {$i < $nf} {incr i} {
    set out1 [atomselect 0 "name CA and not (resid 22 136 ) and not within 20 of (resid 425 322 and protein)" frame $i]
    set out2 [atomselect 1 "name CA and resid [$out1 get resid]" frame $i]
    set tran [measure fit $out2 $out1]
    set move [atomselect 1 "all" frame $i]
    $move move $tran
    set sel1 [atomselect 0 "mass>2 and resid 425 322 and protein" frame $i]
    set sel2 [atomselect 1 "resid 425 322" frame $i]
    puts $fp "$i [measure rmsd $sel1 $sel2]"
}
exit

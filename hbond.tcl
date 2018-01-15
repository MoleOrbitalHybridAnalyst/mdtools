#vmd-text-mode top.psf trj.xyz -e hbond.tcl -args resid1 resid2 xl xh yl yh zl zh startingframe outfile
#carboxylic oxygen must have type OC.*
set nf [molinfo 0 get numframes]
#some paras here
#only waters within 3A be regarded as acceptors
#set dist 3.5
set dist 4
#set angle 30
set angle 90
set start_resid [lindex $argv 0]
set end_resid [lindex $argv 1]
#only waters within this cube be considered
set xmin [lindex $argv 2]
set xmax [lindex $argv 3]
set ymin [lindex $argv 4]
set ymax [lindex $argv 5]
set zmin [lindex $argv 6]
set zmax [lindex $argv 7]
set start_frame [lindex $argv 8]
set fpout [open [lindex $argv 9] w]
set prev_sel [atomselect 0 [concat "resid $start_resid" { and type "OD.*" OE1 OE2}] ]
set end_id [[atomselect 0 [concat "resid $end_resid and " {type "OD.*" OE1 OE2}]] list]
for {set i $start_frame} {$i < $nf} {incr i} {
    $prev_sel frame $i
    set prev_list [$prev_sel list]
    #only same residue of (water and hydronium or ending points and name O.* within dist+1 of prev_sel) is selected here 
    #meta_sel also includes the prev_sel
    set meta_sel [atomselect 0 [concat {name "O.*"} " and ((waters or resname H3O or index $end_id) and within [expr $dist + 1] of index $prev_list and x<$xmax and x>$xmin and y<$ymax and y>$ymin and z<$zmax and z>$zmin)"] frame $i]
    if {[$meta_sel list] == {}} {
        puts $fpout "$i\t0"
        continue
    }
    set meta_sel [atomselect 0  "index $prev_list or same residue as index [$meta_sel list]" frame $i]
    set hbond [measure hbonds $dist $angle $meta_sel]
    if {[lindex $hbond 0] == {}} {
        puts $fpout "$i\t0"
        continue
    }
    #select the residues excluding prev_sel to make sure it is going forwards
    set cnt_sel [atomselect 0 "(same residue as (index [concat [lindex $hbond 0] [lindex $hbond 1]]) and not index $prev_list)" frame $i]
    set cnt_list [$cnt_sel list]
    #loop for 20 times
    for {set j 0} {$j < 20} {incr j} {
        #if new cnt_list is empty, all the hbonds are within old cnt_sel, the old cnt_sel is isolated, break
        if {$cnt_list == {}} {
            puts $fpout "$i\t0"
            break
        }
        #if new current list has ending points break
        set itsct [atomselect 0 "index $cnt_list and index $end_id" frame $i]
        if {[$itsct num]} {
            puts $fpout "$i\t1"
            break
       } 
        #use cnt_sel to construct meta_sel
        #select H2O,H3O+,target residue within dist+1 of cnt_sel, include cnt_sel itself, exclude prev_sel 
        set meta_sel [atomselect 0 [concat {name "O.*"} " and (waters or resname H3O or index $end_id) and within [expr $dist + 1] of index $cnt_list and x<$xmax and x>$xmin and y<$ymax and y>$ymin and z<$zmax and z>$zmin"] frame $i]
    if {[$meta_sel list] == {}} {
        puts $fpout "$i\t0"
        break
    }
        set meta_sel [atomselect 0 "index $cnt_list or ((same residue as index [$meta_sel list]) and not index $prev_list)" frame $i]
        set hbond [measure hbonds $dist $angle $meta_sel]
        if {[lindex $hbond 0] == {}} {
            puts $fpout "$i\t0"
            break
        }
        #now, there is something in hbond, it could just be the cnt atoms (impossible for prev atoms as we exclude them when selecting meta_sel). We'd like to exclude them 
        #push the old cnt atoms to prev atom list
        set prev_list [concat $prev_list $cnt_list]
        set cnt_sel [atomselect 0 "(same residue as (index [concat [lindex $hbond 0] [lindex $hbond 1]]) and not index $cnt_list)" frame $i]
        #update cnt_list
        set cnt_list [$cnt_sel list]
    }
}
close $fpout 
exit

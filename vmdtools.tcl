proc draw_arrow {mol start end {color blue} {radius1 0.15} {radius2 0.25}} {
   set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
   graphics $mol color $color
   graphics $mol cylinder $start $middle radius $radius1
   graphics $mol cone $middle $end radius $radius2
}

# maybe duplicated in measure pbcneighbors
proc pbcwithin {serial range bx by bz} {
   lassign [lindex [[atomselect 0 "serial $serial"] get {x y z}] 0] a b c
   set s [atomselect 0 "within $range of serial $serial"]
   set r2 [expr $range * $range]
   set s [atomselect 0 "index [$s list] or (x-$a-$bx)^2+(y-$b)^2+(z-$c)^2<=$r2"]
   set s [atomselect 0 "index [$s list] or (x-$a+$bx)^2+(y-$b)^2+(z-$c)^2<=$r2"]
   set s [atomselect 0 "index [$s list] or (x-$a)^2+(y-$b+$by)^2+(z-$c)^2<=$r2"]
   set s [atomselect 0 "index [$s list] or (x-$a)^2+(y-$b-$by)^2+(z-$c)^2<=$r2"]
   set s [atomselect 0 "index [$s list] or (x-$a)^2+(y-$b)^2+(z-$c+$bz)^2<=$r2"]
   set s [atomselect 0 "index [$s list] or (x-$a)^2+(y-$b)^2+(z-$c-$bz)^2<=$r2"]
   return [$s list]
}

proc draw_box {mol center vx vy vz xl xh yl yh zl zh {color blue} {width 1}} {
   set dlll [vecadd [vecscale $xl $vx] [vecadd [vecscale $yl $vy] [vecadd [vecscale $zl $vz] $center]]]
   set dhll [vecadd [vecscale $xh $vx] [vecadd [vecscale $yl $vy] [vecadd [vecscale $zl $vz] $center]]]
   set dlhl [vecadd [vecscale $xl $vx] [vecadd [vecscale $yh $vy] [vecadd [vecscale $zl $vz] $center]]]
   set dhhl [vecadd [vecscale $xh $vx] [vecadd [vecscale $yh $vy] [vecadd [vecscale $zl $vz] $center]]]
   set dllh [vecadd [vecscale $xl $vx] [vecadd [vecscale $yl $vy] [vecadd [vecscale $zh $vz] $center]]]
   set dhlh [vecadd [vecscale $xh $vx] [vecadd [vecscale $yl $vy] [vecadd [vecscale $zh $vz] $center]]]
   set dlhh [vecadd [vecscale $xl $vx] [vecadd [vecscale $yh $vy] [vecadd [vecscale $zh $vz] $center]]]
   set dhhh [vecadd [vecscale $xh $vx] [vecadd [vecscale $yh $vy] [vecadd [vecscale $zh $vz] $center]]]
   graphics $mol color $color
   graphics $mol line $dlll $dhll width $width
   graphics $mol line $dlll $dlhl width $width
   graphics $mol line $dlhl $dhhl width $width
   graphics $mol line $dhll $dhhl width $width
   graphics $mol line $dlll $dllh width $width
   graphics $mol line $dhll $dhlh width $width
   graphics $mol line $dlhl $dlhh width $width
   graphics $mol line $dhhl $dhhh width $width
   graphics $mol line $dllh $dhlh width $width
   graphics $mol line $dllh $dlhh width $width
   graphics $mol line $dlhh $dhhh width $width
   graphics $mol line $dhlh $dhhh width $width
}

proc geo_center {select} {
   set w [expr 1.0 / [$select num]]
   set com {0 0 0}
   foreach x [$select get {x y z}] {
      set com [vecadd $com [vecscale $w $x]]
   }
   return $com
}

proc draw_path {mol path_file {color orange} {radius 0.2} {resolution 10}} {
   set fp [open $path_file r]
   set path_filedat [read $fp]
   close $fp
   graphics $mol color $color
   global path_dat
   set path_dat [split $path_filedat "\n"]
   foreach line $path_dat {
      if {[expr {$line eq ""}]} {
         continue
      }
      set splits [split $line]
      graphics $mol sphere [list [lindex $splits 0] [lindex $splits 1] [lindex $splits 2]] radius $radius resolution $resolution
   }
}

proc draw_path_w_cen {mol cen path_offset_file {color orange} {radius 0.2} {resolution 10}} {
   set fp [open $path_offset_file r]
   set path_offset_filedat [read $fp]
   close $fp
   graphics $mol color $color
   set path_offset_dat [split $path_offset_filedat "\n"]
   foreach line $path_offset_dat {
      if {[expr {$line eq ""}]} {
         continue
      }
      set splits [split $line]
      set offset [list [lindex $splits 0] [lindex $splits 1] [lindex $splits 2]]
      graphics $mol sphere [vecadd $cen $offset] radius $radius resolution $resolution
   }
}

# http://www.ks.uiuc.edu/Training/Tutorials/vmd-imgmv/imgmv/tutorial-html/node3.html
proc enable_trace_path {} {
   global vmd_frame
   trace variable vmd_frame(0) w draw_path_counter
}

proc disable_trace_path {} {
   global vmd_frame
   trace vdelete vmd_frame(0) w draw_path_counter
}

proc draw_path_counter { name element op } {
   global vmd_frame
   global path_dat
   graphics 0 delete all
   graphics 0 color orange
   set linecount 0
   foreach line $path_dat {
      if {[expr {$line eq ""}]} {
         continue
      }
      if {$linecount eq $vmd_frame(0)} {
         set splits [split $line]
         graphics 0 sphere [list [lindex $splits 0] [lindex $splits 1] [lindex $splits 2]] radius 0.5 resolution 10
         break
      }
      incr linecount
   }
}

proc show_atom_counter {} {}

# http://www.ks.uiuc.edu/Training/Tutorials/vmd-imgmv/imgmv/tutorial-html/node3.html
proc enable_trace_texts {} {
   global vmd_frame
   trace variable vmd_frame(0) w draw_text_counter
}

proc disable_trace_texts {} {
   global vmd_frame
   trace vdelete vmd_frame(0) w draw_text_counter
}

proc draw_text_counter { name element op } {
   global vmd_frame
   global text_dat
   graphics 0 delete all
   set linecount 0
   foreach line $text_dat {
      if {[expr {$line eq ""}]} {
         continue
      }
      if {$linecount eq $vmd_frame(0)} {
         set splits [split $line]
#         graphics 0 sphere [list [lindex $splits 0] [lindex $splits 1] [lindex $splits 2]] radius 0.5 resolution 10
         if {[lindex $splits 3] > 0.0} {
            graphics 0 color red
         } else {
            graphics 0 color blue
         }
         set aux [list [lindex $splits 0] [lindex $splits 1] [lindex $splits 2]]
         graphics 0 text $aux [lindex $splits 3] size 1.4 thickness 2.8
         break
      }
      incr linecount
   }
}

proc draw_texts {mol text_file {size 1.4} {thickness 2.8}} {
   set fp [open $text_file r]
   set text_filedat [read $fp]
   close $fp
   global text_dat
   set text_dat [split $text_filedat "\n"]
   foreach line $text_dat {
      if {[expr {$line eq ""}]} {
         continue
      }
      set splits [split $line]
      if {[lindex $splits 3] > 0.0} {
         graphics $mol color red
      } else {
         graphics $mol color blue
      }
      set aux [list [lindex $splits 0] [lindex $splits 1] [lindex $splits 2]]
      graphics $mol text $aux [lindex $splits 3] size $size thickness $thickness
   }
}

proc draw_texts_u_serials {mol text_file {size 1.4} {thickness 2.8}} {
   set fp [open $text_file r]
   set text_filedat [read $fp]
   close $fp
   set text_dat [split $text_filedat "\n"]
   foreach line $text_dat {
      if {[expr {$line eq ""}]} {
         continue
      }
      set splits [split $line]
      if {[lindex $splits 1] > 0.0} {
         graphics $mol color red
      } else {
         graphics $mol color blue
      }
      set aux [[atomselect 0 "serial [lindex $splits 0]"] get {x y z}]
      set aux [lindex $aux 0]
      graphics $mol text $aux [lindex $splits 1] size $size thickness $thickness
   }
}

proc draw_arrows_u_serials {mol arrow_file {color red} {radius1 0.1} {radius2 0.15}} {
   set fp [open $arrow_file r]
   set arrow_filedat [read $fp]
   close $fp
   set arrow_dat [split $arrow_filedat "\n"]
   foreach line $arrow_dat {
      if {[expr {$line eq ""}]} {
         continue
      }
      set splits [split $line]
      set start [lindex [[atomselect $mol "serial [lindex $splits 0]"] get {x y z}] 0]
      set v [lindex $splits 1]
      lappend v [lindex $splits 2]
      lappend v [lindex $splits 3]
      set end [vecadd $start $v]
      draw_arrow $mol $start $end $color $radius1 $radius2
   }
}

proc draw_arrows_u_coords {mol arrow_file {color red} {radius1 0.1} {radius2 0.15}} {
   set fp [open $arrow_file r]
   set arrow_filedat [read $fp]
   close $fp
   set arrow_dat [split $arrow_filedat "\n"]
   foreach line $arrow_dat {
      if {[expr {$line eq ""}]} {
         continue
      }
      set splits [split $line]
      set start [lindex $splits 0]
      lappend start [lindex $splits 1]
      lappend start [lindex $splits 2]
      set v [lindex $splits 3]
      lappend v [lindex $splits 4]
      lappend v [lindex $splits 5]
      set end [vecadd $start $v]
      draw_arrow $mol $start $end $color $radius1 $radius2
   }
}

# write id(serial) type mol(resid) x y z vx vy vz
proc write_lammpstrj {sel file_name {timestep 0}} {
   set ids [$sel get serial]
   set types [$sel get type]
   set mols [$sel get resid]
   set xs [$sel get x]
   set ys [$sel get y]
   set zs [$sel get z]
   set vxs [$sel get vx]
   set vys [$sel get vy]
   set vzs [$sel get vz]
   set fp [open $file_name w]
   puts $fp "ITEM: TIMESTEP"
   puts $fp $timestep
   puts $fp "ITEM: NUMBER OF ATOMS"
   set natoms [$sel num]
   puts $fp $natoms
   puts $fp "ITEM: BOX BOUNDS pp pp pp"
   set center [geo_center $sel]
   set box {}
   lappend box [lindex [pbc get -namd] 0 0 0]
   lappend box [lindex [pbc get -namd] 0 1 1]
   lappend box [lindex [pbc get -namd] 0 2 2]
   set lo [vecsub $center [vecscale 0.5 $box]]
   set hi [vecadd $center [vecscale 0.5 $box]]
   for {set i 0} {$i < 3} {incr i} {
      puts -nonewline $fp "[lindex $lo $i] "
      puts $fp [lindex $hi $i]
   }
   puts $fp "ITEM: ATOMS id type mol x y z vx vy vz"
   for {set i 0} {$i < $natoms} {incr i} {
      puts -nonewline $fp "[lindex $ids $i] "
      puts -nonewline $fp "[lindex $types $i] "
      puts -nonewline $fp "[lindex $mols $i] "
      puts -nonewline $fp "[lindex $xs $i] "
      puts -nonewline $fp "[lindex $ys $i] "
      puts -nonewline $fp "[lindex $zs $i] "
      puts -nonewline $fp "[lindex $vxs $i] "
      puts -nonewline $fp "[lindex $vys $i] "
      puts $fp [lindex $vzs $i]
   }
   close $fp
}

# why does not measure have this kind of thing ?
proc measure_distances {sel0 sel1 {frame ""} } {
   set l0 [$sel0 list]
   set l1 [$sel1 list]
   set distances {}
   foreach i0 $l0 {
      foreach i1 $l1 {
         set pair $i0
         lappend pair $i1
         if {$frame eq ""} {
            lappend distances [measure bond $pair]
         } else {
            lappend distances [measure bond $pair frame $frame]
         }
      }
   }
   return $distances
}

# this kind of thing does not exist ???
proc min_value {l} {
   set min [lindex $l 0]
   foreach e $l {
      if {[expr {$e < $min}]} {
         set min $e
      }
   }
   return $min
}
proc max_value {l} {
   set max [lindex $l 0]
   foreach e $l {
      if {[expr {$e > $max}]} {
         set max $e
      }
   }
   return $max
}

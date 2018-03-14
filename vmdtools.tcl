proc draw_arrow {mol start end {color blue} {radius1 0.15} {radius2 0.25}} {
   set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
   graphics $mol color $color
   graphics $mol cylinder $start $middle radius $radius1
   graphics $mol cone $middle $end radius $radius2
}

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
   set path_dat [split $path_filedat "\n"]
   foreach line $path_dat {
      if {[expr {$line eq ""}]} {
         continue
      }
      set splits [split $line]
      graphics $mol sphere [list [lindex $splits 0] [lindex $splits 1] [lindex $splits 2]] radius $radius resolution $resolution
   }
}

proc draw_text {mol text_file {size 1.4} {thickness 2.8}} {
   set fp [open $text_file r]
   set text_filedat [read $fp]
   close $fp
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

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

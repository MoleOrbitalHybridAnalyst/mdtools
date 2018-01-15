echo "times,dx,dy,dz"
box_x=31.075
box_y=31.075
box_z=31.075
perl -ane '$i = @F[1]; $j=@F[2]; $k=@F[3]; if($i>10) {$i-='$box_x';} if($i<-10) {$i+='$box_x';} if($j>10) {$j-='$box_y';} if($j<-10) {$j+='$box_y';} if($k>10) {$k-='$box_z';} if($k<-10) {$k+='$box_z';} print "@F[0],$i,$j,$k\n"' $1

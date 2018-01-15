package require psfgen
set w3id 6267
set w3seg W1
set segP1 [atomselect 0 "segname P1"]
set segO1 [atomselect 0 "segname O1"]
set segO2 [atomselect 0 "segname O2"]
set segO3 [atomselect 0 "segname O3"]
set segO4 [atomselect 0 "segname O4"]
set segW3 [atomselect 0 "segname $w3seg and resid $w3id"]
set segW1 [atomselect 0 "segname W1 and not index [$segW3 list]"]
set segW2 [atomselect 0 "segname W2 and not index [$segW3 list]"]
cd gcv0.000
$segP1 writepdb "chainP1.pdb"
$segO1 writepdb "chainO1.pdb"
$segO2 writepdb "chainO2.pdb"
$segO3 writepdb "chainO3.pdb"
$segO4 writepdb "chainO4.pdb"
$segW3 writepdb "chainW3.pdb"
$segW1 writepdb "chainW1.pdb"
$segW2 writepdb "chainW2.pdb"
cd ..
exit

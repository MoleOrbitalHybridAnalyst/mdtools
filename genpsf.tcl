#in classical 771 and 908 are protonated 
#888 and 876 are disulfide 
package require psfgen
topology /home/chhli/top_zw.rtf

set mutateid 1 

segment P1 {
    pdb chainP1.pdb
    first NTER 
    last CTER
}
coordpdb chainP1.pdb P1

patch GLUP P1:771
#patch GLUP P1:908
patch DISU P1:876 P1:888

regenerate angles dihedrals

segment P2 {
    pdb chainP2.pdb
    first NTER 
    last CTER
}
coordpdb chainP2.pdb P2


segment O1 {
    pdb chainO1.pdb
    first NONE
    last NONE
}
coordpdb chainO1.pdb O1

segment O2 {
    pdb chainO2.pdb
    first NONE
    last NONE
}
coordpdb chainO2.pdb O2

segment O3 {
    pdb chainO3.pdb
    first NONE
    last NONE
}
coordpdb chainO3.pdb O3

segment W1 {
    auto none
    pdb chainW1.pdb
    first NONE
    last NONE
}
coordpdb chainW1.pdb W1

segment W2 {
    auto none
    pdb chainW2.pdb
    first NONE
    last NONE
}
coordpdb chainW2.pdb W2

segment W3 {
    auto none
    pdb chainW3.pdb
    first NONE
    last NONE
}
coordpdb chainW3.pdb W3

segment W4 {
    auto none
    pdb chainW4.pdb
    first NONE
    last NONE
}
coordpdb chainW4.pdb W4

segment W5 {
    auto none
    pdb chainW5.pdb
    first NONE
    last NONE
}
coordpdb chainW5.pdb W5

segment W6 {
    auto none
    pdb chainW6.pdb
    first NONE
    last NONE
}
coordpdb chainW6.pdb W6

segment W7 {
    auto none
    pdb chainW7.pdb
    first NONE
    last NONE
}
coordpdb chainW7.pdb W7

segment W8 {
    auto none
    pdb chainW8.pdb
    first NONE
    last NONE
}
coordpdb chainW8.pdb W8

segment W9 {
    auto none
    pdb chainW9.pdb
    first NONE
    last NONE
}
coordpdb chainW9.pdb W9

segment W10 {
    auto none
    pdb chainW10.pdb
    first NONE
    last NONE
}
coordpdb chainW10.pdb W10

segment W11 {
    auto none
    pdb chainW11.pdb
    first NONE
    last NONE
}
coordpdb chainW11.pdb W11

segment W12 {
    auto none
    pdb chainW12.pdb
    first NONE
    last NONE
}
coordpdb chainW12.pdb W12

segment W13 {
    auto none
    pdb chainW13.pdb
    first NONE
    last NONE
}
coordpdb chainW13.pdb W13

segment W14 {
    auto none
    pdb chainW14.pdb
    first NONE
    last NONE
}
coordpdb chainW14.pdb W14

segment W15 {
    auto none
    pdb chainW15.pdb
    first NONE
    last NONE
}
coordpdb chainW15.pdb W15

segment W16 {
    auto none
    pdb chainW16.pdb
    mutate $mutateid H3O
    first NONE
    last NONE
}
coordpdb chainW16.pdb W16

regenerate angles dihedrals


guesscoord

writepdb check.pdb
writepsf check.psf
exit

mol representation DynamicBonds 1.4 0.1 12.0
mol representation VDW 0.8 100.0
mol material AOShiny
mol showrep 0 0 off
mol addrep 0
mol modselect 1 0 carbon
mol modstyle 1 0 Bonds 0.1 100
mol addrep 0
mol modselect 2 0 (z<-16 or z>16) and (oxygen) and within 8 of carbon and index < 2000
mol modstyle 2 0 VDW
mol selupdate 2 0 1
mol addrep 0
mol modselect 3 0 (x^2+y^2)<16 and (z>-16.8 and z<16.8) and (oxygen or hydrogen) and within 8 of carbon and index < 2000
mol modstyle 3 0 VDW
mol selupdate 3 0 1
mol addrep 0
mol modselect 4 0 hydrogen and within 1.2 of index 3999
mol modstyle 4 0 VDW
mol addrep 0
mol modselect 5 0 index 3999
mol modstyle 5 0 VDW
display projection orthographic
rotate y by 90

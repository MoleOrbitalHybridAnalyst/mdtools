#import sys, commands, re, glob, types
import sys, re, glob, types
from os import popen
from math import *             # any function could be used by set()
import numpy as np

class dump:

  # --------------------------------------------------------------------

  def __init__(self,*list):
    self.snaps = []
    self.nsnaps = self.nselect = 0
    self.names = {}
    self.tselect = tselect(self)
    self.aselect = aselect(self)
    self.atype = "type"
    self.bondflag = 0
    self.bondlist = []
    self.triflag = 0
    self.trilist = []
    self.lineflag = 0
    self.linelist = []
    self.objextra = None

    # flist = list of all dump file names

    words = list[0].split()
    self.flist = []
    for word in words: self.flist += glob.glob(word)
    if len(self.flist) == 0 and len(list) == 1:
      raise StandardError("no dump file specified")
    
    if len(list) == 1:
      self.increment = 0
      self.read_all()
    else:
      self.increment = 1
      self.nextfile = 0
      self.eof = 0

  # --------------------------------------------------------------------

  def read_all(self):

    # read all snapshots from each file
    # test for gzipped files

    for file in self.flist:
      if file[-3:] == ".gz":
        f = popen("%s -c %s" % (gunzip,file),'r')
      else: f = open(file)

      snap = self.read_snapshot(f)
      while snap:
        self.snaps.append(snap)
        print(snap.time, end = " ")
        sys.stdout.flush()
        snap = self.read_snapshot(f)

      f.close()
    print()

    # sort entries by timestep, cull duplicates

    self.snaps.sort(key=self.compare_time)
    self.cull()
    self.nsnaps = len(self.snaps)
    print("read %d snapshots" % self.nsnaps)

    # select all timesteps and atoms

    self.tselect.all()

    # print column assignments
    
    if len(self.names):
      print("assigned columns:",self.names2str())
    else:
      print("no column assignments made")

    # if snapshots are scaled, unscale them

    if (not self.names.has_key("x")) or \
       (not self.names.has_key("y")) or \
       (not self.names.has_key("z")):
      print("dump scaling status is unknown")
    elif self.nsnaps > 0:
      if self.scale_original == 1: 
          self.unscale()
      elif self.scale_original == 0: 
          print("dump is already unscaled")
      else: print("dump scaling status is unknown")

  # --------------------------------------------------------------------
  # read next snapshot from list of files

  def next(self):

    if not self.increment: 
        raise StandardError("cannot read incrementally")

    # read next snapshot in current file using eof as pointer
    # if fail, try next file
    # if new snapshot time stamp already exists, read next snapshot

    while 1:
      f = open(self.flist[self.nextfile],'rb')
      f.seek(self.eof)
      snap = self.read_snapshot(f)
      if not snap:
        self.nextfile += 1
        if self.nextfile == len(self.flist):
            return -1
        f.close()
        self.eof = 0
        continue
      self.eof = f.tell()
      f.close()
      try:
        self.findtime(snap.time)
        continue
      except: 
        break

    # select the new snapshot with all its atoms

    self.snaps.append(snap)
    snap = self.snaps[self.nsnaps]
    snap.tselect = 1
    snap.nselect = snap.natoms
    for i in range(snap.natoms): snap.aselect[i] = 1
    self.nsnaps += 1
    self.nselect += 1

    return snap.time

  # --------------------------------------------------------------------
  # read a single snapshot from file f
  # return snapshot or 0 if failed
  # for first snapshot only:
  #   assign column names (file must be self-describing)
  #   set scale_original to 0/1/-1 for unscaled/scaled/unknown
  #   convert xs,xu to x in names
  
  def read_snapshot(self,f):
    try:
      snap = Snap()
      item = f.readline()
      snap.time = int(f.readline().split()[0])    # just grab 1st field
      item = f.readline()
      snap.natoms = int(f.readline())

      snap.aselect = np.zeros(snap.natoms)

      item = f.readline()
      words = item.split("BOUNDS ")
      if len(words) == 1: snap.boxstr = ""
      else: snap.boxstr = words[1].strip()
      if "xy" in snap.boxstr: snap.triclinic = 1
      else: snap.triclinic = 0
      
      words = f.readline().split()
      if len(words) == 2:
        snap.xlo,snap.xhi,snap.xy = float(words[0]),float(words[1]),0.0
      else:
        snap.xlo,snap.xhi,snap.xy = \
            float(words[0]),float(words[1]),float(words[2])

      words = f.readline().split()
      if len(words) == 2:
        snap.ylo,snap.yhi,snap.xz = float(words[0]),float(words[1]),0.0
      else:
        snap.ylo,snap.yhi,snap.xz = \
            float(words[0]),float(words[1]),float(words[2])

      words = f.readline().split()
      if len(words) == 2:
        snap.zlo,snap.zhi,snap.yz = float(words[0]),float(words[1]),0.0
      else:
        snap.zlo,snap.zhi,snap.yz = \
            float(words[0]),float(words[1]),float(words[2])
          
      item = f.readline()
      if len(self.names) == 0:
        self.scale_original = -1
        xflag = yflag = zflag = -1
        words = item.split()[2:]
        if len(words):
          for i in range(len(words)):
            if words[i] == "x" or words[i] == "xu":
              xflag = 0
              self.names["x"] = i
            elif words[i] == "xs" or words[i] == "xsu":
              xflag = 1
              self.names["x"] = i
            elif words[i] == "y" or words[i] == "yu":
              yflag = 0
              self.names["y"] = i
            elif words[i] == "ys" or words[i] == "ysu":
              yflag = 1
              self.names["y"] = i
            elif words[i] == "z" or words[i] == "zu":
              zflag = 0
              self.names["z"] = i
            elif words[i] == "zs" or words[i] == "zsu":
              zflag = 1
              self.names["z"] = i
            else: self.names[words[i]] = i
          if xflag == 0 and yflag == 0 and zflag == 0: self.scale_original = 0
          if xflag == 1 and yflag == 1 and zflag == 1: self.scale_original = 1
          
      if snap.natoms:
        words = f.readline().split()
        ncol = len(words)
        for i in range(1,snap.natoms):
          words += f.readline().split()
        floats = map(float,words)
        if oldnumeric: atoms = np.zeros((snap.natoms,ncol),np.Float)
        else: atoms = np.zeros((snap.natoms,ncol),np.float)
        start = 0
        stop = ncol
        for i in range(snap.natoms):
          atoms[i] = floats[start:stop]
          start = stop
          stop += ncol
      else: atoms = None
      snap.atoms = atoms
      return snap
    except:
      return 0

  # --------------------------------------------------------------------
  # map atom column names
  
  def map(self,*pairs):
    if len(pairs) % 2 != 0:
      raise StandardError("dump map() requires pairs of mappings")
    for i in range(0,len(pairs),2):
      j = i + 1
      self.names[pairs[j]] = pairs[i]-1

  # --------------------------------------------------------------------
  # delete unselected snapshots

  def delete(self):
    ndel = i = 0
    while i < self.nsnaps:
      if not self.snaps[i].tselect:
        del self.snaps[i]
        self.nsnaps -= 1
        ndel += 1
      else: i += 1
    print("%d snapshots deleted" % ndel)
    print("%d snapshots remaining" % self.nsnaps)

  # --------------------------------------------------------------------
  # scale coords to 0-1 for all snapshots or just one
  # use 6 params as h-matrix to treat orthongonal or triclinic boxes

  def scale(self,*list):
    if len(list) == 0:
      print("Scaling dump ...")
      x = self.names["x"]
      y = self.names["y"]
      z = self.names["z"]
      for snap in self.snaps: self.scale_one(snap,x,y,z)
    else:
      i = self.findtime(list[0])
      x = self.names["x"]
      y = self.names["y"]
      z = self.names["z"]
      self.scale_one(self.snaps[i],x,y,z)

  # --------------------------------------------------------------------

  def scale_one(self,snap,x,y,z):
    if snap.xy == 0.0 and snap.xz == 0.0 and snap.yz == 0.0:
      xprdinv = 1.0 / (snap.xhi - snap.xlo)
      yprdinv = 1.0 / (snap.yhi - snap.ylo)
      zprdinv = 1.0 / (snap.zhi - snap.zlo)
      atoms = snap.atoms
      if atoms != None:
        atoms[:,x] = (atoms[:,x] - snap.xlo) * xprdinv
        atoms[:,y] = (atoms[:,y] - snap.ylo) * yprdinv
        atoms[:,z] = (atoms[:,z] - snap.zlo) * zprdinv
    else:
      xlo_bound = snap.xlo; xhi_bound = snap.xhi
      ylo_bound = snap.ylo; yhi_bound = snap.yhi
      zlo_bound = snap.zlo; zhi_bound = snap.zhi
      xy = snap.xy
      xz = snap.xz
      yz = snap.yz
      xlo = xlo_bound - min((0.0,xy,xz,xy+xz))
      xhi = xhi_bound - max((0.0,xy,xz,xy+xz))
      ylo = ylo_bound - min((0.0,yz))
      yhi = yhi_bound - max((0.0,yz))
      zlo = zlo_bound
      zhi = zhi_bound
      h0 = xhi - xlo
      h1 = yhi - ylo
      h2 = zhi - zlo
      h3 = yz
      h4 = xz
      h5 = xy
      h0inv = 1.0 / h0
      h1inv = 1.0 / h1
      h2inv = 1.0 / h2
      h3inv = yz / (h1*h2)
      h4inv = (h3*h5 - h1*h4) / (h0*h1*h2)
      h5inv = xy / (h0*h1)
      atoms = snap.atoms
      if atoms != None:
        atoms[:,x] = (atoms[:,x] - snap.xlo)*h0inv + \
            (atoms[:,y] - snap.ylo)*h5inv + \
            (atoms[:,z] - snap.zlo)*h4inv
        atoms[:,y] = (atoms[:,y] - snap.ylo)*h1inv + \
            (atoms[:,z] - snap.zlo)*h3inv
        atoms[:,z] = (atoms[:,z] - snap.zlo)*h2inv
        
  # --------------------------------------------------------------------
  # unscale coords from 0-1 to box size for all snapshots or just one
  # use 6 params as h-matrix to treat orthongonal or triclinic boxes

  def unscale(self,*list):
    if len(list) == 0:
      print("Unscaling dump ...")
      x = self.names["x"]
      y = self.names["y"]
      z = self.names["z"]
      for snap in self.snaps: self.unscale_one(snap,x,y,z)
    else:
      i = self.findtime(list[0])
      x = self.names["x"]
      y = self.names["y"]
      z = self.names["z"]
      self.unscale_one(self.snaps[i],x,y,z)

  # --------------------------------------------------------------------

  def unscale_one(self,snap,x,y,z):
    if snap.xy == 0.0 and snap.xz == 0.0 and snap.yz == 0.0:
      xprd = snap.xhi - snap.xlo
      yprd = snap.yhi - snap.ylo
      zprd = snap.zhi - snap.zlo
      atoms = snap.atoms
      if atoms != None:
        atoms[:,x] = snap.xlo + atoms[:,x]*xprd
        atoms[:,y] = snap.ylo + atoms[:,y]*yprd
        atoms[:,z] = snap.zlo + atoms[:,z]*zprd
    else:
      xlo_bound = snap.xlo; xhi_bound = snap.xhi
      ylo_bound = snap.ylo; yhi_bound = snap.yhi
      zlo_bound = snap.zlo; zhi_bound = snap.zhi
      xy = snap.xy
      xz = snap.xz
      yz = snap.yz
      xlo = xlo_bound - min((0.0,xy,xz,xy+xz))
      xhi = xhi_bound - max((0.0,xy,xz,xy+xz))
      ylo = ylo_bound - min((0.0,yz))
      yhi = yhi_bound - max((0.0,yz))
      zlo = zlo_bound
      zhi = zhi_bound
      h0 = xhi - xlo
      h1 = yhi - ylo
      h2 = zhi - zlo
      h3 = yz
      h4 = xz
      h5 = xy
      atoms = snap.atoms
      if atoms != None:
        atoms[:,x] = snap.xlo + atoms[:,x]*h0 + atoms[:,y]*h5 + atoms[:,z]*h4
        atoms[:,y] = snap.ylo + atoms[:,y]*h1 + atoms[:,z]*h3
        atoms[:,z] = snap.zlo + atoms[:,z]*h2
        
  # --------------------------------------------------------------------
  # wrap coords from outside box to inside

  def wrap(self):
    print("Wrapping dump ...")

    x = self.names["x"]
    y = self.names["y"]
    z = self.names["z"]
    ix = self.names["ix"]
    iy = self.names["iy"]
    iz = self.names["iz"]
    
    for snap in self.snaps:
      xprd = snap.xhi - snap.xlo
      yprd = snap.yhi - snap.ylo
      zprd = snap.zhi - snap.zlo
      atoms = snap.atoms
      atoms[:,x] -= atoms[:,ix]*xprd
      atoms[:,y] -= atoms[:,iy]*yprd
      atoms[:,z] -= atoms[:,iz]*zprd

  # --------------------------------------------------------------------
  # unwrap coords from inside box to outside

  def unwrap(self):
    print("Unwrapping dump ...")

    x = self.names["x"]
    y = self.names["y"]
    z = self.names["z"]
    ix = self.names["ix"]
    iy = self.names["iy"]
    iz = self.names["iz"]
    
    for snap in self.snaps:
      xprd = snap.xhi - snap.xlo
      yprd = snap.yhi - snap.ylo
      zprd = snap.zhi - snap.zlo
      atoms = snap.atoms
      atoms[:,x] += atoms[:,ix]*xprd
      atoms[:,y] += atoms[:,iy]*yprd
      atoms[:,z] += atoms[:,iz]*zprd

  # --------------------------------------------------------------------
  # wrap coords to same image as atom ID stored in "other" column
  # if dynamic extra lines or triangles defined, owrap them as well
      
  def owrap(self,other):
    print("Wrapping to other ...")
    
    id = self.names["id"]
    x = self.names["x"]
    y = self.names["y"]
    z = self.names["z"]
    ix = self.names["ix"]
    iy = self.names["iy"]
    iz = self.names["iz"]
    iother = self.names[other]

    for snap in self.snaps:
      xprd = snap.xhi - snap.xlo
      yprd = snap.yhi - snap.ylo
      zprd = snap.zhi - snap.zlo
      atoms = snap.atoms
      ids = {}
      for i in range(snap.natoms):
        ids[atoms[i][id]] = i
      for i in range(snap.natoms):
        j = ids[atoms[i][iother]]
        atoms[i][x] += (atoms[i][ix]-atoms[j][ix])*xprd
        atoms[i][y] += (atoms[i][iy]-atoms[j][iy])*yprd
        atoms[i][z] += (atoms[i][iz]-atoms[j][iz])*zprd
      # should bonds also be owrapped ?
      if self.lineflag == 2 or self.triflag == 2:
        self.objextra.owrap(snap.time,xprd,yprd,zprd,ids,atoms,iother,ix,iy,iz)
          
  # --------------------------------------------------------------------
  # convert column names assignment to a string, in column order
  
  def names2str(self):
    pairs = list(self.names.items())
    values = list(self.names.values())
    ncol = len(pairs)
    str = ""
    for i in range(ncol):
      if i in values: str += pairs[values.index(i)][0] + ' '
    return str

  # --------------------------------------------------------------------
  # sort atoms by atom ID in all selected timesteps by default
  # if arg = string, sort all steps by that column
  # if arg = numeric, sort atoms in single step

  def sort(self,*list):
    if len(list) == 0:
      print("Sorting selected snapshots ...")
      id = self.names["id"]
      for snap in self.snaps:
        if snap.tselect: self.sort_one(snap,id)
    elif type(list[0]) is types.StringType:
      print("Sorting selected snapshots by %s ..." % list[0])
      id = self.names[list[0]]
      for snap in self.snaps:
        if snap.tselect: self.sort_one(snap,id)
    else:
      i = self.findtime(list[0])
      id = self.names["id"]
      self.sort_one(self.snaps[i],id)

  # --------------------------------------------------------------------
  # sort a single snapshot by ID column

  def sort_one(self,snap,id):
    atoms = snap.atoms
    ids = atoms[:,id]
    ordering = np.argsort(ids)
    for i in range(len(atoms[0])):
      atoms[:,i] = np.take(atoms[:,i],ordering)

  # --------------------------------------------------------------------
  # write a single dump file from current selection

  def write(self,file,header=1,append=0):
    if len(self.snaps): 
        namestr = self.names2str()
    if not append: 
        f = open(file,"w")
    else: 
        f = open(file,"a")

    if "id" in self.names: id = self.names["id"]
    else: id = -1
    if "type" in self.names: type = self.names["type"]
    else: type = -1
    
    for snap in self.snaps:
      if not snap.tselect: continue
      print(snap.time, end = " ")
      sys.stdout.flush()

      if header:
        print("ITEM: TIMESTEP", file=f)
        print(snap.time, file=f)
        print("ITEM: NUMBER OF ATOMS", file=f)
        print(snap.nselect, file=f)
        if snap.boxstr: 
            print("ITEM: BOX BOUNDS",snap.boxstr, file=f)
        else: print("ITEM: BOX BOUNDS", file=f)
        if snap.triclinic:
          print(snap.xlo,snap.xhi,snap.xy, file=f) 
          print(snap.ylo,snap.yhi,snap.xz, file=f)
          print(snap.zlo,snap.zhi,snap.yz, file=f)
        else:
          print(snap.xlo,snap.xhi, file=f) 
          print(snap.ylo,snap.yhi, file=f)
          print(snap.zlo,snap.zhi, file=f)
        print("ITEM: ATOMS",namestr, file=f)
      
      atoms = snap.atoms
      nvalues = len(atoms[0])
      for i in range(snap.natoms):
        if not snap.aselect[i]: continue
        line = ""
        for j in range(nvalues):
          if j == id or j == type:
            line += str(int(atoms[i][j])) + " "
          else:
            line += str(atoms[i][j]) + " "
        print(line, file=f)
    f.close()
    print("\n%d snapshots" % self.nselect)

  # --------------------------------------------------------------------
  # write one dump file per snapshot from current selection

  def scatter(self,root):
    if len(self.snaps): namestr = self.names2str()
    for snap in self.snaps:
      if not snap.tselect: continue
      print(snap.time, end = " ")
      sys.stdout.flush()
      
      file = root + "." + str(snap.time)
      f = open(file,"w")
      print("ITEM: TIMESTEP", file=f)
      print(snap.time, file=f)
      print("ITEM: NUMBER OF ATOMS", file=f)
      print(snap.nselect, file=f)
      if snap.boxstr: 
          print("ITEM: BOX BOUNDS",snap.boxstr, file=f)
      else: 
          print("ITEM: BOX BOUNDS", file=f)
      if snap.triclinic:
        print(snap.xlo,snap.xhi,snap.xy, file=f) 
        print(snap.ylo,snap.yhi,snap.xz, file=f)
        print(snap.zlo,snap.zhi,snap.yz, file=f)
      else:
        print(snap.xlo,snap.xhi, file=f) 
        print(snap.ylo,snap.yhi, file=f)
        print(snap.zlo,snap.zhi, file=f)
      print("ITEM: ATOMS",namestr, file=f)
      
      atoms = snap.atoms
      nvalues = len(atoms[0])
      for i in range(snap.natoms):
        if not snap.aselect[i]: continue
        line = ""
        for j in range(nvalues):
          if (j < 2):
            line += str(int(atoms[i][j])) + " "
          else:
            line += str(atoms[i][j]) + " "
        print(line, file=f)
      f.close()
    print("\n%d snapshots" % self.nselect)

  # --------------------------------------------------------------------
  # find min/max across all selected snapshots/atoms for a particular column

  def minmax(self,colname):
    icol = self.names[colname]
    min = 1.0e20
    max = -min
    for snap in self.snaps:
      if not snap.tselect: continue
      atoms = snap.atoms
      for i in range(snap.natoms):
        if not snap.aselect[i]: continue
        if atoms[i][icol] < min: min = atoms[i][icol]
        if atoms[i][icol] > max: max = atoms[i][icol]
    return (min,max)

  # --------------------------------------------------------------------
  # set a column value via an equation for all selected snapshots

  def set(self,eq):
    print("Setting ...")
    pattern = "\$\w*"
    list = re.findall(pattern,eq)

    lhs = list[0][1:]
    if not self.names.has_key(lhs):
      self.newcolumn(lhs)
      
    for item in list:
      name = item[1:]
      column = self.names[name]
      insert = "snap.atoms[i][%d]" % (column)
      eq = eq.replace(item,insert)
    ceq = compile(eq,'','single')

    for snap in self.snaps:
      if not snap.tselect: continue
      for i in range(snap.natoms):
        if snap.aselect[i]: 
            exec(ceq)
          
  # --------------------------------------------------------------------
  # set a column value via an input vec for all selected snapshots/atoms

  def setv(self,colname,vec):
    print("Setting ...")
    if not self.names.has_key(colname):
      self.newcolumn(colname)
    icol = self.names[colname]

    for snap in self.snaps:
      if not snap.tselect: continue
      if snap.nselect != len(vec):
        raise StandardError("vec length does not match # of selected atoms")
      atoms = snap.atoms
      m = 0
      for i in range(snap.natoms):
        if snap.aselect[i]:
          atoms[i][icol] = vec[m]
          m += 1
          
  # --------------------------------------------------------------------
  # clone value in col across selected timesteps for atoms with same ID

  def clone(self,nstep,col):
    istep = self.findtime(nstep)
    icol = self.names[col]
    id = self.names["id"]
    ids = {}
    for i in range(self.snaps[istep].natoms):
      ids[self.snaps[istep].atoms[i][id]] = i
    for snap in self.snaps:
      if not snap.tselect: continue
      atoms = snap.atoms
      for i in range(snap.natoms):
        if not snap.aselect[i]: continue
        j = ids[atoms[i][id]]
        atoms[i][icol] = self.snaps[istep].atoms[j][icol]

  # --------------------------------------------------------------------
  # values in old column are spread as ints from 1-N and assigned to new column

  def spread(self,old,n,new):
    iold = self.names[old]
    if not self.names.has_key(new): self.newcolumn(new)
    inew = self.names[new]

    min,max = self.minmax(old)
    print("min/max = ",min,max)

    gap = max - min
    invdelta = n/gap
    for snap in self.snaps:
      if not snap.tselect: continue
      atoms = snap.atoms
      for i in range(snap.natoms):
        if not snap.aselect[i]: continue
        ivalue = int((atoms[i][iold] - min) * invdelta) + 1
        if ivalue > n: ivalue = n
        if ivalue < 1: ivalue = 1
        atoms[i][inew] = ivalue

  # --------------------------------------------------------------------
  # return vector of selected snapshot time stamps

  def time(self):
    vec = self.nselect * [0]
    i = 0
    for snap in self.snaps:
      if not snap.tselect: continue
      vec[i] = snap.time
      i += 1
    return vec

  # --------------------------------------------------------------------
  # extract vector(s) of values for atom ID n at each selected timestep

  def atom(self,n,*list):
    if len(list) == 0:
      raise StandardError("no columns specified")
    columns = []
    values = []
    for name in list:
      columns.append(self.names[name])
      values.append(self.nselect * [0])
    ncol = len(columns)
    
    id = self.names["id"]
    m = 0
    for snap in self.snaps:
      if not snap.tselect: continue
      atoms = snap.atoms
      for i in range(snap.natoms):
        if atoms[i][id] == n: break
      if atoms[i][id] != n:
        raise StandardError("could not find atom ID in snapshot")
      for j in range(ncol):
        values[j][m] = atoms[i][columns[j]]
      m += 1

    if len(list) == 1: return values[0]
    else: return values
  
  # --------------------------------------------------------------------
  # extract vector(s) of values for selected atoms at chosen timestep

  def vecs(self,n,*list):
    snap = self.snaps[self.findtime(n)]
    
    if len(list) == 0:
      raise StandardError("no columns specified")
    columns = []
    values = []
    for name in list:
      columns.append(self.names[name])
      values.append(snap.nselect * [0])
    ncol = len(columns)

    m = 0
    for i in range(snap.natoms):
      if not snap.aselect[i]: continue
      for j in range(ncol):
        values[j][m] = snap.atoms[i][columns[j]]
      m += 1

    if len(list) == 1: return values[0]
    else: return values

  # --------------------------------------------------------------------
  # add a new column to every snapshot and set value to 0
  # set the name of the column to str

  def newcolumn(self,str):
    ncol = len(self.snaps[0].atoms[0])
    self.map(ncol+1,str)
    for snap in self.snaps:
      atoms = snap.atoms
      if oldnumeric: newatoms = np.zeros((snap.natoms,ncol+1),np.Float)
      else: newatoms = np.zeros((snap.natoms,ncol+1),np.float)
      newatoms[:,0:ncol] = snap.atoms
      snap.atoms = newatoms

  # --------------------------------------------------------------------
  # sort snapshots on time stamp

  def compare_time(self,a,b):
    if a.time < b.time:
      return -1
    elif a.time > b.time:
      return 1
    else:
      return 0

  # --------------------------------------------------------------------
  # delete successive snapshots with duplicate time stamp

  def cull(self):
    i = 1
    while i < len(self.snaps):
      if self.snaps[i].time == self.snaps[i-1].time:
        del self.snaps[i]
      else:
        i += 1
  
  # --------------------------------------------------------------------
  # iterate over selected snapshots

  def iterator(self,flag):
    start = 0
    if flag: start = self.iterate + 1
    for i in range(start,self.nsnaps):
      if self.snaps[i].tselect:
        self.iterate = i
        return i,self.snaps[i].time,1
    return 0,0,-1
  
  # --------------------------------------------------------------------
  # return list of atoms to viz for snapshot isnap
  # if called with flag, then index is timestep, so convert to snapshot index
  # augment with bonds, tris, lines if extra() was invoked
  
  def viz(self,index,flag=0):
    if not flag: isnap = index
    else:
      times = self.time()
      n = len(times)
      i = 0
      while i < n:
        if times[i] > index: break
        i += 1
      isnap = i - 1

    snap = self.snaps[isnap]

    time = snap.time
    box = [snap.xlo,snap.ylo,snap.zlo,snap.xhi,snap.yhi,snap.zhi]
    id = self.names["id"]
    type = self.names[self.atype]
    x = self.names["x"]
    y = self.names["y"]
    z = self.names["z"]

    # create atom list needed by viz from id,type,x,y,z
    # need Numeric/Numpy mode here
    
    atoms = []
    for i in range(snap.natoms):
      if not snap.aselect[i]: continue
      atom = snap.atoms[i]
      atoms.append([atom[id],atom[type],atom[x],atom[y],atom[z]])

    # create list of bonds from static or dynamic bond list
    # then generate bond coords from bondlist
    # alist = dictionary of atom IDs for atoms list
    # lookup bond atom IDs in alist and grab their coords
    # try is used since some atoms may be unselected
    #   any bond with unselected atom is not added to bonds
    # need Numeric/Numpy mode here

    bonds = []
    if self.bondflag:
      if self.bondflag == 1: bondlist = self.bondlist
      elif self.bondflag == 2:
        tmp1,tmp2,tmp3,bondlist,tmp4,tmp5 = self.objextra.viz(time,1)
      alist = {}
      for i in range(len(atoms)): alist[int(atoms[i][0])] = i
      for bond in bondlist:
        try:
          i = alist[bond[2]]
          j = alist[bond[3]]
          atom1 = atoms[i]
          atom2 = atoms[j]
          bonds.append([bond[0],bond[1],atom1[2],atom1[3],atom1[4],
                        atom2[2],atom2[3],atom2[4],atom1[1],atom2[1]])
        except: continue

    # create list of tris from static or dynamic tri list
    # if dynamic, could eliminate tris for unselected atoms

    tris = []
    if self.triflag:
      if self.triflag == 1: tris = self.trilist
      elif self.triflag == 2:
        tmp1,tmp2,tmp3,tmp4,tris,tmp5 = self.objextra.viz(time,1)
      
    # create list of lines from static or dynamic tri list
    # if dynamic, could eliminate lines for unselected atoms

    lines = []
    if self.lineflag:
      if self.lineflag == 1: lines = self.linelist
      elif self.lineflag == 2:
        tmp1,tmp2,tmp3,tmp4,tmp5,lines = self.objextra.viz(time,1)

    return time,box,atoms,bonds,tris,lines
  
  # --------------------------------------------------------------------

  def findtime(self,n):
    for i in range(self.nsnaps):
      if self.snaps[i].time == n: return i
    raise StandardError("no step %d exists" % n)

  # --------------------------------------------------------------------
  # return maximum box size across all selected snapshots

  def maxbox(self):
    xlo = ylo = zlo = None
    xhi = yhi = zhi = None
    for snap in self.snaps:
      if not snap.tselect: continue
      if xlo == None or snap.xlo < xlo: xlo = snap.xlo
      if xhi == None or snap.xhi > xhi: xhi = snap.xhi
      if ylo == None or snap.ylo < ylo: ylo = snap.ylo
      if yhi == None or snap.yhi > yhi: yhi = snap.yhi
      if zlo == None or snap.zlo < zlo: zlo = snap.zlo
      if zhi == None or snap.zhi > zhi: zhi = snap.zhi
    return [xlo,ylo,zlo,xhi,yhi,zhi]

  # --------------------------------------------------------------------
  # return maximum atom type across all selected snapshots and atoms

  def maxtype(self):
    icol = self.names["type"]
    max = 0
    for snap in self.snaps:
      if not snap.tselect: continue
      atoms = snap.atoms
      for i in range(snap.natoms):
        if not snap.aselect[i]: continue
        if atoms[i][icol] > max: max = atoms[i][icol]
    return int(max)

  # --------------------------------------------------------------------
  # grab bonds/tris/lines from another object
  # if static, grab once, else store obj to grab dynamically

  def extra(self,arg):

    # data object, grab bonds statically
    
    if type(arg) is types.InstanceType and ".data" in str(arg.__class__):
      self.bondflag = 0
      try:
        bondlist = []
        bondlines = arg.sections["Bonds"]
        for line in bondlines:
          words = line.split()
          bondlist.append([int(words[0]),int(words[1]),
                           int(words[2]),int(words[3])])
        if bondlist:
          self.bondflag = 1
          self.bondlist = bondlist
      except:
        raise StandardError("could not extract bonds from data object")

    # cdata object, grab tris and lines statically
    
    elif type(arg) is types.InstanceType and ".cdata" in str(arg.__class__):
      self.triflag = self.lineflag = 0
      try:
        tmp,tmp,tmp,tmp,tris,lines = arg.viz(0)
        if tris:
          self.triflag = 1
          self.trilist = tris
        if lines:
          self.lineflag = 1
          self.linelist = lines
      except:
        raise StandardError("could not extract tris/lines from cdata object")

    # mdump object, grab tris dynamically
    
    elif type(arg) is types.InstanceType and ".mdump" in str(arg.__class__):
      self.triflag = 2
      self.objextra = arg

    # bdump object, grab bonds dynamically
    
    elif type(arg) is types.InstanceType and ".bdump" in str(arg.__class__):
      self.bondflag = 2
      self.objextra = arg

    # ldump object, grab lines dynamically
    
    elif type(arg) is types.InstanceType and ".ldump" in str(arg.__class__):
      self.lineflag = 2
      self.objextra = arg

    # tdump object, grab tris dynamically
    
    elif type(arg) is types.InstanceType and ".tdump" in str(arg.__class__):
      self.triflag = 2
      self.objextra = arg

    else:
      raise StandardError("unrecognized argument to dump.extra()")
      
  # --------------------------------------------------------------------

  def compare_atom(self,a,b):
    if a[0] < b[0]:
      return -1
    elif a[0] > b[0]:
      return 1
    else:
      return 0  

# --------------------------------------------------------------------
# one snapshot

class Snap:
  pass

# --------------------------------------------------------------------
# time selection class

class tselect:

  def __init__(self,data):
    self.data = data
    
  # --------------------------------------------------------------------

  def all(self):
    data = self.data
    for snap in data.snaps:
      snap.tselect = 1
    data.nselect = len(data.snaps)
    data.aselect.all()
    print("%d snapshots selected out of %d" % (data.nselect,data.nsnaps))

  # --------------------------------------------------------------------

  def one(self,n):
    data = self.data
    for snap in data.snaps:
      snap.tselect = 0
    i = data.findtime(n)
    data.snaps[i].tselect = 1
    data.nselect = 1
    data.aselect.all()
    print("%d snapshots selected out of %d" % (data.nselect,data.nsnaps))

  # --------------------------------------------------------------------

  def none(self):
    data = self.data
    for snap in data.snaps:
      snap.tselect = 0
    data.nselect = 0
    print("%d snapshots selected out of %d" % (data.nselect,data.nsnaps))

  # --------------------------------------------------------------------

  def skip(self,n):
    data = self.data
    count = n-1
    for snap in data.snaps:
      if not snap.tselect: continue
      count += 1
      if count == n:
        count = 0
        continue
      snap.tselect = 0
      data.nselect -= 1
    data.aselect.all()
    print("%d snapshots selected out of %d" % (data.nselect,data.nsnaps))
  
  # --------------------------------------------------------------------

  def test(self,teststr):
    data = self.data
    snaps = data.snaps
    cmd = "flag = " + teststr.replace("$t","snaps[i].time")
    ccmd = compile(cmd,'','single')
    for i in range(data.nsnaps):
      if not snaps[i].tselect: continue
      exec(ccmd)
      if not flag:
        snaps[i].tselect = 0
        data.nselect -= 1
    data.aselect.all()
    print("%d snapshots selected out of %d" % (data.nselect,data.nsnaps))

# --------------------------------------------------------------------
# atom selection class

class aselect:

  def __init__(self,data):
    self.data = data

  # --------------------------------------------------------------------

  def all(self,*args):
    data = self.data
    if len(args) == 0:                           # all selected timesteps
      for snap in data.snaps:
        if not snap.tselect: continue
        for i in range(snap.natoms): snap.aselect[i] = 1
        snap.nselect = snap.natoms
    else:                                        # one timestep
      n = data.findtime(args[0])
      snap = data.snaps[n]
      for i in range(snap.natoms): snap.aselect[i] = 1
      snap.nselect = snap.natoms

  # --------------------------------------------------------------------

  def test(self,teststr,*args):
    data = self.data

    # replace all $var with snap.atoms references and compile test string
    
    pattern = "\$\w*"
    list = re.findall(pattern,teststr)
    for item in list:
      name = item[1:]
      column = data.names[name]
      insert = "snap.atoms[i][%d]" % column
      teststr = teststr.replace(item,insert)
    cmd = "flag = " + teststr
    ccmd = compile(cmd,'','single')

    if len(args) == 0:                           # all selected timesteps
      for snap in data.snaps:
        if not snap.tselect: continue
        for i in range(snap.natoms):
          if not snap.aselect[i]: continue
          exec(ccmd)
          if not flag:
            snap.aselect[i] = 0
            snap.nselect -= 1
      for i in range(data.nsnaps):
        if data.snaps[i].tselect:
          print("%d atoms of %d selected in first step %d" % \
                (data.snaps[i].nselect,data.snaps[i].natoms,data.snaps[i].time))
          break
      for i in range(data.nsnaps-1,-1,-1):
        if data.snaps[i].tselect:
          print("%d atoms of %d selected in last step %d" % \
                (data.snaps[i].nselect,data.snaps[i].natoms,data.snaps[i].time))
          break

    else:                                        # one timestep
      n = data.findtime(args[0])
      snap = data.snaps[n]
      for i in range(snap.natoms):
        if not snap.aselect[i]: continue
        exec(ccmd)
        if not flag:
          snap.aselect[i] = 0
          snap.nselect -= 1

#!/usr/local/bin/python3.5

EVB_PAR_NAME = 'evb.par'
EVB_TYPE_NAME = 'evb.type'
EVB_CFG_NAME = 'evb.cfg'
outname = 'evb.check'


defs = set()
def_stack = []
defined = set()
alias_def = dict()

def assert_num(lst, n):
	if len(lst) < n:
		raise RuntimeError("list length should be {:d}, but is {:d}".format(\
				n, len(lst)))
	elif len(lst) > n and lst[n][0] != ':':
		raise RuntimeError("list length should be {:d}, but is {:d}".format(\
				n, len(lst)))

fcfg = open(str(EVB_CFG_NAME), 'r')
for l in fcfg:
	lspl = l.split()
	if len(lspl) == 0:
		continue
	if lspl[0][0] == ':':
		continue
	if lspl[0] == '#define':
		assert_num(lspl, 2)
		defined.add(lspl[1])
fcfg.close()

ftype = open(str(EVB_TYPE_NAME), 'r')
for l in ftype:
	lspl = l.split()
	if len(lspl) == 0:
		continue
	if lspl[0][0] == ':':
		continue
	if lspl[0] == '#define':
		assert_num(lspl, 3)
		alias_def[lspl[1]] = lspl[2]
ftype.close()

fpar = open(str(EVB_PAR_NAME), 'r')
fout = open(str(outname), 'w')
tagStack = [True]
segname = []
# TODO: check the defined names
# TODO: check the numbers
# TODO: check if any missing atom, bond, angle, dih, and imp
definedNames = ['Hydronium', 'DA_Gaussian', 'VII']
for l in fpar:
	lspl = l.split()
	if len(lspl) == 0:
		fout.write('\n')
		continue
	if lspl[0][0] == ':':
		#fout.write(l)
		continue
	if lspl[0] == '#ifdef':
		assert_num(lspl, 2)
		def_stack.append(lspl[1])
		tagStack.append(tagStack[-1] and (lspl[1] in defined))
	elif lspl[0] == '#elif' or lspl[0] == '#elseif':
		assert_num(lspl, 2)
		front = def_stack.pop()
		tagStack.pop()
		defs.add(front)
		def_stack.append(lspl[1])
		tagStack.append(tagStack[-1] and (lspl[1] in defined))
	elif lspl[0] == '#endif':
		assert_num(lspl, 1)
		front = def_stack.pop()
		tagStack.pop()
		defs.add(front)
	elif tagStack[-1]:
		#lwrite = []
		for term in lspl:
			if term[0] == ':':
				break
			if term in alias_def:
				fout.write("{:>8s} ".format(alias_def[term]))
				#lwrite.append(alias_def[term])
			else:
				fout.write("{:>8s} ".format(term))
				#lwrite.append(term)
		fout.write('\n')



assert(len(def_stack) == 0)
print(defs)
print(defined)
print(alias_def)

fpar.close()
fout.close()

# read state.dat eigen_vector.dat
# ad hoc:
#   serial = (residue - 2) * 3 + 29
#   glu is always GLU 
#   trj is printed every STRIDE steps
#   output to ci_evb.dat

from collections import defaultdict

STRIDE = 2000
GLU = 19

if __name__=="__main__":
    raw_eigs=[]
    with open("eigen_vector.dat","r") as fp:
        for line in fp:
            splits = line.split()
            if(int(splits[0])%STRIDE!=0):
                continue
            splits = splits[1:]
            raw_eig = [float(e) for e in splits]
            raw_eigs.append(raw_eig)
    raw_states=[]
    with open("state.dat","r") as fp:
        for line in fp:
            splits = line.split()
            if(int(splits[0])%STRIDE!=0):
                continue
            splits = splits[1:]
            raw_state = [int(e) for e in splits]
            raw_states.append(raw_state)
    assert len(raw_eigs) == len(raw_states)
    fp = open("ci_evb.dat","w")
    for (index,(s,e)) in enumerate(zip(raw_states,raw_eigs)):
        result = defaultdict(float)
        if(len(s)!=len(e)):
            print("incompatible eig and state at step %d"%(STRIDE*index))
        for (si,ei) in zip(s,e):
            # translate residue to heavy atom serial:
            if(si == 1):
                serial = GLU
            else:
                serial = (si - 2) * 3 + 29
            # accumulate the eigs with the same heavy atom:
            result[serial] += ei
            #result[si] += ei
        # print result:
        print("%d "%(STRIDE*index),file=fp,end="")
        for _ in result:
            print("%d %f "%(_,result[_]),end="",file=fp)
        print("",file=fp)

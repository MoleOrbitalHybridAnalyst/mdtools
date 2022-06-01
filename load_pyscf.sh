export PYTHONPATH=~/packages/pyscf_jojo:~/miniconda3/lib/python3.9/site-packages/:/home/chhli/share
module load gcc/9.2.0 mkl-2017.0.098
export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_def.so:$MKLROOT/lib/intel64/libmkl_sequential.so:$MKLROOT/lib/intel64/libmkl_core.so

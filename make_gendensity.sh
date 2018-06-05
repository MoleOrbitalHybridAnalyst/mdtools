#g++ gendensity.cpp -o gendensity -lpdb -std=c++0x -L/home/chhli/local/lib -I/home/chhli/share/pdb -Wall
#g++ gendensity.cpp -o gendensity -lpdb -L/home/chhli/local/lib -I/home/chhli/share/pdb/include -Wall
module load intelmpi/2017.up4+intel-17.0 gcc/6.1
mpicxx gendensity.cpp -o gendensity -lpdb -L/home/chhli/local/lib -I/home/chhli/share/pdb/include -Wall -O3
#mpicxx gendensity.cpp -o gendensity -lpdb -L/home/chhli/local/lib -I/home/chhli/share/pdb/include -Wall -g

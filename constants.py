# energy in kcal/mol
kcal = 1
kj = 0.23900573614
hartree = 627.509608030
eV = 23.060541945

# lenght in A
m = 1e10
A = 1
Bohr = 0.529177249

# time in ps
s = 1e12
ps = 1
fs = 1e-3
ns = 1e3
us = 1e6

NA = 6.0221409e+23                           # in mol^-1

# mass in kcal/mol /A^2 * ps^2
# 1 kg = 1e-3 kJ / (m^2/s^2) = 1e-3 kj/mol/NA / (m^2/s^2) * NA = 1e-3 kj/(m^2/s^2) * NA
kg = 1e-3 * kj / (m**2/s**2) * NA 

# physical constants
hbar = 0.01517872841498232                   # in kcal/mol * ps
kb = 0.001985875                             # in kcal/mol/K
c = 299792458 * m / s                        # in A / ps
mp = 1.67262192369E-27 * kg                  # in kcal/mol / A^2 * ps^2

# derived units
au_t = hbar / hartree  # https://manual.cp2k.org/trunk/CP2K_INPUT/MOTION/PRINT/VELOCITIES.html

# C++ program to build 2D heterostructures

(1) to compile 
g++ -o match.x main.cpp aims.cpp coincide.cpp

(20 to use 

./match.x $prefix1 $prefix2 $theta0 $theta1 $dtheta -$flag1 $value1 -$flag2 $value2 .....

Required input parameters:
"prefix1" and "prefix2": they are the names of input geometries. Geometries must be in FHI-AIMS format, and the filename should be "prefix1.in" and "prefix2.in"
The matched supercells are searched between angles "theta0" and "theta1", with "dtheta" as the step. They are in "degrees".

Optional input parameters:
Nmax: defaut 10 supercell can not be larger than "Nmax" * "Nmax"
dist: default 3, interlayer distance between the two structure
tolerance: default 0.02, maximum strain in lattice vectors.
angle_tolerance: default 0.05, maximum angle strain epsilon_21.
nat_max: default no limit, maximum number of atoms in the supercell 
n_sol_output: default no limit, number of configurations to output. (from small to large supercell size)
max_tilt: default 0, the angle between lattice vectors in supercell are constraint to be between "max_tilt" and "180-max_tilt". This value is to used to avoid very "thin" supercells.

When you match two hexagonal(square) lattices, you can always find hexagonal(square) solutions. If this cases, 'max_tilt" doesn't work. The output solutions are constraint to be hexagonal(square). In other cases, recommended value of "max_tilt" is 45 degree.


Example:
./match.x mos2 c 0 30 0.1 -Nmax 10 -tolerance 0.04 -n_sol_output 5

# C++ program to build 2D heterostructures

(1) to compile 
g++ -o match.x main.cpp aims.cpp coincide.cpp

(20 to use 

./match.x $prefix1 $prefix2 $theta0 $theta1 $dtheta -$flag1 $value1 -$flag2 $value2 .....

Required input parameters:
"prefix1" and "prefix2": they are the names of input geometries. Geometries must be in FHI-AIMS format, and the filename should be "prefix1.in" and "prefix2.in"
The matched supercells are searched between angles "theta0" and "theta1", with "dtheta" as the step.

Optional input parameters:
Nmax: defaut 10 supercell can not be larger than "Nmax" * "Nmax"
dist: default 3, interlayer distance between the two structure
tolerance: default 2%, maximum strain in lattice vectors.
angle_tolerance: default 0.05, maximum angle strain epsilon_21.
nat_max: default no limit, maximum number of atoms in the supercell 
n_sol_output: default no limit, number of configurations to output. (from small to large supercell size)
max_tilt: default 0, the angle between lattice vectors in supercell are constraint to be between "max_tilt" and "180-max_tilt". This value is to used to avoid very "thin" supercells.


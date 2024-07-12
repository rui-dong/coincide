#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdlib>
#include <math.h>
#include "aims.h"

#define PI 3.1415926
#define N_VECTOR_MAX 400000        // max number of matched vector pairs
#define N_SOL_FACTOR 1         // number of solution per rotation angle
#define BIG_NUM 65535
#define EPS 1e-4
#define DEBUG_MODE 1

int brav_lattice ( double , double , double , double , double ) ;

class coincide {
    //private
    int ** solution_tmp;
    double ** norm_tmp;
    int n_sol;
    int brav_1, brav_2; // 1 SQR 2 REC 3 HEX 4 CRT 5 OBL
    double cos_angle,angle;

    void eliminate_linear ( void );
    int find_min_area ( int );
    void sort_solution ( void );
    
  public:
    int fix_const; // FIX constraint 
	double angle_const;
    double angle_start, angle_end, angle_step;
    int num_angle,num_solution;
    double norm_a1, norm_a2, norm_b1, norm_b2, angle1, angle2;
    int N;
    int N_atom_limit;
    double tolerance, angle_tolerance, max_tilt;
    double a[2][2],b[2][2];
    int nat1,nat2;
    int ** solution;
    int ** solution_final;
    double ** strain;

    coincide ( aims, aims, double, double, double );
    void loop ( void );
    void solution_output ( char * );
    void strain_calculator ( int );

    FILE * debugFile, * outputFile;
 
};



#ifndef AIMS_H
#define AIMS_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <cstdlib>
#include <math.h>

#define NAT_MAX 5000
#define EPS 1e-4
#define PI 3.1415926
#define DEBUG_AIMS 0

using namespace std;

void alatt2cellp (double [3][3] , double [6] ) ;
void cellp2alatt (double [6] , double [3][3] ) ;
double vec_norm( double * , int ) ;
double vec_dot( double * , double * , int ) ;
void average_cell( double [3][3] , double [3][3] , double [3][3] ) ;

class aims {
    int i,j,k,l,m,n;
  public:
    int N, M, sys_type;
    char prefix [64];
    char prefix_out [64];

    int nat,nlatt;
    double xyz[NAT_MAX][3], xyz_frac[NAT_MAX][3];
    char name[NAT_MAX][10];
    int frac_flag[NAT_MAX]; // mixed "atom" and "atom_frac" are allowed
    double alatt[3][3];
    double cell_p[6];
  
    int if_frac; // 1: output fractional coordinates
    int fix_axis[3]; // 1: fix cell dimension in this direction
    
    aims (char *);
    aims ();
    virtual ~aims() {};
    int output (char *);
    void xyz2frac ( void );
    void frac2xyz ( void );
    void multiple( int, int, int, int );
    void simple_multiple( int, int, int);
    void supercell ( int [4] );
    void displace ( double [3] );
};


int aims_merge ( aims * , aims , aims , int , int , int) ; 
int aims_merge_pointer ( aims * , aims *, aims *, int , int , int) ; 





#endif

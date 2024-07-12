#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <cstdlib>
#include <math.h>

#include "aims.h"
#include "coincide.h"
/*
#define N_PER_FILE 10000
#define n_sol_output 5
*/
/*
	Quite good version for automatic generation of a small number of configurations
*/


main( int argc, char *argv[]){
  
  char * prefix;
  int i,j,k,if_angle=0;
  double cellp[6];
  double alatt[3][3]={{1, 2, 3},{4,5,6},{7,8,9}};
  double angle_start, angle_end, angle_step,angle;
  double tolerance, angle_tolerance, dist, max_tilt;
  int nat_max, fix_hex, N_max, n_sol_output;
  int m[4]={-1,1,2,1};
  char fileout[64];
  int count;
  double vector[3];

  // two candidate layers
  prefix=argv[1]; aims ref1 ( prefix );
  prefix=argv[2]; aims ref2 ( prefix );
  // rotation angle 
  angle_start=atof(argv[3]);
  angle_end=atof(argv[4]);
  angle_step=atof(argv[5]);
  // optional parameters
  // initial value
  N_max=0;
  dist=0;
  nat_max=0;
  fix_hex=-1;
  tolerance=0; angle_tolerance=0;
  max_tilt=0;
  n_sol_output=-1;
  i=6;
  while ( i<argc ) {
    if ( *argv[i]==45 ) {
      if ( strcmp(argv[i]+1,"Nmax")==0 ) {
        N_max=atoi(argv[++i]);
        printf("find input value \"Nmax=%d\".\n",N_max);
      } else if ( strcmp(argv[i]+1,"dist")==0 ) {
        dist=atof(argv[++i]);
		printf("find input value \"distance=%f\".\n",dist);
      } else if ( strcmp(argv[i]+1,"tolerance")==0 ) {
        tolerance=atof(argv[++i]);
		printf("find input value \"tolerance=%f\".\n",tolerance);
      } else if ( strcmp(argv[i]+1,"angle_tolerance")==0 ) {
        angle_tolerance=atof(argv[++i]);
		printf("find input value \"angle_tolerance=%f\".\n",angle_tolerance);
      } else if ( strcmp(argv[i]+1,"nat_max")==0 ) {
        nat_max=atoi(argv[++i]);
		printf("find input value \"nat_max=%d\".\n",nat_max);
       } else if ( strcmp(argv[i]+1,"n_sol_output")==0 ) {
        n_sol_output=atoi(argv[++i]);
		printf("find input value \"n_sol_output=%d\".\n",n_sol_output);
      //} else if ( strcmp(argv[i]+1,"fix_hex")==0 ) {
      //  fix_hex=atof(argv[++i]);
      //  printf("find input value \"fix_hex=%d\".\n",fix_hex);
      } else if ( strcmp(argv[i]+1,"max_tilt")==0 ) {
        max_tilt=atof(argv[++i]);
        printf("find input value \"max_tilt=%f\".\n",max_tilt);
      } else {
        printf(" input arguments are in wrong format!!\n");
        abort();
      }
    } else {
      printf(" input arguments are in wrong format!!\n");
      abort();
    }
    i++;
  }


  
  // inter-layer distances, lower "1" and lift "2"
  double z_min,z_max,z_tot;
  double com1[3],com2[3];
  for ( i=0;i<3;i++ ) { com1[i]=0; com2[i]=0; } 
  for ( j=0;j<3;j++ ) { // find the "center of mass" of each layer
    for ( i=0;i<ref1.nat;i++ ) { com1[j] += ref1.xyz[i][2]; }
    com1[j] /= ref1.nat;
    for ( i=0;i<ref2.nat;i++ ) { com2[j] += ref2.xyz[i][2]; }
    com2[j] /= ref2.nat;
  }
  // shift
  vector[0]=0; vector[1]=0;
  if ( dist!=0 ) { // pre-set distance between the out-most atoms
    z_max=-5000; z_min=5000;
    for ( i=0;i<ref1.nat;i++ ) { 
      if ( z_max < ref1.xyz[i][2]) { z_max = ref1.xyz[i][2]; }
    }
    vector[2]=-z_max-dist/2;
    ref1.displace(vector);
    for ( i=0;i<ref2.nat;i++ ) { 
      if ( z_min > ref2.xyz[i][2]) { z_min = ref2.xyz[i][2]; }
    }
    vector[2]=-z_min+dist/2;
    ref2.displace(vector);
    // find the "thickness" of junction, to calculate supercell size along z
    z_max=-5000; z_min=5000;
    for ( i=0;i<ref2.nat;i++ ) { 
      if ( z_max < ref2.xyz[i][2]) { z_max = ref2.xyz[i][2]; }
    }
    for ( i=0;i<ref1.nat;i++ ) { 
      if ( z_min > ref1.xyz[i][2]) { z_min = ref1.xyz[i][2]; }
    }
    // set the lattice constant and re-calculate the fractional coord
	z_tot=z_max-z_min+15 ;
	if (ref1.alatt[2][2] < z_tot) { ref1.alatt[2][2]=z_tot; } 
	if (ref2.alatt[2][2] < z_tot) { ref2.alatt[2][2]=z_tot; } 
    ref1.xyz2frac(); ref2.xyz2frac();

  } else { // move the layer to the "middle "
    vector[2]=-com1[2]-ref1.alatt[2][2]/2;
    ref1.displace(vector);
    vector[2]=-com2[2]+ref2.alatt[2][2]/2;
    ref2.displace(vector);
    // z of the supercell = sum of the component
    ref1.alatt[2][2]+=ref2.alatt[2][2]; ref1.xyz2frac();
    ref2.alatt[2][2]=ref1.alatt[2][2]; ref2.xyz2frac();
  }

  //ref1.output("tmp.1");
  //ref2.output("tmp.2");


  // create the "working" layers
  aims layer0, layer1, layer2;
  layer1=ref1;
  layer2=ref2;

  // create the job of lattice matching
  printf("==== creating lattice matching job of %s and %s\n",argv[1],argv[2]);
  printf("==== rotation angle from %f to %f, with stepsize=%f\n",angle_start, angle_end, angle_step);
  coincide match (layer1, layer2, angle_start, angle_end, angle_step);
  // overwrite the parameters (if input was found)
  if ( N_max!=0 ) {  match.N=N_max; }
  if ( nat_max!=0 ) { match.N_atom_limit=nat_max; }
  if ( fix_hex!=-1 ) { match.fix_const=fix_hex; }
  if ( tolerance!=0 ) { match.tolerance=tolerance; }
  if ( angle_tolerance!=0 ) { match.angle_tolerance=angle_tolerance; }
  if ( max_tilt!=0 ) { match.max_tilt=max_tilt; }
  // echoing the parameters 
  printf("==== parameters for lattice matching ====\n");
  printf("== Maximum lattice size for searching: %d\n",match.N);
  printf("== Maximum total number of atom: %d\n",match.N_atom_limit);
  printf("== Tolerance for vector pairs: %f\n",match.tolerance);
  printf("== Angle tolerance for vector pairs: %f\n",match.angle_tolerance);
  printf("== Level of fixing hexagonal cell: %d\n",match.fix_const);
  
  //
  printf("==== starting the loop ====\n");
  match.loop();
  printf("==== loop finished ====\n");

  sprintf(fileout,"solution_final.log");
  match.solution_output(fileout);

  if ( n_sol_output == -1 ) {
    n_sol_output = match.num_solution;
    printf("== output all configurations generated by the code: %d in total\n",n_sol_output);
  } else {
    printf("== output the first %d configurations generated by the code.\n",n_sol_output);
  }
  
  for ( i=0;i<n_sol_output;i++ ) {
    if ( match.solution_final[i][0] == 0 && match.solution_final[i][1] == 0 ) {
      continue;
    }
    printf("processing # %d at i_angle=%d : ",i,match.solution_final[i][8]);

    layer1=ref1;
    layer2=ref2;

    //get surface 1
    printf("== layer 1: ");
    for ( j=0;j<4;j++ ) {
      m[j]=match.solution_final[i][j];
      printf("%d ",m[j]);
    }
    layer1.supercell(m);

    //get surface 2
    printf("== layer 2: ");
    for ( j=0;j<4;j++ ) {
      m[j]=match.solution_final[i][j+4];
      printf("%d ",m[j]);
    }
    printf("==\n");
    layer2.supercell(m);

    aims_merge(&layer0,layer1,layer2,1,1,0); //frac-frac-avg_cell
    angle=match.angle_step*match.solution_final[i][8]+match.angle_start;

    sprintf(fileout,"config_%d.in",i+1);
    //sprintf(fileout,"config_%d.in",argv[1],argv[2],i+1);
    layer0.output(fileout);
  } 


}










